import os
import requests
import pysam
from Bio import Entrez
from endonucleases import restriction_enzymes
from logic import cut_blunt, cut_sticky
from file_handler import reverse_complement

# TEMP_DIR relativní ke cwd
TEMP_DIR = "temp_files"
os.makedirs(TEMP_DIR, exist_ok=True)

# --- EVA metadata ---
def get_eva_species_entry(organism: str, assembly_name: str) -> dict:
    url = "https://www.ebi.ac.uk/eva/webservices/rest/v1/meta/species/list?loaded=true"
    r = requests.get(url)
    r.raise_for_status()
    entries = []
    for batch in r.json().get('response', []):
        entries.extend(batch.get('result', []))
    for e in entries:
        if organism in (e.get('taxonomyScientificName'),
                        e.get('taxonomyCommonName'),
                        e.get('taxonomyCode'),
                        e.get('taxonomyEvaName')):
            if assembly_name in (e.get('assemblyName'),
                                 e.get('assemblyCode'),
                                 e.get('assemblyAccession')):
                return e
    raise ValueError(f"Species/assembly not found: {organism}/{assembly_name}")

# --- Download VCF & CSI ---
def download_vcf_and_csi(entry: dict) -> tuple[str, str]:
    sci = entry['taxonomyScientificName']
    species_dir = sci.lower().replace(' ', '_')
    asm_dir = entry['assemblyName'].replace('-', '').replace(' ', '')
    tax_id = entry['taxonomyId']
    acc = entry['assemblyAccession']
    vcf_fn = f"{tax_id}_{acc}_current_ids.vcf.gz"
    csi_fn = vcf_fn + ".csi"
    base_url = (
        f"https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_7/"
        f"by_species/{species_dir}/{asm_dir}/"
    )

    vcf_path = os.path.join(TEMP_DIR, vcf_fn)
    csi_path = os.path.join(TEMP_DIR, csi_fn)

    def dl(url, path):
        r = requests.get(url, stream=True)
        r.raise_for_status()
        if 'text/html' in r.headers.get('Content-Type', ''):
            raise ValueError(f"Invalid download: {url}")
        with open(path, 'wb') as f:
            for chunk in r.iter_content(8192):
                f.write(chunk)

    if not os.path.exists(vcf_path):
        dl(base_url + vcf_fn, vcf_path)
    if not os.path.exists(csi_path):
        dl(base_url + csi_fn, csi_path)

    return vcf_path, csi_path

# --- Parse VCF by CSI ---
def load_variations_from_vcf_gz(vcf_gz: str, chrom: str, start: int, end: int) -> list:
    csi = vcf_gz + ".csi"
    tb = pysam.TabixFile(vcf_gz, index=csi)
    if chrom not in tb.contigs and f"chr{chrom}" in tb.contigs:
        chrom = f"chr{chrom}"
    vars = []
    for rec in tb.fetch(chrom, start, end):
        cols = rec.split('\t')
        pos, ref, alt = int(cols[1]), cols[3], cols[4]
        if len(ref) == 1 and all(len(a) == 1 for a in alt.split(',')):
            for a in alt.split(','):
                vars.append({'position': pos, 'alleles': f"{ref}/{a}"})
    return vars

# --- IUPAC mapping ---
iupac_map = {
    frozenset(['A']): 'A', frozenset(['C']): 'C',
    frozenset(['G']): 'G', frozenset(['T']): 'T',
    frozenset(['A','G']): 'R', frozenset(['C','T']): 'Y',
    frozenset(['C','G']): 'S', frozenset(['A','T']): 'W',
    frozenset(['G','T']): 'K', frozenset(['A','C']): 'M',
    frozenset(['C','G','T']): 'B', frozenset(['A','G','T']): 'D',
    frozenset(['A','C','T']): 'H', frozenset(['A','C','G']): 'V',
    frozenset(['A','C','G','T']): 'N'
}

def convert_variations_to_iupac(seq: str, vars: list, region_start: int) -> str:
    s = list(seq)
    for v in vars:
        idx = v['position'] - region_start
        if 0 <= idx < len(s):
            alleles = v['alleles'].split('/')
            key = frozenset(a.upper() for a in alleles)
            s[idx] = iupac_map.get(key, 'N')
    return ''.join(s)

# --- FASTA fetch (assemblyAccession fallback logic) ---
def fetch_fasta_sequence(organism: str,
                         chromosome: str,
                         start: int,
                         end: int,
                         assembly: str = None) -> str:
    Entrez.email = Entrez.email or "your.email@domain.com"
    def esearch_term(term):
        h = Entrez.esearch(db="nuccore", term=term, retmax=1)
        ids = Entrez.read(h).get("IdList", [])
        h.close()
        return ids

    queries = []
    if assembly:
        queries.extend([
            f'"{assembly}"[Assembly Accession] AND chromosome {chromosome}[Title] AND RefSeq[Filter]',
            f'"{assembly}"[Assembly Accession] AND chr{chromosome}[Title] AND RefSeq[Filter]',
            f'"{assembly}"[Assembly Accession] AND chromosome_{chromosome}[Title] AND RefSeq[Filter]',
            f'"{assembly}"[Assembly Accession] AND chr_{chromosome}[Title] AND RefSeq[Filter]',
            f'"{organism}"[Organism] AND "{assembly}"[Assembly Accession] AND RefSeq[Filter]',
        ])
    queries.extend([
        f'"{organism}"[Organism] AND chromosome {chromosome}[Title] AND RefSeq[Filter]',
        f'"{organism}"[Organism] AND chr{chromosome}[Title] AND RefSeq[Filter]',
        f'"{organism}"[Organism] AND chromosome_{chromosome}[Title] AND RefSeq[Filter]',
        f'"{organism}"[Organism] AND chr_{chromosome}[Title] AND RefSeq[Filter]',
    ])

    for term in queries:
        ids = esearch_term(term)
        if ids:
            ef = Entrez.efetch(db="nuccore",
                               id=ids[0],
                               rettype="fasta",
                               retmode="text",
                               seq_start=start,
                               seq_stop=end)
            lines = ef.read().splitlines()
            ef.close()
            return "".join(lines[1:])
    raise ValueError(f"NCBI RefSeq accession not found for {organism} chr{chromosome}"
                     + (f", assembly={assembly}" if assembly else ""))

def analyze_variant_enzyme_cutting(variant_sequences):
    results = {}
    for enzyme, enzyme_info in restriction_enzymes.items():
        enzyme_results = {}
        valid_cut = False
        for seq_id, seq in variant_sequences.items():
            try:
                overhang = enzyme_info.get('overhang', 0)
                if overhang > 0:
                    ffrags, rfrags = cut_sticky(seq, enzyme)
                    frag_data = {
                        "forward": {
                            "fragment_count": len(ffrags),
                            "fragments": ffrags,
                            "fragment_lengths": [len(x) for x in ffrags]
                        },
                        "reverse": {
                            "fragment_count": len(rfrags),
                            "fragments": rfrags,
                            "fragment_lengths": [len(x) for x in rfrags]
                        }
                    }
                    if len(ffrags) > 1 or len(rfrags) > 1:
                        valid_cut = True
                else:
                    bfrags = cut_blunt(seq, enzyme)
                    frag_data = {
                        "forward": {
                            "fragment_count": len(bfrags),
                            "fragments": bfrags,
                            "fragment_lengths": [len(x) for x in bfrags]
                        }
                    }
                    if len(bfrags) > 1:
                        valid_cut = True
                enzyme_results[seq_id] = frag_data
            except Exception as e:
                enzyme_results[seq_id] = {"error": str(e)}
        if valid_cut:
            results[enzyme] = enzyme_results

    enzyme_candidates = {}
    for enzyme, en_res in results.items():
        candidate_list = []
        for strand in ["forward", "reverse"]:
            counts = []
            for seq_id, data in en_res.items():
                if strand in data and "fragment_count" in data[strand]:
                    counts.append(data[strand]["fragment_count"])
            if counts and (max(counts) - min(counts) == 1):
                candidate_list.append(min(counts))
        if candidate_list:
            enzyme_candidates[enzyme] = min(candidate_list)

    if enzyme_candidates:
        min_candidate = min(enzyme_candidates.values())
        optimal_enzymes = []
        for enzyme, candidate in enzyme_candidates.items():
            if candidate == min_candidate:
                optimal_enzymes.append((enzyme, (candidate, candidate+1), results[enzyme]))
    else:
        min_candidate = None
        optimal_enzymes = []

    return results, optimal_enzymes, min_candidate


# --- Main ---
def run_variant_analysis(email: str, organism: str, assembly: str,
                         chromosome: str, region: str, snp_pos: int):
    Entrez.email = email
    start, end = map(int, region.split('-'))

    # 1) FASTA (cache)
    fasta_fn = f"ref_{organism.replace(' ', '_')}_{chromosome}_{start}_{end}.fasta"
    fasta_path = os.path.join(TEMP_DIR, fasta_fn)
    if not os.path.exists(fasta_path):
        seq = fetch_fasta_sequence(organism, chromosome, start, end, assembly)
        with open(fasta_path, 'w') as fw:
            fw.write('>ref\n' + seq)
    else:
        with open(fasta_path) as f:
            seq = ''.join(l.strip() for l in f if not l.startswith('>'))

    # 2) VCF
    entry = get_eva_species_entry(organism, assembly)
    vcfgz, _ = download_vcf_and_csi(entry)
    variants = load_variations_from_vcf_gz(vcfgz, chromosome, start, end)

    # 3) temp-VCF
    temp_vcf = os.path.join(TEMP_DIR, f"region_{chromosome}_{start}_{end}.vcf")
    with open(temp_vcf, 'w') as fw:
        fw.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for v in variants:
            r, a = v['alleles'].split('/')[:2]
            fw.write(f"{chromosome}\t{v['position']}\t.\t{r}\t{a}\t.\tPASS\t.\n")

    # 4) IUPAC
    iupac_seq = convert_variations_to_iupac(seq, variants, start)

    # 5) ref/alt
    match = next((v for v in variants if v['position']==snp_pos), None)
    if match:
        r, a = match['alleles'].split('/')[:2]
    else:
        r=a='N'
    base = list(iupac_seq)
    idx = snp_pos - start
    ref_seq = base[:] 
    alt_seq = base[:]
    ref_seq[idx]=r; alt_seq[idx]=a
    var_seqs = {'reference': ''.join(ref_seq),
                'alternative': ''.join(alt_seq)}

    # 6) enzyme
    all_res, opt, candidate = analyze_variant_enzyme_cutting(var_seqs)

    return {
        'sequence': seq,
        'variants': variants,
        'iupac_sequence': iupac_seq,
        'variant_sequences': var_seqs,
        'enzyme_results': all_res,
        'optimal': opt,
        'candidate': candidate,
        'fasta_file': fasta_path,
        'temp_vcf': temp_vcf
    }
