from datetime import datetime
from file_handler import reverse_complement
from endonucleases import restriction_enzymes
import re

degenerate_map = {
    'A': 'A',
    'T': 'T',
    'C': 'C',
    'G': 'G',
    'N': '[ATGC]',
    'W': '[AT]',
    'R': '[AG]',
    'Y': '[CT]',
    'S': '[GC]',
    'K': '[GT]',
    'M': '[AC]',
    'B': '[CGT]',
    'D': '[AGT]',
    'H': '[ACT]',
    'V': '[ACG]'
}

def validate_enzyme_sequence(sequence):
    valid_bases = set('ATCGNWRYSKMBDHV')
    invalid_chars = set(char.upper() for char in sequence) - valid_bases
    if invalid_chars:
        raise ValueError(f"Invalid characters found in sequence: {invalid_chars}")
    return True

def reverse_complement_enzyme_sequence(sequence):
    complement_dict = {
        'A': 'T', 'T': 'A',
        'C': 'G', 'G': 'C',
        'N': 'N',
        'W': 'W', 
        'R': 'Y', 'Y': 'R',
        'S': 'S',
        'K': 'M', 'M': 'K',
        'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B'
    }
    return ''.join(complement_dict.get(base.upper(), base) for base in reversed(sequence))

def convert_to_regex(sequence):
    regex = "".join(degenerate_map.get(base, f"[{base}]") for base in sequence)
    return regex

def recognition_sequence_blunt(sequence, enzyme):
    recognition_sequence = enzyme["sequence"]
    cut_position = enzyme["cut"]
    
    regex = convert_to_regex(recognition_sequence)
    matches = [(m.start(), m.start() + cut_position) for m in re.finditer(regex, sequence)]     # re.finditer vrací nepřekrývající se shody
    cut_sites = [match[1] for match in matches]
    return cut_sites

def cut_blunt(sequence, enzyme_name):
    enzyme = restriction_enzymes[enzyme_name]
    cut_sites = recognition_sequence_blunt(sequence, enzyme)
    
    fragments = []
    start = 0
    for cut_site in cut_sites:
        fragments.append(sequence[start:cut_site])
        start = cut_site
    fragments.append(sequence[start:])  # fragment za posledním řezem
    
    return fragments

def complement_bases_enzyme(seq):
    comp = {
        'A': 'T', 'T': 'A',
        'C': 'G', 'G': 'C',
        'N': 'N',
        'W': 'W', 'R': 'Y', 'Y': 'R',
        'S': 'S',
        'K': 'M', 'M': 'K',
        'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B'
    }
    return "".join(comp.get(base.upper(), base) for base in seq)

def recognition_sequence_sticky(sequence, enzyme):
    # Forward strand:
    recognition_seq = enzyme["sequence"]
    forward_regex = convert_to_regex(recognition_seq)
    forward_matches = list(re.finditer(forward_regex, sequence))
    forward_cut_sites = [m.start() + enzyme["cut"] for m in forward_matches]
    forward_cut_sites = sorted(set(forward_cut_sites))
    # Reverse strand:
    rev_recognition_seq = reverse_complement_enzyme_sequence(enzyme["sequence"])
    reverse_regex = convert_to_regex(rev_recognition_seq)
    rev_sequence = reverse_complement(sequence)
    reverse_matches = list(re.finditer(reverse_regex, rev_sequence))
    reverse_cut_sites = [m.start() + enzyme["cut"] for m in reverse_matches]
    reverse_cut_sites = sorted(set(reverse_cut_sites))
    
    return forward_cut_sites, reverse_cut_sites

def cut_sticky(sequence, enzyme_name):
    enzyme = restriction_enzymes[enzyme_name]
    forward_sites, reverse_sites = recognition_sequence_sticky(sequence, enzyme) # Získáme cut pozice pro forward a reverse strand samostatně
    # Cut forward
    forward_fragments = []
    start = 0
    for site in forward_sites:
        forward_fragments.append(sequence[start:site])
        start = site
    forward_fragments.append(sequence[start:])
    # Cut reverse
    rev_seq = reverse_complement(sequence)
    rev_fragments = []
    start = 0
    for site in reverse_sites:
        rev_fragments.append(rev_seq[start:site])
        start = site
    rev_fragments.append(rev_seq[start:])
    
    return forward_fragments, rev_fragments

def process_sequence(sequence, enzymes):
    results = {}
    for enzyme_name in enzymes:
        sites = recognition_sequence_blunt(sequence, restriction_enzymes[enzyme_name])
        results[enzyme_name] = sites  # enzyme is a string
    return results