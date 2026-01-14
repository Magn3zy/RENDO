import os
from endonucleases import restriction_enzymes
from logic import cut_sticky, cut_blunt

def analyze_variant_enzyme_cutting(variant_sequences):
    results = {}
    
    for enzyme, enzyme_info in restriction_enzymes.items():
        enzyme_results = {}
        valid_cut = False  # nastaveno na True, pokud alespoň u jedné sekvence dojde k dělení (více než 1 fragment)
        for seq_id, seq in variant_sequences.items():
            overhang = enzyme_info.get('overhang', 0)
            try:
                if overhang > 0:
                    
                    forward_fragments, reverse_fragments = cut_sticky(seq, enzyme)
                    frag_data = {
                        "forward": {
                            "fragment_count": len(forward_fragments),
                            "fragments": forward_fragments,
                            "fragment_lengths": [len(frag) for frag in forward_fragments]
                        },
                        "reverse": {
                            "fragment_count": len(reverse_fragments),
                            "fragments": reverse_fragments,
                            "fragment_lengths": [len(frag) for frag in reverse_fragments]
                        }
                    }
                    enzyme_results[seq_id] = frag_data
                    if len(forward_fragments) > 1 or len(reverse_fragments) > 1:
                        valid_cut = True
                else:
                    
                    fragments = cut_blunt(seq, enzyme)
                    frag_data = {
                        "forward": {
                            "fragment_count": len(fragments),
                            "fragments": fragments,
                            "fragment_lengths": [len(frag) for frag in fragments]
                        }
                    }
                    enzyme_results[seq_id] = frag_data
                    if len(fragments) > 1:
                        valid_cut = True
            except Exception as e:
                enzyme_results[seq_id] = {"error": str(e)}
        if valid_cut:
            results[enzyme] = enzyme_results

    
    enzyme_candidates = {}  
    for enzyme, enzyme_results in results.items():
        candidate_list = []
        for strand in ["forward", "reverse"]:
            counts = []
            for seq_id, res in enzyme_results.items():
                if strand in res and "fragment_count" in res[strand]:
                    counts.append(res[strand]["fragment_count"])
            if counts:
                if max(counts) - min(counts) == 1:
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
