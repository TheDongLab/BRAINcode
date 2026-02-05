#!/usr/bin/env python3
# ~/donglab/pipelines/scripts/rnaseq/snap_junctions_to_exons.py

import sys
from bisect import bisect_left

if len(sys.argv) != 4:
    print("Usage: python snap_junctions_to_exons.py <refflat> <input_junction> <output_junction>")
    sys.exit(1)

refflat_file = sys.argv[1]
input_junction = sys.argv[2]
output_junction = sys.argv[3]

def get_closest(sorted_list, val):
    """Find closest value in sorted list using binary search"""
    pos = bisect_left(sorted_list, val)
    if pos == 0: 
        return sorted_list[0]
    if pos == len(sorted_list): 
        return sorted_list[-1]
    before, after = sorted_list[pos-1], sorted_list[pos]
    return before if after - val > val - before else after

# Build reference exon boundary dictionary
ref_bounds = {}
with open(refflat_file, 'r') as f:
    for line in f:
        cols = line.strip().split('\t')
        if len(cols) < 11: 
            continue
        chrom = cols[2]
        if chrom not in ref_bounds: 
            ref_bounds[chrom] = set()
        
        # Add all exon starts
        for s in cols[9].split(','): 
            if s: 
                ref_bounds[chrom].add(int(s))
        
        # Add all exon ends
        for e in cols[10].split(','): 
            if e: 
                ref_bounds[chrom].add(int(e))

# Sort sets for binary search efficiency
sorted_ref = {k: sorted(list(v)) for k, v in ref_bounds.items()}

print(f"[INFO] Loaded exon boundaries for {len(sorted_ref)} chromosomes")

# Snap junctions to nearest exon boundaries
snapped_count = 0
total_count = 0

with open(input_junction, 'r') as fin, open(output_junction, 'w') as fout:
    for line in fin:
        cols = line.strip().split('\t')
        chrom, s, e = cols[0], int(cols[1]), int(cols[2])
        total_count += 1
        
        if chrom in sorted_ref:
            # Snap start/end to nearest known exon boundary within 500bp
            snap_s = get_closest(sorted_ref[chrom], s)
            snap_e = get_closest(sorted_ref[chrom], e)
            
            if abs(snap_s - s) < 500:
                s = snap_s
                snapped_count += 1
            if abs(snap_e - e) < 500:
                e = snap_e
                snapped_count += 1
        
        cols[1], cols[2] = str(s), str(e)
        fout.write("\t".join(cols) + "\n")

print(f"[INFO] Snapped {snapped_count} boundaries out of {total_count * 2} total ({snapped_count / (total_count * 2) * 100:.1f}%)")
print(f"[DONE] Written snapped junctions to {output_junction}")
