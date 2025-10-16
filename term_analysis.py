#!/usr/bin/env python3
from Bio import SeqIO
from collections import Counter
import numpy as np

# To run :
# chmod +x term_analysis.py
# ./term_analysis.py /path/to/assembled_genome.fasta

def analyze_terminal_composition(fasta_file, window=1000, num_windows=10):
    """Compare composition across genome including terminals"""
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        seq_len = len(seq)
        
        print(f"\n=== Compositional Analysis: {record.id} ===")
        
        # Calculate GC content in windows
        gc_contents = []
        positions = []
        
        # Start region
        for i in range(num_windows):
            window_seq = seq[i*window:(i+1)*window]
            gc = (window_seq.count('G') + window_seq.count('C')) / len(window_seq) * 100
            gc_contents.append(gc)
            positions.append(i*window)
        
        # Middle region
        mid = seq_len // 2
        for i in range(-num_windows//2, num_windows//2):
            start = mid + i*window
            window_seq = seq[start:start+window]
            gc = (window_seq.count('G') + window_seq.count('C')) / len(window_seq) * 100
            gc_contents.append(gc)
            positions.append(start)
        
        # End region
        for i in range(num_windows):
            start = seq_len - (num_windows-i)*window
            window_seq = seq[start:start+window]
            gc = (window_seq.count('G') + window_seq.count('C')) / len(window_seq) * 100
            gc_contents.append(gc)
            positions.append(start)
        
        # Calculate variance
        gc_variance = np.var(gc_contents)
        gc_mean = np.mean(gc_contents)
        
        print(f"Average GC: {gc_mean:.2f}%")
        print(f"GC variance: {gc_variance:.4f}")
        
        # Compare terminal regions
        start_gc = np.mean(gc_contents[:num_windows])
        end_gc = np.mean(gc_contents[-num_windows:])
        terminal_diff = abs(start_gc - end_gc)
        
        print(f"Start GC: {start_gc:.2f}%")
        print(f"End GC: {end_gc:.2f}%")
        print(f"Terminal difference: {terminal_diff:.2f}%")
        
        # Interpretation
        if terminal_diff < 1.0 and gc_variance < 2.0:
            print("✓ Uniform composition suggests CIRCULAR genome")
        elif terminal_diff > 3.0:
            print("✗ Terminal composition differs - likely LINEAR or degraded")
        else:
            print("? Inconclusive - moderate variation")

if __name__ == "__main__":
    import sys
    analyze_terminal_composition(sys.argv[1])