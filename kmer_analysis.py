#!/usr/bin/env python3
from Bio import SeqIO
from collections import Counter

# To run :
# chmod +x term_analysis.py
# ./kmer_analysis.py /path/to/assembled_genome.fasta

def kmer_terminal_analysis(fasta_file, k=15, terminal_size=10000):
    """Compare k-mer frequencies at terminals"""

    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        seq_len = len(seq)

        print(f"\n=== K-mer Analysis: {record.id} ===")

        # Get terminal regions
        start_region = seq[:terminal_size]
        end_region = seq[-terminal_size:]

        # Count k-mers
        def get_kmers(sequence, k):
            kmers = []
            for i in range(len(sequence) - k + 1):
                kmers.append(sequence[i:i+k])
            return Counter(kmers)

        start_kmers = get_kmers(start_region, k)
        end_kmers = get_kmers(end_region, k)

        # Find shared k-mers
        shared = set(start_kmers.keys()) & set(end_kmers.keys())

        # Calculate similarity
        total_start = len(start_kmers)
        total_end = len(end_kmers)
        shared_count = len(shared)

        similarity = shared_count / min(total_start, total_end) * 100

        print(f"Terminal region size: {terminal_size:,} bp")
        print(f"K-mer size: {k}")
        print(f"Start unique k-mers: {total_start:,}")
        print(f"End unique k-mers: {total_end:,}")
        print(f"Shared k-mers: {shared_count:,}")
        print(f"Similarity: {similarity:.2f}%")

        # Check for exact overlaps
        overlap_found = False
        for kmer in shared:
            if start_kmers[kmer] > 5 and end_kmers[kmer] > 5:
                # High frequency in both - possible overlap
                overlap_found = True
                break

        if similarity > 70 or overlap_found:
            print("✓ High terminal similarity suggests CIRCULAR")
        elif similarity < 30:
            print("✗ Low terminal similarity suggests LINEAR")
        else:
            print("? Moderate similarity - inconclusive")

if __name__ == "__main__":
    import sys
    kmer_terminal_analysis(sys.argv[1])
