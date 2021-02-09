"""Helper functions for computing features of gRNAs"""
import numpy as np


def compute_gc_content(seq):
    sum_gc = seq.count('G') + seq.count('C')
    return float(sum_gc)/len(seq)


def reverse_compliment(seq):
    tab = str.maketrans("ACGT", "TGCA")
    return seq.translate(tab)[::-1]


def detect_palindrome(seq):
    rc = reverse_compliment(seq)
    return seq == rc


def compute_weight(seq):
    from Bio import SeqUtils
    """Compute the Molecular Weight of the sgRNAs"""
    rna = seq.replace('T', 'U')
    return int(SeqUtils.molecular_weight(rna, seq_type='RNA'))


def compute_entropy(seq, k=1):
    """Compute the Shannon Entropy of the spacer
     for all subsequences of length k (k=1 by default)
     """
    assert k < 20
    total = len(seq)/k
    subseqs = set([seq[i:i+k] for i in range(len(seq)-k+1)])
    entropy = -1 * sum([(seq.count(kmer)/total) * np.log(seq.count(kmer)/total) for kmer in subseqs])
    return entropy/np.log(len(subseqs))


def detect_start(seq):
    return "ATG" in seq


def compute_fold_change(t1, t2):
    return -np.log2(t2/t1)
