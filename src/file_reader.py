# src/file_reader.py
from Bio import SeqIO

def read_fasta(path: str) -> str:
    """Read first sequence from FASTA (or concat multiple) and return sequence string."""
    seqs = []
    for record in SeqIO.parse(path, "fasta"):
        seqs.append(str(record.seq))
    return ''.join(seqs)
