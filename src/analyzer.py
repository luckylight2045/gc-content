# src/analyzer.py

import re

def find_motif(seq: str, motif: str) -> list[int]:
    """
    Return list of positions (1-based) where motif occurs in sequence.
    """
    seq = seq.upper()
    motif = motif.upper()
    return [m.start()+1 for m in re.finditer(motif, seq)]

def calculate_gc(seq: str) -> float:
    """Return GC% (rounded to 2 decimals) for sequence (letters only)."""
    seq = seq.upper()
    # keep only A,T,G,C
    bases = [b for b in seq if b in ("A", "T", "G", "C")]
    total = len(bases)
    if total == 0:
        return 0.0
    g = bases.count("G")
    c = bases.count("C")
    return round(((g + c) / total) * 100, 2)

def gc_skew(seq: str, window: int = 100, step: int = 10) -> list[tuple[int,float]]:
    """
    Return GC skew (G-C)/(G+C) in sliding windows.
    Returns list of tuples: (start position, skew value)
    """
    seq = seq.upper()
    seq = "".join([b for b in seq if b in ("A","T","G","C")])
    result = []
    for start in range(0, max(1, len(seq)-window+1), step):
        w = seq[start:start+window]
        g = w.count("G")
        c = w.count("C")
        skew = (g - c)/(g + c) if (g + c) != 0 else 0
        result.append((start, skew))
    return result


def find_orfs(seq: str) -> list[dict]:
    """
    Find ORFs in a DNA sequence.
    Returns a list of dictionaries: {'start': int, 'end': int, 'length': int, 'sequence': str}
    """
    seq = seq.upper()
    stop_codons = {"TAA", "TAG", "TGA"}
    orfs = []
    for frame in range(3):
        i = frame
        while i < len(seq) - 2:
            codon = seq[i:i+3]
            if codon == "ATG":
                start = i
                i += 3
                while i < len(seq) - 2 and seq[i:i+3] not in stop_codons:
                    i += 3
                end = i + 3 if i < len(seq) - 2 else len(seq)
                orfs.append({"start": start+1, "end": end, "length": end-start, "sequence": seq[start:end]})
            else:
                i += 3
    return orfs


def translate_dna(seq: str) -> str:
    """Translate DNA to protein (single ORF, no introns)."""
    codon_table = {
        # Standard genetic code
        'ATA':'I','ATC':'I','ATT':'I','ATG':'M',
        'ACA':'T','ACC':'T','ACG':'T','ACT':'T',
        'AAC':'N','AAT':'N','AAA':'K','AAG':'K',
        'AGC':'S','AGT':'S','AGA':'R','AGG':'R',
        'CTA':'L','CTC':'L','CTG':'L','CTT':'L',
        'CCA':'P','CCC':'P','CCG':'P','CCT':'P',
        'CAC':'H','CAT':'H','CAA':'Q','CAG':'Q',
        'CGA':'R','CGC':'R','CGG':'R','CGT':'R',
        'GTA':'V','GTC':'V','GTG':'V','GTT':'V',
        'GCA':'A','GCC':'A','GCG':'A','GCT':'A',
        'GAC':'D','GAT':'D','GAA':'E','GAG':'E',
        'GGA':'G','GGC':'G','GGG':'G','GGT':'G',
        'TCA':'S','TCC':'S','TCG':'S','TCT':'S',
        'TTC':'F','TTT':'F','TTA':'L','TTG':'L',
        'TAC':'Y','TAT':'Y','TAA':'*','TAG':'*',
        'TGC':'C','TGT':'C','TGA':'*','TGG':'W',
    }
    seq = seq.upper()
    protein = ""
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        protein += codon_table.get(codon, "X")
    return protein


def validate_sequence(seq: str) -> bool:
    """Return True if sequence contains only A, T, G, C, N (optional), else False."""
    allowed = {"A", "T", "G", "C", "N"}
    return all(base.upper() in allowed for base in seq if base.isalpha())


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    seq = seq.upper()
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(complement.get(base, base) for base in reversed(seq))


def sliding_window_gc(seq: str, window: int = 100, step: int = 10):
    """
    Return list of (position, gc%) for sliding window.
    position = start index (0-based) of the window.
    """
    seq = seq.upper()
    # keep only A,T,G,C
    seq = ''.join([b for b in seq if b in ("A", "T", "G", "C")])
    results = []
    for start in range(0, max(1, len(seq) - window + 1), step):
        window_seq = seq[start:start + window]
        gc = calculate_gc(window_seq)
        results.append((start, gc))
    return results

def nucleotide_counts(seq: str) -> dict:
    """Return a dictionary with counts of A, T, G, C and total length."""
    seq = seq.upper()
    bases = [b for b in seq if b in ("A", "T", "G", "C")]
    return {
        "A": bases.count("A"),
        "T": bases.count("T"),
        "G": bases.count("G"),
        "C": bases.count("C"),
        "Length": len(bases)
    }

