# -*- coding: utf-8 -*-

from typing import List


def read_seq_from_file (path_to_file: str) -> str:
    """Reads a sequence from a multi-line text file

    Args:
        path_to_file (str): Path to the file containing the sequence.

    Returns:
        str: Sequence from file
    """
    file = open(path_to_file, "r")
    lines = file.readlines()
    seq = ""
    for line in lines:
        seq += line.replace("\n", "")
    file.close()
    return seq

def write_seq_to_file(seq: str, path_to_file: str) -> None:
    """Writes sequence in file.

    Args:
        seq (str): Sequence to write
        path_to_file (str): Path to the file
    """
    file = open(path_to_file, "w")
    file.write(seq)
    file.close()
    return None

def validate_dna(seq: str) -> bool:
    """Checks if DNA sequence is valid. Returns True if sequence is valid, returns False if not.

    Args:
        seq (str): Sequence to be validated as DNA sequence or not
    """
    upper_seq = seq.upper()
    validator = upper_seq.count("A") + upper_seq.count("C") + upper_seq.count("T") + upper_seq.count("G")
    if validator == len(upper_seq): return True
    else: return False

def frequency (seq: str) -> dict:
    """Checks the frequency of symbols on a sequence and returns a dictionary with the keys being the symbols 
    and the values being the frequence of said symbol

    Args:
        seq (str): String sequence to be checked.

    Returns:
        dict: Dictionary with symbols as keys and values as frequencies.
    """
    dic = {}
    for sym in seq.upper():
        if sym in dic: dic[sym] += 1
        else: dic[sym] = 1
    return dic

def gc_content (seq: str) -> int:
    """Function that returns percentage of Guanine (G) and Cytosine (C) nucleotides in a DNA sequence

    Args:
        seq (str): DNA sequence_

    Returns:
        int: Guanine (G) and Cytosine (C) percentage
    """
    gc_count = 0
    for sym in seq.upper():
        if sym == "G" or sym == "C": gc_count += 1
    return gc_count / len(seq)

def gc_content_subseq (seq: str, size: int) -> List[int]:
    """Returns GC content of sub sequences of specified size

    Args:
        seq (str): Sequence
        size (int): Size of sub-sequences

    Returns:
        List[]: List of GC content of all sub-sequences of size size
    """
    result = []
    for i in range(0, len(seq)-size+1, size):
        subseq = seq[i:i+size]
        gc_ = gc_content(subseq)
        result.append(gc_)
    return result

def transcription(seq: str) -> str:
    """Computes the RNA sequence from the DNA sequence provided

    Args:
        seq (str): DNA sequence

    Returns:
        str: RNA sequence
    """
    assert validate_dna(seq), "Invalid DNA sequence"
    return seq.upper().replace("T", "U")

def reverse_complement(seq: str) -> str:
    """Computes the reverse complement of a DNA sequence provided

    Args:
        seq (str): DNA sequence

    Returns:
        str: Reverse complement sequence of DNA sequence provided
    """
    assert validate_dna(seq), "Invalid DNA sequence"
    complement = ""
    for sym in seq.upper():
        if sym == "A":
            complement += "T"
        elif sym == "T":
            complement += "A"
        elif sym == "C":
            complement += "G"
        elif sym == "G":
            complement += "C"
    return complement[::-1]

def translate_codon(codon: str) -> str:
    """Retrieves the codon corresponding aminoacid using a dictionary with the standard genetic code 

    Args:
        codon (str): DNA codon

    Returns:
        str: Aminoacid
    """
    DNA_Codons = {
    # 'M' - START, '_' - STOP
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "_", "TAG": "_", "TGA": "_"}

    if codon in DNA_Codons: return DNA_Codons[codon]
    else: return None

def translate_seq (seq: str, start_pos: int = 0) -> str:
    """Translates DNA sequence into an aminoacid sequence.

    Args:
        seq (str): DNA sequence

    Returns:
        str: Aminoacid sequence
    """
    amino_seq = ""
    codons = []
    seq = seq.upper()

    start = start_pos

    while start < len(seq):
        codons.append(seq[start:start+3])
        start += 3

    for codon in codons:
        if translate_codon(codon) == "_": break

        if translate_codon(codon) == None:
            continue
        else: amino_seq += translate_codon(codon)
    return amino_seq

def reading_frames(dna_seq: str) -> List[str]:
    """Computes de six reading frames of a DNA sequence including the reverse complement ones.

    Args:
        dna_seq (str): DNA sequence

    Returns:
        List[str]: List of all reading frames
    """
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    res = []
    res.append(translate_seq(dna_seq, 0)) # first index read
    res.append(translate_seq(dna_seq, 1)) # second index read
    res.append(translate_seq(dna_seq, 2)) # third index read
    rcomplement = reverse_complement(dna_seq)
    res.append(translate_seq(rcomplement, 0))
    res.append(translate_seq(rcomplement, 1))
    res.append(translate_seq(rcomplement, 2))
    return res

def all_proteins_rf(aa_seq: str) -> List[str]:
    """Computes all possible proteins in an aminoacid sequence.

    Args:
        aa_seq (str): Aminoacid sequence

    Returns:
        List[str]: List of possible proteins.
    """
    aa_seq = aa_seq.upper()
    proteins = []
    current_protein = []
    for aa in aa_seq:
        if aa == "M":
            current_protein.append("")
        for i in range(len(current_protein)):
            current_protein[i] += aa
        
        else:
            if aa == "_":
                if current_protein:
                    for p in current_protein:
                        proteins.append(p)
                    current_protein = []
    return proteins




