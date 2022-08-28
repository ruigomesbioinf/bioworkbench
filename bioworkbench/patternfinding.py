# -*- coding: utf-8 -*-

from typing import List

def search_first_occ(seq: str, pattern: str) -> int:
    """Return the position of the pattern or -1 if the pattern is not found.

    Args:
        seq (str): Sequence
        pattern (str): Pattern to search

    Returns:
        int: Integer of first position of the patter in the sequence
    """
    found = False
    i = 0
    while i <= len(seq)-len(pattern) and not found:
        j = 0
        while j < len(pattern) and pattern[j]==seq[i+j]:
            j += 1
        if j == len(pattern): found = True
        else: i += 1
    if found == True: return i
    else: return -1

def search_all_occurrences(seq: str, pattern: str) -> List[int]:
    """Returns all the positions of the patterns in the sequence

    Args:
        seq (str): Sequence
        pattern (str): Pattern to find

    Returns:
        List[int]: List of all the indexes
    """
    res = []
    for i in range(len(seq)-len(pattern)+1):
        j = 0
        while j < len(pattern) and pattern[j] == seq[i+j]:
            j += 1
        if j == len(pattern): res.append(i)

    return res

seqDNA = "ATAGAATAGATAATAGTC"
print ( search_all_occurrences(seqDNA, "AAT") )
