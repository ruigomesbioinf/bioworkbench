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

class BoyerMoore:
    def __init__(self, alphabet: str, pattern: str) -> None:
        self.alphabet = alphabet
        self.pattern = pattern
        self.preprocess()

    def preprocess(self):
        self.preprocess_bcr()
        self.preprocess_gsr()

    def preprocess_bcr(self):
        self.occ = {}
        for symb in self.alphabet:
            self.occ[symb] = -1
        for j in range(len(self.pattern)):
            c = self.pattern[j]
            self.occ[c] = j
        return self.occ

    def preprocess_gsr(self):
        self.f = [0] * (len(self.pattern) + 1)
        self.s = [0] * (len(self.pattern) + 1)
        i = len(self.pattern)
        j = len(self.pattern) + 1
        self.f[i] = j
        
        while i>0:
            while j <= len(self.pattern) and self.pattern[i-1] != self.pattern[j-1]:
                if self.s[j] == 0: self.s[j] = j-i
                j = self.f[j]
            i -= 1
            j -= 1
            self.f[i] = j
        j = self.f[0]
        for i in range(len(self.pattern)):
            if self.s[i] == 0: self.f[i] = j
            if i == j: j = self.f[j]

    def search_pattern(self, sequence: str) -> List[int]:
        res = []
        i = 0
        while i <= len(sequence) - len(self.pattern):
            j = len(self.pattern) - 1
            while j >= 0 and self.pattern[j] == sequence[j+i]: j -= 1
            if j < 0:
                res.append(i)
                i += self.s[0]
            else:
                c = sequence[j+i]
                i += max(self.s[j+1], j - self.occ[c])
        return res



def test():
    bm = BoyerMoore("ACTG", "ACCA")
    print(bm.search_pattern("ATAGAACCAATGAACCATGATGAACCATGGATACCCAACCACC"))

test()