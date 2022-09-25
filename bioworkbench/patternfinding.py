# -*- coding: utf-8 -*-

from typing import List
import re

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

def find_pattern_re(seq: str, pat: str):
    """Find pattern in string using regular expressions

    Args:
        seq (str): sequence
        pat (str): pattern to search on sequence (as regular expression)
    """
    mo = re.search(pat, seq)
    if mo != None:
        return mo.span()[0]
    else: return -1
    
def find_all_occurrences_re(seq: str, pat: str):
    """Returns all ocurrences of pattern in sequence using regular expressions to do it.

    Args:
        seq (str): sequence 
        pat (str): pattern to search on sequence (as regular expression)
    """
    mo = re.finditer(pat, seq)
    res = []
    for x in mo:
        res.append(x.span()[0])
    return res

def find_prosite(seq: str, pattern: str):
    """A function that turns a Prosite pattern into a regular expression for pattern searching

    Args:
        seq (str): Sequence where the search will happen
        pattern (str): Prosite pattern
    """
    regexp = pattern.replace("-", "")
    regexp = regexp.replace("x", ".")
    regexp = regexp.replace("{", "(")
    regexp = regexp.replace("}", ")")
    mo = re.search(regexp, seq)
    if mo != None:
        return mo.span()[0]
    else:
        return -1
    
def iub_to_re(iub: str):
    """Function that converts IUPAC extended alphabet patterns to regular expressions for string search

    Args:
        iub (str): IUPAC extended alphabet pattern
    """
    dic = {
        "A": "A",
        "C": "C",
        "G": "G",
        "T": "T",
        "U": "U",
        "K": "[GT]",
        "M": "[AC]",
        "R": "[AG]",
        "S": "[CG]",
        "W": "[AT]",
        "Y": "[CT]",
        "B": "[CGT]",
        "D": "[AGT]",
        "H": "[ACT]",
        "V": "[ACG]",
        "N": "[ACGT]"
    }
    
    site = iub.replace("^", "")
    regexp = ""
    
    for char in site:
        regexp += dic[char]
        
    return regexp

def cut_positions(enzyme: str, sequence: str):
    """Function that detects where a given enzyme will cut a given DNA sequence

    Args:
        enzyme (str): Restriction enzyme
        sequence (str): Sequence where given enzyme will cut
    """ 
    cut_pos = enzyme.find("^")
    regexp = iub_to_re(enzyme)
    
    matches = re.finditer(regexp, sequence)
    locs = []
    
    for m in matches:
        locs.append(m.start() + cut_pos)
        
    return locs

def cut_subsequences(locs: List[int], sequence: str):
    """Function that computes the resulting sub-sequences after restriction enzyme cuts (restriction map).

    Args:
        locs (List[int]): A list of locations of cutting positions
        sequence (str): DNA sequence
    """
    result = []
    positions = locs
    positions.insert(0, 0)
    positions.append(len(sequence))
    for i in range(len(positions) - 1):
        result.append(sequence[positions[i]:positions[i+1]])
        
    return result
    

def test():
    pos = cut_positions("G^ATTC", "GTAGAAGATTCTGAGATCGATTC")
    res = cut_subsequences(pos, "GTAGAAGATTCTGAGATCGATTC")
    print(pos)
    print(res)
    
test()


class BoyerMoore:
    """Heuristic algorithm for pattern searching on strings. This algorithm is based on two
    rules to make it work, which are implemented on this class
    """
    def __init__(self, alphabet: str, pattern: str) -> None:
        self.alphabet = alphabet
        self.pattern = pattern
        self.preprocess()

    def preprocess(self):
        """A function called by the class constructor that does the pre-+rocessing of both rules.
        """
        self.preprocess_bcr()
        self.preprocess_gsr()

    def preprocess_bcr(self):
        """Implementation of the Bad Character Rule, that ensures that we can advance our pattern
        to the next occurence (in the pattern) of the symbol in the sequence at the position of 
        the mismatch.
        """
        self.occ = {}
        for symb in self.alphabet:
            self.occ[symb] = -1
        for j in range(len(self.pattern)):
            c = self.pattern[j]
            self.occ[c] = j
        return self.occ

    def preprocess_gsr(self):
        """Implementation of the Good Suffix Rule, that ensures that when a mismatch occurs we can move
        to the next instance in the pattern of the part (suffix) that matched before (in the right) the mismatch.
        """
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
            if self.s[i] == 0: self.s[i] = j
            if i == j: j = self.f[j]

    def search_pattern(self, sequence: str) -> List[int]:
        """Functions that allows the user to search over target sequences for a given pattern
        using the boyer-moore rules.
        """
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

class Automata:
    """Class that implements the Deterministic Finite Automata for pattern searching
    """
    def __init__(self, alphabet: str, pattern: str) -> None:
        self.alphabet = alphabet
        self.numstates = len(pattern) + 1
        self.transition_table = {}
        self.build_transition_table(pattern)
        
    def build_transition_table(self, pattern: str) -> None:
        """Function called by the constructor that creates the transition table

        Args:
            pattern (str): pattern of search
        """
        for x in range(self.numstates):
            for char in self.alphabet:
                prefix = pattern[0:x] + char
                self.transition_table[(x, char)] = self.overlap(prefix, pattern)
                
    def overlap(self, s1: str, s2: str) -> int:
        """Function that provides the maximum overlap between two sequences s1 and s2.

        Args:
            s1 (str): sequence 1
            s2 (str): sequence 2

        Returns:
            int: maximum overlap integer
        """
        max_overlap = min(len(s1), len(s2))
        for i in range(max_overlap, 0, -1):
            if s1[-i:] == s2[:i]:
                return i
        return 0
    
    def print_automata (self):
        """Function that allows the user to print information about the automata and the transition table
        """
        print(f"States: {self.numstates}")
        print(f"Alphabet: {self.alphabet}")
        print("Transition table: ")
        for k in self.transition_table.keys():
            print(f"{k[0]}, {k[1]} -> {self.transition_table[k]}")
            
    def next_state(self, current, symbol):
        return self.transition_table.get((current, symbol))
    
    def apply_seq(self, seq: str):
        """Function that computes a list of states the DFA goes through when processinf one sequence.
        """
        q = 0
        res = [q]
        for char in seq:
            q = self.next_state(q, char)
            res.append(q)
        return res
    
    def occurrences_pattern(self, seq: str) -> List[int]:
        """Computes a list of the positions of the pattern in the sequence

        Args:
            seq (str): Sequence

        Returns:
            List[int]: List of positions
        """
        q = 0
        res = []
        for i in range(len(seq)):
            q = self.next_state(q, seq[i])
            if q == self.numstates-1:
                res.append(i - self.numstates + 2)
                
        return res
        
    