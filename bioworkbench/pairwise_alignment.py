# -*- coding: utf-8 -*-

# dotplots

def create_mat (nrows: int, ncols: int):
    """Creates a matrix filled with zeros with a given number of columns and rows

    Args:
        nrows (int): number of rows
        ncols (int): number of columns
    """
    mat = []
    for i in range(nrows):
        mat.append([])
        for j in range(ncols):
            mat[i].append(0)
            
    return mat

def dotplot (seq1: str, seq2: str):
    """Creates a matrix filled with zeros using the create_mat function and then place
    ones in the appropriate places where the characters in the two sequences are equal. 

    Args:
        seq1 (str): Sequence 1  
        seq2 (str): Sequence 2
    """
    mat = create_mat(len(seq1), len(seq2))
    for i in range(len(seq1)):
        for j in range(len(seq2)):
            if seq1[i] == seq2[j]:
                mat[i][j] = "*"
            else:
                mat[i][j] = " "
    return mat

def print_dotplot(mat: "create_mat", s1: str, s2: str):
    """Prints a visual friendly version of the dotplot

    Args:
        mat (create_mat): A matrix object
        s1 (str): Sequence 1
        s2 (str): Sequence 2
    """
    print(" ", s2)
    for i in range(len(mat)):
        print(s1[i], " ", *mat[i], sep="")
        

        
if __name__ == "__main__":
    seq1 = "ATGCGAC"
    seq2 = "CTGGGAC"
    a = dotplot(seq1, seq2)
    print_dotplot(a, seq1, seq2)