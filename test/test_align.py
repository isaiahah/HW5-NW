# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    NW = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    alignment = NW.align(seq1, seq2)

    # Check output
    assert alignment[0] == -10 - 1 + NW.sub_dict[('M', 'M')] + NW.sub_dict[('Q', 'Q')] + NW.sub_dict[('R', 'R')]
    assert alignment[1] == "MYQR"
    assert alignment[2] == "M-QR"

    # Check matrices
    inf = float('inf')
    assert np.all(NW._align_matrix == np.array([[0, -inf, -inf, -inf], 
                                                [-inf, 5, -11, -13],
                                                [-inf, -12, 4, -8],
                                                [-inf, -12, -1, 5],
                                                [-inf, -14, -6, 4]]))
    assert np.all(NW._gapA_matrix == np.array([[-10, -11, -12, -13],
                                               [-inf, -inf, -6, -7],
                                               [-inf, -inf, -23, -7],
                                               [-inf, -inf, -23, -12],
                                               [-inf, -inf, -25, -17]]))
    assert np.all(NW._gapB_matrix == np.array([[-10, -inf, -inf, -inf],
                                               [-11, -inf, -inf, -inf],
                                               [-12, -6, -22, -24],
                                               [-13, -7, -7, -19],
                                               [-14, -8, -8, -6]]))
    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    NW = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    alignment = NW.align(seq3, seq4)

    # Check output
    assert alignment[0] == 17
    assert alignment[1] == "MAVHQLIRRP"
    assert alignment[2] == "M---QLIRHP"
