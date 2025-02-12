# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB
        
        # Initialize matrix private attributes for use in alignment
        # create matrices for alignment scores, gaps, and backtracing
        inf = float('inf') # Shout out to the IEEE floating point standard implementation
        mat_shape = (len(self._seqA) + 1, len(self._seqB) + 1)

        # Init matricies of alignment scores
        self._align_matrix = np.zeros(mat_shape)
        self._gapA_matrix = np.zeros(mat_shape)
        self._gapB_matrix = np.zeros(mat_shape)

        # Init matrices for backtrace procedure
        # 0: means start state, reserved for (0, 0)
        # 1: Go to alignment matrix
        # 2: Go to backA matrix
        # 3: Go to backB matrix
        # -1: This cell should not be reached
        self._back = np.zeros(mat_shape)
        self._back_A = np.zeros(mat_shape)
        self._back_B = np.zeros(mat_shape)

        self._align_matrix[0, 0] = 0
        self._gapA_matrix[0, 0] = self.gap_open
        self._gapB_matrix[0, 0] = self.gap_open
        for i in range(1, mat_shape[0]):
            # Match matrix: -inf, not applicable
            self._align_matrix[i, 0] = -inf
            self._back[i, 0] = -1
            # gapA: -inf, not applicable
            self._gapA_matrix[i, 0] = -inf
            self._back_A[i, 0] = -1
            # gapB: extend gap
            self._gapB_matrix[i, 0] = self._gapB_matrix[i-1, 0] + self.gap_extend
            self._back_B[i, 0] = 3
        for j in range(1, mat_shape[1]):
            # Match matrix: -inf, not applicable
            self._align_matrix[0, j] = -inf
            self._back[0, j] = -1
            # gapA: extend gap
            self._gapA_matrix[0, j] = self._gapA_matrix[0, j - 1] + self.gap_extend
            self._back_A[0, j] = 2
            # gapB: -inf, not applicable
            self._gapB_matrix[0, j] = -inf
            self._back_B[0, j] = -1


        # Calculate global alignment here
        for i in range(1, mat_shape[0]):
            for j in range(1, mat_shape[1]):
                # Update match matrix: match from previous match, extension on A, or extension on B
                # Uses highroad alignment: defaults to from gapA, then match, then gapB
                aa_pair = (self._seqA[i-1], self._seqB[j-1])
                self._align_matrix[i, j] = self._gapA_matrix[i-1, j-1] + self.sub_dict[aa_pair]
                self._back[i ,j] = 2
                if (self._align_matrix[i-1, j-1] + self.sub_dict[aa_pair]) > self._align_matrix[i, j]:
                    self._align_matrix[i, j] = self._align_matrix[i-1, j-1] + self.sub_dict[aa_pair]
                    self._back[i ,j] = 1
                if (self._gapB_matrix[i-1, j-1] + self.sub_dict[aa_pair]) > self._align_matrix[i, j]:
                    self._align_matrix[i, j] = self._gapB_matrix[i-1, j-1] + self.sub_dict[aa_pair]
                    self._back[i ,j] = 3

                # Update gapA matrix: insert gap from previous match or continue gap
                self._gapA_matrix[i, j] = self._align_matrix[i, j-1] + self.gap_open + self.gap_extend
                self._back_A[i, j] = 1
                if (self._gapA_matrix[i, j-1] + self.gap_extend) > self._gapA_matrix[i, j]:
                    self._gapA_matrix[i, j] = self._gapA_matrix[i, j-1] + self.gap_extend
                    self._back_A[i, j] = 2
                
                # Update gapB matrix: insert gap from previous match or continue gap
                self._gapB_matrix[i, j] = self._align_matrix[i-1, j] + self.gap_open + self.gap_extend
                self._back_B[i, j] = 1
                if (self._gapB_matrix[i-1, j] + self.gap_extend) > self._gapB_matrix[i, j]:
                    self._gapB_matrix[i, j] = self._gapB_matrix[i-1, j] + self.gap_extend
                    self._back_B[i, j] = 3
        
        # print(self._align_matrix)
        # print(self._gapA_matrix)
        # print(self._gapB_matrix)
        # print(self._back)
        # print(self._back_A)
        # print(self._back_B)
        		    
        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Determine the initial starting point (and alignment score)
        curr_i = len(self._seqA)
        curr_j = len(self._seqB)
        # Determine best alignment score (and starting end backtrace)
        self.alignment_score = self._align_matrix[curr_i, curr_j]
        curr_matrix = 1
        if (self._gapA_matrix[curr_i, curr_j] > self.alignment_score):
            self.alignment_score = self._gapA_matrix[curr_i, curr_j]
            curr_matrix = 2
        if (self._gapB_matrix[curr_i, curr_j] > self.alignment_score):
            self.alignment_score = self._gapB_matrix[curr_i, curr_j]
            curr_matrix = 3
        
        # Backtrace the alignment
        while curr_i > 0 or curr_j > 0:
            if curr_matrix not in (1, 2, 3):
                raise Exception(f"Reached a value other than 1, 2, or 3 at {curr_i}, {curr_j}")
            if curr_matrix == 1:
                self.seqA_align = self._seqA[curr_i-1] + self.seqA_align
                self.seqB_align = self._seqB[curr_j-1] + self.seqB_align
                curr_matrix = self._back[curr_i, curr_j]
                curr_i -= 1
                curr_j -= 1
            if curr_matrix == 2:
                self.seqA_align = '-' + self.seqA_align
                self.seqB_align = self._seqB[curr_j-1] + self.seqB_align
                curr_matrix = self._back_A[curr_i, curr_j]
                curr_i -= 0
                curr_j -= 1
            if curr_matrix == 3:
                self.seqA_align = self._seqA[curr_i-1] + self.seqA_align
                self.seqB_align = '-' + self.seqB_align
                curr_matrix = self._back_B[curr_i, curr_j]
                curr_i -= 1
                curr_j -= 0

        return (self.alignment_score, self.seqA_align, self.seqB_align)


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
