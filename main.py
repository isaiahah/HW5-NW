# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    NW = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    alignments = [NW.align(gg_seq, hs_seq), NW.align(mm_seq, hs_seq), 
                  NW.align(br_seq, hs_seq), NW.align(tt_seq, hs_seq)]
    score_species = [(alignments[0][0], "Gallus gallus"), 
                     (alignments[1][0], "Mus musculus"), 
                     (alignments[2][0], "Balaeniceps rex"), 
                     (alignments[3][0], "Tursiops truncatus")]
    score_species = sorted(score_species, reverse=True)

    # print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    for i, pair in enumerate(score_species):
        print(f"Rank {i + 1}: Species {pair[0]} with score {pair[1]}")
    

if __name__ == "__main__":
    main()
