Problem Summaries:

{A1}

When aligning two sequences, it is important to choose an alignment scoring scheme (substitution matrix cost and gap penalty function) that best captures the types of evolutionary events that are common in the type of sequences being aligned. 

Here, the Needleman-Wunsch dynamic programming algorithm is modified to identify the optimal alignment considering the multi-gap-free alignment problem. An alignment is multi-gap-free if it does not contain consecutive gaps in the same sequence. 
(implementation in python)

{A3}

(a) Inference based on Vibrio cholerae genome sequence & annotation: average lengths of intergenic regions and genic regions; nucleotide/codon frequencies for intergenic/genic regions.
Bash script with bedtools.

(b)&(c) Viterbi Algorithm Implementation & Test results with Vibrio vulnificus sequence.

(d)Fraction calculations of perfectly matched predictions, matched only at the beginning & only at the end.
