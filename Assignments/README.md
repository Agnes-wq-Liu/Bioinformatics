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

Fractions using annotation file as denominator:
match on both coordinates is: 0.536635284486.
match on start only is: 0.00777732296357.
match on end only is: 0.266475644699.
match on neither coordinates is: 0.189111747851.

Fractions using prediction file as denominator:
match on both coordinates is: 0.592592592593.
match on start only is: 0.00316169828365.
match on end only is: 0.289069557362.
match on neither coordinates is: 0.115176151762.

(e)prediction accuracy analysis: 

Annotated gene properties:

There are more forms of start codons and overlapping genes in the annotated file that might lead to completely un-detected gene regions. (supported by results of A3_1e_assessment.py on first sequence seen in vulnificus fasta file)
 
Also the parameters are generated from the species of Vibrio cholerae. Although this species is in the same genus as species of interest, there is a good chance for some parameters be significantly deviated from those of Vibrio vulnificus (eg, a specific start codon frequency being too low in cholerae than that of vulnificus), causing significant prediction biase. This is highly likely since there exists more start codons than my prediction parameters.

In addition, the annotation file has some genes with overlapping regions. Note some genes have start position smaller than previous stop position value, indicating overlapping of the genes, thus there might be some predictions that are partially correct.

My predicted gene properties:

By logic of the implementation of Viterbi Algorithm in 1b, all my genes start strictly with the inferred 3 start codons (ATG, TTG, GTG) and end strictly with the inferred 3 stop codons (TAA, TAG, TGA), which is not the case for annotated genes.

