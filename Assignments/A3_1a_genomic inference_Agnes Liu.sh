#/bin/bash
#preprocessing of the data: only considering CDS on + strand
sed '/CDS\|gff-version 3/!d' Vibrio_cholerae.GFC_11.37.gff3 > tmp.gff3
sed '/+/!d' tmp.gff3 > new_annot.gff3

#(i) average intergenic lengths
bedtools maskfasta -fi Vibrio_cholerae.GFC_11.dna.toplevel.fa -bed new_annot.gff3 -mc ' ' -fo intergenic.fa
sed '/A\|T\|C\|G/!d' intergenic.fa > new_inter.fa
tr -d '\n' < new_inter.fa > new_2.fa
tr -s ' ' < new_2.fa > new_3.fa
tr -s ' ' '\n' < new_3.fa > new_4.fa
#now im getting a file with each line a intergenic section
#do the similar operation as before
INT_CNT=$(< "new_4.fa" wc -l)
TOTL_INT_L=$(< "new_4.fa" wc -c)
AVG_INT_L=$((TOTL_INT_L/INT_CNT))
echo "(i) avg intergenic length" > 1a_output.txt
echo "intergenic lengths add up to: $TOTL_INT_L" >> 1a_output.txt
echo "intergenic section count=$INT_CNT" >> 1a_output.txt
echo "the average intergenic length is: $AVG_INT_L" >> 1a_output.txt

#(ii) avg gene lengths
#storing the line count (gene count)
GENE_CNT=$(< "new_annot.gff3" wc -l)
bedtools getfasta -fi Vibrio_cholerae.GFC_11.dna.toplevel.fa -bed new_annot.gff3 -fo genes.fa
sed '/DN/d' genes.fa > processed_genes.fa
TOTAL_GENE_L=$(< "processed_genes.fa" wc -c)
AVG_GENE_L=$((TOTAL_GENE_L/GENE_CNT))
echo \ >> 1a_output.txt
echo "(ii) avg gene length" >> 1a_output.txt
echo "gene lengths add up to: $TOTAL_GENE_L" >> 1a_output.txt
echo "gene count=$GENE_CNT" >> 1a_output.txt
echo "the average gene length is: $AVG_GENE_L" >> 1a_output.txt

#(iii) nucleotide frequency in intergenic regions
echo \ >> 1a_output.txt
echo "(iii) nucleotide frequencies for intergenic regions" >> 1a_output.txt
tr -d '\n' < new_4.fa > merged_inter.fa
freq_A=$(tr -cd 'A' < merged_inter.fa | wc -c)
A=$(echo $freq_A\/$TOTL_INT_L | bc -l)
echo "A: $A" >> 1a_output.txt
freq_T=$(tr -cd 'T' < merged_inter.fa | wc -c)
T=$(echo $freq_T\/$TOTL_INT_L | bc -l)
echo "T: $T" >> 1a_output.txt
freq_G=$(tr -cd 'G' < merged_inter.fa | wc -c)
G=$(echo $freq_G\/$TOTL_INT_L | bc -l)
echo "G: $G" >> 1a_output.txt
freq_C=$(tr -cd 'C' < merged_inter.fa | wc -c)
C=$(echo $freq_C\/$TOTL_INT_L | bc -l)
echo "C: $C" >> 1a_output.txt

#infer the codon frequency of genic regions
#first, let me clean up my files
rm new_inter.fa new_2.fa new_3.fa intergenic.fa tmp.gff3 genes.fa merged_inter.fa 
echo \ >> 1a_output.txt
echo "(iv) codon frequencies for genic regions" >> 1a_output.txt
#calculate start codon frequencies
echo "start codons" >> 1a_output.txt
cut -c-3 processed_genes.fa > sC.txt
totalstart=$(< "sC.txt" wc -l)
cp sC.txt tmp.txt
for ((i = 0 ; i < 3 ; i++)); do
  thisCodon=$(head -1 tmp.txt)
  sed -i "/$thisCodon/d" tmp.txt
  sed "/$thisCodon/d" sC.txt > new.txt
  tmpCount=$(< "new.txt" wc -l)
  rm new.txt
  thisCount=$((totalstart-tmpCount))
  this_freq=$(echo $thisCount\/$totalstart | bc -l)
  echo "$thisCodon: $this_freq" >> 1a_output.txt
done
rm tmp.txt

#stop codons
echo "stop codons" >> 1a_output.txt
rev processed_genes.fa > sttC.txt
cut -c-3 sttC.txt > stC1.txt
rev stC1.txt > stC.txt
totalstop=$(< "stC.txt" wc -l)
cp stC.txt tmp.txt
for ((i = 0 ; i < 3 ; i++)); do
  thisCodon=$(head -1 tmp.txt)
  sed -i "/$thisCodon/d" tmp.txt
  sed "/$thisCodon/d" stC.txt > new.txt
  tmpCount=$(< "new.txt" wc -l)
  rm new.txt
  thisCount=$((totalstop-tmpCount))
  this_freq=$(echo $thisCount\/$totalstop | bc -l)
  echo "$thisCodon: $this_freq" >> 1a_output.txt
done

echo "middle codon" >> 1a_output.txt
#finally, compute the middle codon frequencies
sed -i 's/.\{3\}/& /g' processed_genes.fa
cut -d ' ' -f 2- processed_genes.fa > tmp.txt
rev tmp.txt > newtmp.txt
rm tmp.txt
cut -d ' ' -f 3- newtmp.txt > tmp.txt
rev tmp.txt > middle_codons.txt
rm tmp.txt newtmp.txt
tr ' ' '\n' < middle_codons.txt > lineformat.txt

totalCodon=$(< "lineformat.txt" wc -l)
cp lineformat.txt tmp.txt
for ((i = 0 ; i < 61 ; i++)); do
  thisCodon=$(head -1 tmp.txt)
  sed -i "/$thisCodon/d" tmp.txt
  sed "/$thisCodon/d" lineformat.txt > new.txt
  tmpCount=$(< "new.txt" wc -l)
  rm new.txt
  thisCount=$((totalCodon-tmpCount))
  this_freq=$(echo $thisCount\/$totalCodon | bc -l)
  echo "$thisCodon: $this_freq" >> 1a_output.txt
done

#finally, remove all unecessary files
rm intergenic.txt tmp.txt lineformat.txt middle_codons.txt processed_genes.fa new_4.fa new_annot.gff3 stC* sC* sttC*
