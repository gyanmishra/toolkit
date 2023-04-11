
Source : Biostar Post
#To build a BED file of intervals for your genes:
#wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/danRer7/database/refGene.txt.gz | gunzip -c > danRer7.refGene.txt
#grep -f genes.txt danRer7.refGene.txt | awk -v OFS="\t" '{ print $3,$5,$6,$13,$9,$4 }' | sort-bed - > danRer7.genes.bed

#To generate TSS regions, separate by strand and get the start position of the gene:

awk -v OFS="\t" '($6=="+"){ print $1, $2, ($2+1), $4, $5, $6; }' danRer7.genes.bed > danRer7.genes.start.for.bed
awk -v OFS="\t" '($6=="-"){ print $1, $3, ($3+1), $4, $5, $6; }' danRer7.genes.bed > danRer7.genes.start.rev.bed

Then apply a 5kb asymmetric padding around the start positions (you can pick a different window to define the regulatory region upstream of the TSS):

bedops --range -5000:0 --everything danRer7.genes.start.for.bed > danRer7.genes.start.for.TSS5k.bed
bedops --range 0:5000 --everything danRer7.genes.start.rev.bed > danRer7.genes.start.rev.TSS5k.bed

#Take the union of the two files:

bedops --everything danRer7.genes.start.for.TSS5k.bed danRer7.genes.start.rev.TSS5k.bed > danRer7.genes.TSS5k.bed

Separate out the padded TSSs per gene:

awk '{ print $0 >> "danRer7.TSS5k."$4".bed"; }' danRer7.genes.TSS5k.bed


Generate FASTA for each gene using the bed2faidx script. Here's how you would do this for one gene, say "LOC567192":

bed2faidx --fastaDir=/path/to/danRer7/indexed/fasta < danRer7.TSS5k.LOC567192.bed > danRer7.TSS5k.LOC567192.fa

