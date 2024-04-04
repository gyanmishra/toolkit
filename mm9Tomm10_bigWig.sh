# The below code has been taken from the BioStar post https://www.biostars.org/p/81185/

# Convert bigwig to bedGraph
ls *.bw | parallel --verbose  'bigWigToBedGraph {} {.}.mm9.bedGraph'

# Liftover from mm9 to mm10
ls *.mm9.bedGraph | parallel --verbose  \
'liftOver {} mm9ToMm10.over.chain {= s/.mm9.bedGraph// =}.mm10.bedGraph {.}.unMapped'

# Sort the bed file
ls *.mm10.bedGraph | parallel --verbose 'sort -k1,1 -k2,2n {} > {.}.sort.bedGraph'

# Fix the overlapping intervals
# The bed file from liftover might have overlapping intervals. You will hit error if directly using bedGraphToBigWig for file conversion. 
# In this case, you need to split the overlapping intervals and assign mean signal to a new bed file.
ls *.mm10.sort.bedGraph | \
parallel --verbose awk -vOFS=\"\\t\" \'{ print \$1, \$2, \$3, \".\", \$4 }\'  {} \> {.}.bed
ls *.mm10.sort.bed | \
parallel --verbose 'bedops --partition {} | bedmap --echo --mean --delim "\t" - {} >  {.}.split.bedGraph'

# Convert bedGraph to bigwig 
ls *.mm10.sort.split.bedGraph | \
parallel --verbose 'bedGraphToBigWig {} mm10.chrom.sizes {.}.bw'
