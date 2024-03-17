**1. `bamtoFASTQ.sh` : Extract original fastq files from aligned BAM files**

**2. `GenerateQualimap_report.py` : Generate table form RNAseq alignmnent report generated by [QualiMap](http://qualimap.conesalab.org/)**
```bash
$ qualimap rnaseq \
       	-bam out.bam \
	-gtf  gencode.vM10.annotation.gtf \
	-outdir qc_qualimap/out \
	-p non-strand-specific \
   	--java-mem-size=4G
$ python GenerateQualimap_report.py -d <qc_qualimap>
```

**3. `find_TSS_of_gene.sh` Extract TSS of gene from USCS refGene.**
The example shown in the script  is for zebrafish but can be implemented for other organism.

**4. [Converting mouse to human gene names and vice versa with biomaRt package](https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/)**

**5. Attach R version of conda environment to jupyter notebook.**
```bash
#Install nb_conda_kernels in the conda environment where jupyter is installed.
$ conda activate name_of_environment
$ conda install nb_conda_kernels
# Install IRkernel package in R
install.packages('IRkernel')
$ Rscript -e \ 'IRkernel::installspec(name="ir33",displayname="R 3.6.1")'
```

