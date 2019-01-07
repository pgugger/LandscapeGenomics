# *Bioinformatics for Landscape Genomics*
Exercises for January 2019 workshop taught at LANASE,  UNAM, Morelia, Mexico. Topics include Linux basics, Illumina data quality assessment, SNP calling with Stacks *de novo* pipeline, and analyzing and presenting data with landscape analyses such as GF, LFMM, or other approaches. 

## Required software

### On local computer

Windows: [PuTTY](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html) and FTP client (e.g., [FileZilla](https://filezilla-project.org/) or [WinSCP](https://winscp.net/eng/download.php))

MacOs: FTP client (e.g., [FileZilla](https://filezilla-project.org/))

### On server

[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  (see "Requirements" section for dependencies)

[fastq_screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)

[plink](https://www.cog-genomics.org/plink2)  (stable version)

[samtools](http://www.htslib.org/) 

[Stacks](http://catchenlab.life.illinois.edu/stacks/)

[vcftools](https://vcftools.github.io/examples.html)

[R]( https://www.r-project.org/) (optionally with [RStudio](https://www.rstudio.com/))

Then in R, run the following commands to install necessary packages (and any dependencies):
	
	install.packages(c("raster", "rgdal", "LEA", "vegan", "gdm"))
	install.packages("gradientForest", repos="http://R-Forge.R-project.org")
	install.packages("extendedForest", repos="http://R-Forge.R-project.org")
	if (!requireNamespace("BiocManager"))
	  install.packages("BiocManager")
	BiocManager::install()
	BiocManager::install("qvalue")
