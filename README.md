# *Bioinformatics for Landscape Genomics*
Exercises for workshop organized by year, taught at UNAM-Morelia, Mexico unless otherwise noted in the folder name. Topics include Linux basics, Illumina data quality assessment, SNP calling with Stacks *de novo* pipeline, and analyzing data with landscape analyses such as GDM, GF, LFMM, or other approaches. I recommend using the most recent year as it reflects the latest versions of software.

## Required software

### On local computer

Windows: [PuTTY](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html) and FTP client (e.g., [FileZilla](https://filezilla-project.org/) or [WinSCP](https://winscp.net/eng/download.php))

MacOs: FTP client (e.g., [FileZilla](https://filezilla-project.org/))

### On server

[Stacks](http://catchenlab.life.illinois.edu/stacks/)

[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

[samtools](http://www.htslib.org/) 

[fastq_screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)

[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  (see "Requirements" section for dependencies)

[plink](https://www.cog-genomics.org/plink2)  (stable version)

[vcftools](https://vcftools.github.io/examples.html)

[admixture](http://software.genetics.ucla.edu/admixture/download.html)

[R]( https://www.r-project.org/) (optionally with [RStudio](https://www.rstudio.com/))

Then in R, run the following commands to install necessary packages (and any dependencies):
	
	install.packages(c("raster", "rgdal", "LEA", "vegan", "gdm"))
	install.packages("gradientForest", repos="http://R-Forge.R-project.org")
	install.packages("extendedForest", repos="http://R-Forge.R-project.org")
	if (!requireNamespace("BiocManager"))
	  install.packages("BiocManager")
	BiocManager::install()
	BiocManager::install(c("qvalue","DESeq2","WGCNA"))
