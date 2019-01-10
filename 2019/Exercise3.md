## EXERCISE 3: Landscape genomics - identifying loci under natural selection with LFMM

In this exercise, we will learn one approach [latent factor mixed modeling (LFMM)](https://academic.oup.com/mbe/article/30/7/1687/972098) to detect specific loci with strong environmental associations suggestive of natural selection for local adaptation. Population genetic variation is influenced by demography, mating patterns, and natural selection, each of which is shaped by the environment. Theoretically, demographic changes and gene flow ("neutral" processes) affect genetic variation throughout the genome, whereas natural selection affects a relatively small number of genes. To disentangle which parts of the genome are under natural selection by the environment, we can use genome-wide data sets to identify SNPs/alleles that are extremely associated with environmental variables after accounting for the background association of genetic variation with the environment due to demographic and mating patterns. These "outliers" are candidate loci that could be under the influence of natural selection for local adaptation along the environmental gradient. We will use latent factor mixed models (LFMM), which is a method that can handle individual-based data (rather than "population"-based) and is considered to be a powerful method for detecting loci under selection. For your own data, you could also consider multivariate ordination such as redundancy analysis (RDA) (*e.g.* [Forester *et al.* 2018](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.14584)) or [BayEnv2](https://bitbucket.org/tguenther/bayenv2_public/src) or [BayeScan](http://cmpg.unibe.ch/software/BayeScan/) when you have population-based data.

Refer to [Frichot *et al.* (2013)](https://academic.oup.com/mbe/article/30/7/1687/972098) for details of LFMM and to the LEA R package [publication](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12382) and [vignette](http://membres-timc.imag.fr/Olivier.Francois/LEA/files/LEA_github.pdf) for details of syntax and usage.

### Preparing the LFMM input

LFMM requires two input files: SNP data in 012 format with missing data coded as 9, and a table with environmental variables of interest in which rows are samples and column are variables. Neither can have headers or row names (the original order will be preserved in the output).

#### SNP data

Let's start by generating the SNP input file from the VCF file that you generated yesterday with Stacks. First, we will use `vcftools`, which is generally very useful for manipulating and filtering VCFs, to output the SNPs as genotypes encoded 0, 1, or 2 for homozygous, heterozygous, and other homozygous, respectively. 

	cd ~/Workshop/GBS_Data/Stacks_Output
	vcftools --vcf populations.snps.vcf --012 --out snp

Take a quick look at the output (file names with "012") to understand what each one represents. Then, the easiest thing to do would be to open the files and copy, paste, and edit them the way you want in a text editor or Excel. To stay in the command line environment and acheive the same goal, we can instead run the following commands: 

	sed 's/-1/9/g' snp.012 | cut -f2- > snp.lfmm
	
To summarize, we replaced (`sed`) any missing data encoded as -1 in of snp.012 to be encoded as 9 and then `cut` to keep from the second column to last column. View the resulting file with `less -S` and then move the file `snp.lfmm` to the LFMM folder.

#### Climate data

We will use [WorldClim](http://www.worldclim.org/) climate data, specifically temperature seasonality (bio4 = Tseas), minimum temperature of coldest month (bio6 = Tmin), precipitation seasonality (bio15 = Pseas), and precipitation of driest quarter (bio17 = Pdry). I put these climate data layers in the folder `Climate_Data` as GeoTiff format. Now, we will prepare the climate data for both today's and some of tomorrow's exercises. 

Enter R by typing `R` in the command line to enter the R command line or accessing the server through the RStudio GUI (may not be available depending on server set up).

First, set your working directory and load the required packages and the climate data layers.
	
	setwd("~/Workshop")
	library(raster)
	library(rgdal)
	clim.list <- dir("./Climate_Data/", full.names=T, pattern='.tif')  #makes list of file paths for each layer
	clim.layer <-  stack(clim.list)  #stacks the layers into a single object

We can view the climate data layers as maps simply with the `plot` function.

	pdf("clim.layer.pdf")
	plot(clim.layer)
	dev.off()

These commands save the plots as a PDF. If you are logged into the computer cluster, you will need to retrieve the PDF with `scp` (*e.g.*, FileZilla) to view it on your laptop. If you are running R locally on your laptop's graphical interface, then you can just run the `plot` commands without the `pdf` and `dev.off` commands, which will cause the plots to appear in a window. You can do this for all future plots in this exercise, but you will have to manually save them if you want to keep them.

Finally, we need to load our sample coordinates and extract climate data at each of the points.

	sample.coord <-read.table("sample.coord.txt", header=T, stringsAsFactors=F)
	sample.coord
	
	crs.wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"  #defines the spatial projection system that the points are in (usually WGS84)
	sample.coord.sp <- SpatialPointsDataFrame(sample.coord[,c('Longitude','Latitude')], proj4string=CRS(crs.wgs), data=sample.coord)

	clim.points <- extract(clim.layer, sample.coord.sp)  #extracts the data for each point (projection of climate layer and coordinates must match)
	clim.points <- cbind(sample.coord, clim.points)  #combines the sample coordinates with the climate data points
	write.table(clim.points, "~/Workshop/GF/clim.points", sep="\t", quote=F, row.names=F)  #save the table for use tomorrow with gradient forest (GF)
	clim.points 
	
	clim.env <- clim.points[, 4:7]
	colnames(clim.env) <- NULL
	clim.env
	write.table(clim.env, "~/Workshop/LFMM/clim.env", sep="\t", quote=F, row.names=F) #For LFMM, we will use this version without sample names, latitude, longitude, or column names.
	
Note that instead of using climate variables directly, we could use principal components of the climate variables to generate a set of uncorrelated, synthetic climate variables. When dealing with correlated climate variables (as is typically the case), this strategy may be preferrable, but for today's exercise we will simply use the climate variables themselves.
	
### Assessing population structure

LFMM accounts for overall/background associations of genetic variation with environmental variation using *latent factors* to model unobserved variation. A key step in LFMM is determining the number of latent factors to include, as this can effect the power of the test. It is advised to start by using the number of population clusters (*K*) inferred from a program like Structure or Admixture as the initial value, and then consider values slightly higher and lower. We will estimate *K* using a built-in function called `snmf`, which is conveniently part of the same LEA package for R that implements LFMM. If you want to see an example using Admixture, you can look at [a previous version of this exercise](https://github.com/pgugger/LandscapeGenomics/blob/master/2016/Exercise5.md).

Enter R, set up the work environment, and run `snmf` to estimate *K*, considering *K* from 1-4:

	setwd("~/Workshop/LFMM")
	library(LEA)
	
	project = NULL
	project = snmf("snp.lfmm", K = 1:4, entropy = TRUE, repetitions = 10, project = "new")
	pdf("sNMF.pdf")
	plot(project, col = "blue", pch = 19, cex = 1.2)
	dev.off()

The "best" value of *K* has the lowest cross-entropy value (*y*-axis). What is the best value according to your analysis?

It is also worth mentioning that the `snmf` as implemented in LEA offers the potential to impute missing data and make Structure-type plots. See the software manual for details, but here is the code to make a plot for K = 2:

	best = which.min(cross.entropy(project, K = 2))
	pdf("sNMF.barchart.pdf")
	barchart(project, K = 2, run = best, border = NA, space = 0, col = c("red", "blue"), xlab = "Individuals", ylab = "Ancestry proportions") -> bp
	axis(1, at = 1:length(bp$order), labels = bp$order, las=1, cex.axis = .3)
	dev.off()

####Short tangent on population structure-based outlier analysis 

Some of you may also be interested to know that this structure analysis (`snmf`) enables an an outlier analysis similar to *F*st-outlier analysis but the populations are based on the clusters (*K*) that we just inferred (assuming *K* > 1). We have not talked much about these outlier approaches, but they operate under a similar premise to environmental outlier analyses, except that we look for loci that are extremely differentiated among populations (that presumably have different habitats), rather than loci that are extremely correlated with environmental gradients. In case you would like to implement this analysis, you can run the following code (optional):

	p = snmf.pvalues(project, entropy = TRUE, ploidy = 2, K = 2)  #I chose K=2 here simply as an example and it may not match your conclusion on the "best" K from above
	pvalues = p$pvalues
	pdf("Differentiation_Outliers_K2.pdf")
	plot(-log10(pvalues), pch = 19, col = "blue", cex = .7, xlab = "SNP (ordered by contig arbitrarily)")
	dev.off()

This type of plot is a so-called *Manhattan plot* summarizing the significance (-log10(*P*)) of differentiation among clusters on the *y*-axis. The *x*-axis is usually ordered by chromosomal positions if you know the sequence of your organism's genome, but here we do not, so they are ordered arbitrarily in the same order that Stacks ordered them in the catalog. If you have time later, you could compare the outliers in this analysis to those that you get with the LFMM environmental associations below.

### Latent factor mixed modeling (LFMM)

Finally, we are ready to run LFMM. We can run LFMM simultaneously for all the climate variables contained in our climate data table (default) or one at a time. Run the following basic command

	project = NULL
	project = lfmm("snp.lfmm", "clim.env", K = 1, repetitions = 3, CPU = 16, iterations = 1000, burnin = 500, project = "new")

The command specifies the file name with the SNP data, the "best" number of clusters that we inferred above, the number of repetitions of the model to run, the number of processors to use, and the number of iterations and burn-in. Each can be adjusted according to your project, but you would likely want to run 5-10 repetitions and increase the iterations and burnin (perhaps 10-fold).

When the analysis finishes, we need to combine the data from the three repetitions and compute new calibrated *P*-values. To do that, first extract the *z*-scores for all repetitions for a given climate variable, then take the median. This can be done using the LEA function `z.scores` and the base function `apply` to take the median for each locus. Here is an example for the association tests with Pdry:

	z.pdry = z.scores(project, K = 1, d = 1)
	z.pdry <- apply(z.pdry, 1, median)

Next, we need to calculate lambda (the "genomic inflation factor"), which is commonly used for calibration of *P*-values. However, it is often considered too conservative, so some suggest using a value lower than lambda for the calibration. Lambda is calculated from the median of the median *z*-scores (from above) and a chi-squared distribution for each set of associations:
	
	lambda.pdry = median(z.pdry^2)/qchisq(0.5, df = 1)
	lambda.pdry

The calibrated or "adjusted" *P*-values are then calculated as follows:

	p.pdry.adj = pchisq(z.pdry^2/lambda.pdry, df = 1, lower = FALSE)
	
Now, repeat this correction procedure with the other three climate variables. Note that the value of `d` changes below to get the results for the second climate variable, which is Pseas. Be sure to change d accordingly to retrieve results for each climate variable.
	
	z.pseas = z.scores(project, K = 1, d = 2)
	z.pseas <- apply(z.pseas, 1, median)
	lambda.pseas = median(z.pseas^2)/qchisq(0.5, df = 1)
	p.pseas.adj = pchisq(z.pseas^2/lambda.pseas, df = 1, lower = FALSE)
	
	#etc.

To confirm that the model is behaving well with the *K* we chose and the adjustments to *P*, we need to inspect histograms of the *P*-values. The "best" *K* and proper calibration value will lead to histograms of *P*-values that are flat, except perhaps with an elevated frequency of very low *P*-values, representing the outliers. We can make all the histograms in a multi-paneled plot very simply with the base `hist` function.

	pdf("LFMM_P_Histograms.pdf")
	par(mfrow = c(4,1))
	hist(p.pdry.adj, col = "blue", main = "Pdry", xlab='')
	hist(p.pseas.adj, col = "blue", main = "Pseas", xlab='')
	hist(p.tmin.adj, col = "blue", main = "Tmin", xlab='')
	hist(p.tseas.adj, col = "blue", main = "Tseas", xlab=expression(italic(P)))
	dev.off()

How do these look? Refer back to today's lecture and the LEA/LFMM manual for guidance. If they suggest an overly conservative (right skew) or overly liberal (left skew) model/calibration, then we could repeat the analysis with a lower or higher value of *K*, respectively, or if there is slight right skew we might consider substituting a value lower than lambda to manually calibrate. Keep in mind that we are working with a very small sample size in this tutorial, so the patterns may not be typical.

Once we are convinced that the model is behaving well, we can move on. But, there is one final adjustment we need to make. We need to correct for multiple testing. We performed thousands of statistical tests (one per locus per climate variable), so many tests will appear significant by chance. The most common method of multiple testing correction is the *false discovery rate* (FDR) method of Benjamini and Hochberg (instead of Bonferroni correction, for example). In this process, we will adjust the *P*-values to *Q*-values. This correction can be easily implemented with the library `qvalue`.

	library(qvalue)
	q.pdry<-qvalue(p.pdry.adj)$qvalues
	q.pseas<-qvalue(p.pseas.adj)$qvalues
	q.tmin<-qvalue(p.tmin.adj)$qvalues
	q.tseas<-qvalue(p.tseas.adj)$qvalues

How does the number of significant tests based on *Q*-values (e.g., `sum(q.pdry<0.05)`) compare to the *P*-values (e.g., `sum(p.pdry.adj<0.05)`)?

A common way to visually summarize large numbers of association tests is using Manhattan plots (as we saw earlier with the outliers based on differentiation). All we need to do is plot -log10(Q) for each of the sets of association tests

	pdf("LFMM_Manhattan.pdf")
	par(mfrow = c(4,1))
	plot(-log10(q.pdry), pch = 19, col = "blue", cex = .7, xlab = '')
	plot(-log10(q.pseas), pch = 19, col = "blue", cex = .7, xlab = '')
	plot(-log10(q.tmin), pch = 19, col = "blue", cex = .7, xlab = '')
	plot(-log10(q.tseas), pch = 19, col = "blue", cex = .7, xlab = "SNP (ordered by contig arbitrarily)")
	dev.off()

-log10(0.01) = 2, so you can see that there are many extremely high values greater than 2, representing low *Q*-values less than 0.01. How many of the significant ones (*Q* < 0.05) are also of large effect? To answer this question, you can look at the set of significant SNPs that also have very high or very low *z*-scores.

	sum(q.pdry<0.01 & abs(z.pdry)>2)
	
You can also look for SNPs that have significant relationships with multiple climate variables, *e.g.*

	sum(q.pdry<0.01 & abs(z.pdry)>2 & q.pseas<0.01 & abs(z.pseas)>2)

Explore the results further in this way for this and the other climate variables. Finally, you might want to combine all the *z* and *Q*-values into a single table and then save for use in other software.

	lfmm.results <- cbind(z.pdry, q.pdry, z.pseas, q.pseas, z.tmin, q.tmin, z.tseas, q.tseas)
	head(lfmm.results)  #Note that the SNP locus numbers and positions are absent.
	
We can add the locus number and SNP positions indicated in the file `~/Worskhop/GBS_Data/Stacks_Output/snp.012.pos` that we generated earlier when we converted the SNPs to 012 format. Note that the SNP order should be the same because neither the file format conversion nor LFMM change the order.

	snp.names <-read.table("~/Worskhop/GBS_Data/Stacks_Output/snp.012.pos", header=F)
	colnames(snp.names) <- c("locus", "position")
	lfmm.results <- cbind(snp.names, lfmm.results)
	head(lfmm.results)  #Now we have a clear table with locus names and LFMM results
	
	write.table(lfmm.results, "lfmm.results", sep="\t", quote=F, row.names=F)
	
### Functional annotation

Determining the function or genomic context of interesting loci is challenging in non-model systems that lack genomic resources. Nonetheless, you may find some  functional information through [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) searches at [NCBI](https://www.ncbi.nlm.nih.gov/) or similar sequence databases. For some study systems, you may find closely related organisms that have annotated genome or transcriptome sequences, to which you can align your sequences directly and extract relevant annotation information. This workshop focuses on non-model systems that largely lack such resources, so let's focus on BLAST searches. 

First, examine one or two significant SNPs more closely by finding the sequence of the locus it comes from in the FASTA file that was previously generated in the Stacks (see `~/Worskhop/GBS_Data/Stacks_Output/populations.loci.fa`). Before exiting R, extract portions of the results table with significant SNPs and pick according to your interests, *e.g.*:

	lfmm.results[ which(lfmm.results$q.pdry < 0.01 & abs(lfmm.results$z.pdry) > 2), ]  #significant relationship with one specific climate variable
	lfmm.results[ which(lfmm.results$q.pdry < 0.01 & abs(lfmm.results$z.pdry) > 2 & lfmm.results$q.pseas < 0.01 & abs(lfmm.results$z.pseas) > 2), ]  #Significant with more two climate variables
	#You can also make "or" statements with "|" as the operator

Then, find the locus you selected in the FASTA file. Note that a locus numbered "14" in `lfmm.results` is written as "CLocus_14" in the FASTA file. Finally, copy and paste the relevant sequence into [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch) to see if it strongly matches any known sequences and whether any functional information is available. Because of the low sample size in this exercise, the results are unlikely to be meaningful, but the approach can be readily applied your research projects.

If you finish early, try scaling up the BLAST search to include *all* the significant loci (or a large subset of interest). *Hint*: you can generate a list of significant loci in R as suggested above (perhaps keep just the column with locus numbers), save the list to file, and then use `grep -A -f` (see `man grep` for details) in the Linux command line to extract the list of relevant sequences in the FASTA file and save that subset as a new FASTA file. You will likely want to adjust the locus names in the file from R to better match those in the FASTA file (*e.g.*, 14 *versus* CLocus_14). You should be able to then BLAST the entire FASTA file at once. 



