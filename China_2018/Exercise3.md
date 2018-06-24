## EXERCISE 3: Landscape genomic analysis

In this exercise, we will learn the basics of [gradient forest (GF)](https://doi.org/10.1890/11-0252.1) to quantify, describe, and map the associations of spatial and environmental variables with genetic variation, and [latent factor mixed modeling (LFMM)](https://academic.oup.com/mbe/article/30/7/1687/972098) to detect specific loci with strong environmental associations suggestive of natural selection for local adaptation. Refer to Fitzpatrick & Keller ([2015 *Ecology Letters*](http://dx.doi.org/10.1111/ele.12376)) for more detail about concepts and application of GF, and refer to the R package manuals and vignettes for details and examples of syntax and usage (*e.g.*, [GF](http://gradientforest.r-forge.r-project.org/biodiversity-survey.pdf)). Refer to Frichot *et al.* (2013) for details of LFMM and to the LEA R package vignette for details of syntax and usage (*e.g.*, [LEA](http://membres-timc.imag.fr/Olivier.Francois/LEA/files/LEA_github.pdf).

We will start with GF. GF is a multivariate, machine learning approach that considers potential non-linear patterns in the data. Other similar approaches that you might consider are [generalized dissimilarity modeling (GDM)](https://github.com/pgugger/LandscapeGenomics/blob/master/Exercise4.md), which is distance-based multivariate approach that also explores nonlinear relationships, or [redundancy analysis (RDA)](https://github.com/pgugger/LandscapeGenomics/blob/master/2015/Exercise4.mdown), which is a constrained ordination method that assumes linear relationships among variables.

Typically, GF is used with allele frequency data based on sample sites ("populations"), but we will use individual-based data and treat each of the 16 samples as data points. Otherwise, we wouldn't have enough data points to run the models.

### Gradient forest (GF)
#### Preparing the SNP data
First, we need to take the SNP output from Stacks and convert it to a convenient format for R. First, we will use `vcftools`, which is generally very useful for manipulating and filtering VCFs, to output the SNPs as genotypes encoded 0, 1, or 2 for homozygous, heterozygous, and other homozygous, respectively. 

	cd ~/Workshop/GBS_Data/Stacks_Output
	vcftools --vcf populations.snps.vcf --012 --out snp

Take a quick look at the output (file names with "012") to understand what each one represents. Then, the easiest thing to do would be to open the files and copy, paste, and edit them the way you want in a text editor or Excel. To stay in the command line environment and acheive the same goal, we can instead run the following commands: 

	cut -f2- snp.012 | sed 's/-1/NA/g' >snp.temp
	tr -d '\t' <snp.012.pos | tr '\n' '\t' | sed 's/[[:space:]]*$//' >header
	paste <(echo "ID" | cat - snp.012.indv) <(echo "" | cat header - snp.temp) > snp.forR
	rm header snp.temp

Don't worry about the details of these commands, but someday you may find it useful to learn them. To briefly summarize: we `cut` to keep from the second column to last column of snp.012, then replaced (`sed`) the missing data encoded as -1 to be encoded as NA (assuming we allowed missing data); then rearranged the list of SNP positions into a single tab-separated row to create a header; then `paste`d the SNP IDs with the header and SNP data. View the resulting file with `less -S`. When finished copy it to the `GF_GDM` folder with `cp snp.forR ~/Workshop/GF/`.

Now let's move to R by simply typing `R` in the command line. You are now in the command line environment of R, not the Linux bash shell itself. Some of you may be able to connect with the graphical interface of RStudio, or you could download the files and use the R graphical interface on your own computer, but in either case the commands will be the same.
	
Set the directory you will be working in.

	setwd("~/Workshop/GF")

Load the SNP data.

	snp <- read.table("snp.forR", header = T, row.names = 1)

#### Preparing the climate data
We will use WorldClim climate data, specifically temperature seasonality (bio4 = Tseas), minimum temperature of coldest month (bio6 = Tmin), precipitation seasonality (bio15 = Pseas), and precipitation of driest quarter (bio17 = Pdry). I put these climate data layers in the folder `Climate_Data` as GeoTiff format. First, we will load the required packages and the data.
	
	library(raster)
	library(rgdal)
	clim.list <- dir("./Climate_Data/", full.names=T, pattern='.tif')  #makes list of file paths for each layer
	clim.layer <-  stack(clim.list)  #stacks the layers into a single object

Now, we will crop the climate data layers to just the area of interest. I defined an extent with minimum and maximum longitude and minimum and maximum latitude, allowing us to then use the `crop` function.

	extent <- c(-104, -96, 18, 22) 
	clim.layer.crop <- crop(clim.layer, extent)

We can compare before and after by simply plotting the climate data.

	pdf("clim.layer.BandA.pdf")
	plot(clim.layer)
	plot(clim.layer.crop)
	dev.off()

These commands save the plots as a PDF. If you are logged into the computer cluster, you will need to retrieve the PDF with `scp` (*e.g.*, WinSCP) to view it on your laptop. If you are running R locally on your laptop' graphical interface, then you can just run the `plot` commands without the `pdf` and `dev.off` commands, which will cause the plots to appear in a window. You can do this for all future plots in this exercise, but you will have to manually save them if you want to keep them.

Finally, we need to load our sample coordinates and extract climate data at each of the points.

	sample.coord <-read.table("sample.coord.txt", header=T, stringsAsFactors=F)
	sample.coord
	
	crs.wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"  #defines the spatial projection system that the points are in (usually WGS84)
	sample.coord.sp <- SpatialPointsDataFrame(sample.coord[,c('Longitude','Latitude')], proj4string=CRS(crs.wgs), data=sample.coord)

	clim.points <- extract(clim.layer.crop, sample.coord.sp)  #extracts the data for each point (projection of climate layer and coordinates must match)
	clim.points <- cbind(sample.coord, clim.points)  #combines the sample coordinates with the climate data points
	write.table(clim.points, "clim.points", sep="\t", quote=F, row.names=F)  #save the table for later use
	clim.points


#### Gradient forest (GF) analysis

We will model the associations of spatial and climate variables with allele frequencies (genotypes) of individuals. You have already seen how we got the climate and genetic data. For the spatial variables, one could use latitude and longitude, but a more sophisticated approach might be to use PCNMs or MEMs (principal coordinates of neighbor matrices or Moran's eigenvector maps). These approaches generate a set of uncorrelated spatial variables. We do not have time to go into these methods in detail, but here is the code to generate the PCNM spatial variables.

	library(vegan)
	coord <- clim.points[,c("Longitude","Latitude")]
	pcnm <- pcnm(dist(coord))  #this generates the PCNMs, you could stop here if you want all of them
	keep <- round(length(which(pcnm$value > 0))/2)
	pcnm.keep <- scores(pcnm)[,1:keep]  #keep half of positive ones as suggested by some authors
	pcnm.keep

Now we are ready to begin GF. Load the library and then create a file that contains only the climate and PCNM spatial variables (no lat/lon).

	library(gradientForest)
	env.gf <- cbind(clim.points[,c("Pdry", "Pseas", "Tmin", "Tseas")], pcnm.keep)

In GF, a maximum number of splits can be defined following the developers suggestion

	maxLevel <- log2(0.368*nrow(env.gf)/2)

Running GF involves only one command:

	gf <- gradientForest(cbind(env.gf, snp), predictor.vars=colnames(env.gf), response.vars=colnames(snp), ntree=500, maxLevel=maxLevel, trace=T, corr.threshold=0.50)

The input is the combined climate, spatial and SNP data as input (`cbind(env.gf, snp)`), and the subsequent parts of the command define which variables are predictors and response variables, as well as a number of other parameters that I have left as suggested. When it finishes, there will be warnings about having less than five values for response variables, which is because we have only three: 0, 1, or 2. You can ignore them.

We can plot bar graphs depicting the importance of each spatial and climate variable.

	pdf("GF_VariableImportance.pdf")
	plot(gf, plot.type = "O")
	dev.off()

Which variables are most important?

We can also plot the "turnover functions" showing how allelic composition changes along the spatial or environmental gradients. The shapes are nonlinear and large jumps show steep genetic changes along certain portions of the environmental gradient. The height that the function acheives on the right side of the plot is the total importance and should match the barplot. First, organize the variables by importance and then plot:

	by.importance <- names(importance(gf))

	pdf("GF_TurnoverFunctions.pdf")
	plot(gf, plot.type = "C", imp.vars = by.importance, show.species = F, common.scale = T, cex.axis = 1, cex.lab = 1.2, line.ylab = 1, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 2, 2, 2), omi = c(0.2, 0.3, 0.2, 0.4)))
	dev.off()

Do you see any interesting patterns? It appears genetic variation changes abruptly for Pseas values of 60 and then acheives high cumulative importance.

We can also make plots of turnover functions for individual loci:

	pdf("GF_TurnoverFunctions_bySNP.pdf")
	plot(gf, plot.type = "C", imp.vars = by.importance, show.overall = F, legend = T, leg.posn = "topleft", leg.nspecies = 5, cex.lab = 0.7, cex.legend = 0.4, cex.axis = 0.6, line.ylab = 0.9, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0, 0.3, 0, 0)))
	dev.off()

Each line within each panel represents allelic change at a single SNP. Notice that in each panel some SNPs show very steep changes along the environmental gradient. One might consider these SNPs as candidates for involvement in local adaptation along the gradient, as we will discuss tomorrow. This approach to "outlier" detection is still being tested, but Fitzpatrick & Keller (2015) show a promising example.

Several other plots and tables can be output from `gradientForest`. I recommend generating these with your own data sets, but we will skip them today.

One of the useful features of these models is that the results can be readily projected spatially across a landscape that was not necessarily completely sampled. The mapping can be a little tricky, but first we must extract the climate data for **all** the cells in the climate data layer and remove missing data.

	clim.land <- extract(clim.layer.crop, 1:ncell(clim.layer.crop), df = TRUE)
	clim.land <- na.omit(clim.land)

Then, we can use the `predict` function to "predict" the allelic turnover across the whole landscape using our fitted model (`gf`) and the climate values on the landscape (`clim.land`).

	pred <- predict(gf, clim.land[,-1])  #note the removal of the cell ID column with [,-1])

These predictions then need to be converted to a color scale for mapping. One way is to use principal components analysis (PCA) on the predictions and use the first three axes to define red, green, blue color scales, respectively. After the values are defined in color space, they can be stacked and mapped. 

	PCs <- prcomp(pred, center=T, scale.=F)
	r <- PCs$x[, 1]
	g <- PCs$x[, 2]
	b <- PCs$x[, 3]
	r <- (r - min(r))/(max(r) - min(r)) * 255
	g <- (g - min(g))/(max(g) - min(g)) * 255
	b <- (b - min(b))/(max(b) - min(b)) * 255
	mask<-clim.layer.crop$Pseas
	mask[]<-as.numeric(mask[]>0)
	rastR <- rastG <- rastB <- mask
	rastR[clim.land$ID] <- r
	rastG[clim.land$ID] <- g
	rastB[clim.land$ID] <- b
	rgb.rast <- stack(rastR, rastG, rastB)

Finally, here is how to make the map, along with the original data points added for reference:

	pdf("GF_Map.pdf")
	plotRGB(rgb.rast, bgalpha=0)
	points(clim.points$Longitude, clim.points$Latitude)
	dev.off()

The colors represent genetic variation as predicted based on the modeled relationships with environmental and spatial variables. Similar colors are more similar genetically.

### Latent factor mixed modeling (LFMM)

Now, we will focus on detecting specific loci that have strong associations with environment after adjusting for the overall association we observed in the previous exercise. Population genetic variation is influenced by demography, mating patterns, and natural selection, each of which is shaped by the environment. Theoretically, demographic changes and gene flow ("neutral" processes) affect genetic variation throughout the genome, whereas natural selection affects a relatively small number of genes. To disentangle which parts of the genome are under natural selection by the environment, we can use genome-wide data sets to identify SNPs/alleles that are extremely associated with environmental variables after accounting for the background association of genetic variation with the environment due to demographic and mating patterns. These "outliers" are candidate loci that could be under the influence of natural selection for local adaptation along the environmental gradient. We will use latent factor mixed models (LFMM), which is a multivariate method that can handle individual-based data (rather than "population"-based) and are considerd to be a powerful method for detecting loci under selection. For your own data, you could also consider multivariate ordination such as RDA (*e.g.* [Forester *et al.* 2018](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.14584)) or linear mixed models such as EMMA or GEMMA when you have individual-based data ([example here](https://github.com/pgugger/LandscapeGenomics/blob/master/2015/Exercise4.mdown)), or [BayEnv2](http://gcbias.org/bayenv/) and [BayeScan](http://cmpg.unibe.ch/software/BayeScan/) when you have population-based data.

### Preparing the LFMM input

LFMM requires two input files: SNP data in 012 format with missing data coded as 9, and a table with environmental variables of interest in which rows are samples and column are variables. Neither can have headers or row names (the order will be preserved in the output).

Let's start by generating the SNP input file using the 012-formatted file that we generated above with `vcftools` (in `~/Workshop/GBS_Data/Stacks_Output/`). Replace the default missing data symbol of -1 with 9 and remove the first column that has the sample numbers.

	sed 's/-1/9/g' snp.012 | cut -f2- > snp.lfmm

For the climate data, we can cut the relevant columns and remove the header from `clim.points`, which we saved above in during the GF analysis. View the first view lines of the file to see the column headings (`head`), and then `cut` the relevant columns and remove the header row with `tail`:

	head clim.points
	cut -f4-7 clim.points | tail -n+2 > clim.env
	
Copy the new SNP and climate data files to `~/Workshop/LFMM`. Note that instead of using climate variables directly, we could use principal components of the climate variables to generate a set of uncorrelated, synthetic climate variables. When dealing with correlated climate variables (as is typically the case), this strategy may be preferrable, but for today's exercise we will simply use the climate variables themselves as we did in previous analyses.
	
### Assessing population structure

LFMM accounts for overall or background associations of genetic variation with environmental variation using *latent factors* to model unobserved variation. A key step in LFMM is determining the number of latent factors to include, as this can effect the power of the test. It is advised to start by using the number of population clusters (*K*) inferred from a program like Structure or Admixture as the initial value, and then consider values slightly higher and lower. We will estimate *K* using a built-in function called `snmf`, which is conveniently part of the same LEA package for R that implements LFMM. If you want to see an example using Admixture, you can look at [a previous version of this exercise](https://github.com/pgugger/LandscapeGenomics/blob/master/2016/Exercise5.md).

Enter R, set up the work environment, and run `snmf` to estimate *K*, considering *K* from 1-4.

	setwd("~/Workshop/LFMM")
	library(LEA)
	
	project = NULL
	project = snmf("snp.lfmm", K = 1:4, entropy = TRUE, repetitions = 10, project = "new")
	pdf("SNMF.pdf")
	plot(project, col = "blue", pch = 19, cex = 1.2)
	dev.off()

The "best" value of *K* has the lowest cross-entropy value (*y*-axis). What is the best value according to your analysis?

At the risk of going off on a tangent, some of you may also be interested to know that this structure analysis (`snmf`) enables an an outlier analysis similar to *F*st-outlier analysis but the populations are based on the clusters (*K*) that we just inferred (assuming *K* > 1). We have not talked much about these outlier approaches, but they operate under a similar premise to environmental outlier analyses, except that we look for loci that are extremely differentiated among populations (that presumably have different habitats), rather than loci that are extremely correlated with environmental gradients. In case you would like to implement this analysis, you can run the following code:

	p = snmf.pvalues(project, entropy = TRUE, ploidy = 2, K = 2)
	pvalues = p$pvalues
	pdf("Differentiation_Outliers_K2.pdf")
	plot(-log10(pvalues), pch = 19, col = "blue", cex = .7, xlab = "SNP (ordered by contig arbitrarily)")
	dev.off()

This type of plot is a so-called *Manhattan plot* summarizing the significance (-log10(*P*)) of differentiation among clusters on the *y*-axis and ordering the loci arbitrarily (in this case) along the *x*-axis. If you have time at the end, you could compare the outliers in this analysis to those that you get with the LFMM environmental associations below.

It is also worth mentioning that the `snmf` as implemented in LEA offers the potential to impute missing data and make Structure-type plots. See the software manual for details.

### Latent factor mixed modeling (LFMM)

We can run LFMM simultaneously for all the climate variables contained in our climate data table (default) or one at a time. Run the following basic command

	project = NULL
	project = lfmm("snp.lfmm", "clim.env", K = 1, repetitions = 3, CPU = 8, iterations = 1000, burnin = 500, project = "new")

The command specifies the file name with the SNP data, the "best" number of clusters that we inferred above, the number of repetitions of the model to run, the number of processors to use, and the number of iterations and burn-in. Each can be adjusted accordingly, but you would likely want to run 5-10 repetitions and increase the iterations and burnin (perhaps 10-fold).

When the analysis finishes, we need to combine the data from the three repetitions and compute new calibrated *P*-values. To do that, first extract the *z*-scores for all repetitions for a given climate variable, then take the median. This can be done using the LEA function `z.scores` and the base function `apply` to take the median for each locus. Here is an example for the association tests with Pdry:

	z.pdry = z.scores(project, K = 1, d = 1)
	z.pdry <- apply(z.pdry, 1, median)

Next, we need to calculate lambda (the "genomic inflation factor"), which is commonly used for calibration of *P*-values. However, it is often considered too conservative, so some suggest using a value lower than lambda for the calibration. Lambda is calculated from the median of the median *z*-scores (from above) and a chi-squared distribution for each set of associations:
	
	lambda.pdry = median(z.pdry^2)/qchisq(0.5, df = 1)
	lambda.pdry

The calibrated or "adjusted" *P*-values are then calculated as follows:

	p.pdry.adj = pchisq(z.pdry^2/lambda.pdry, df = 1, lower = FALSE)
	
Now, repeat this correction procedure with the other three climate variables.
	
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

How do these look? Refer back to today's lecture and the LEA/LFMM manual for guidance. If they suggest an overly conservative (right skew) or overly liberal (left skew) model/calibration then we could repeat the analysis with a lower or higher value of *K*, respectively, or if there is slight right skew we might consider substituting a value lower than lambda to manually calibrate. Keep in mind that we are working with a very small sample size in this tutorial, so the patterns may not be typical.

Once we are convinced that the model is behaving well, we can move on. But, there is one final adjustment we need to make. We need to correct for multiple testing. We performed thousands of statistical tests (one per locus per climate variable), so many tests will appear significant by chance. The most common method of multiple testing correction is the *false discovery rate* (FDR) method of Benjamini and Hochberg (instead of Bonferroni correction, for example). In this process, we will adjust the *P*-values to *Q*-values. This correction can be easily implemented with the library `qvalue`.

	library(qvalue)
	q.pdry<-qvalue(p.pdry.adj)$qvalues
	q.pseas<-qvalue(p.pseas.adj)$qvalues
	q.tmin<-qvalue(p.tmin.adj)$qvalues
	q.tseas<-qvalue(p.tseas.adj)$qvalues

How does the number of significant tests based on *Q*-values (e.g., `sum(q.pdry<0.05)`) compare to the *P*-values (e.g., `sum(p.pdry.adj<0.05)`)?

A common way to visually summarize large numbers of association tests is using a Manhattan plots, as we saw earlier with the outliers based on differentiation. All we need to do is plot -log10(Q) for each of the sets of association tests

	pdf("LFMM_Manhattan.pdf")
	par(mfrow = c(4,1))
	plot(-log10(q.pdry), pch = 19, col = "blue", cex = .7, xlab = '')
	plot(-log10(q.pseas), pch = 19, col = "blue", cex = .7, xlab = '')
	plot(-log10(q.tmin), pch = 19, col = "blue", cex = .7, xlab = '')
	plot(-log10(q.tseas), pch = 19, col = "blue", cex = .7, xlab = "SNP (ordered by contig arbitrarily)")
	dev.off()

-log10(0.01) = 2, so you can see that there are many extremely high values greater than 2, representing low *Q*-values less than 0.01. How many of the significant ones (*Q* < 0.05) are also of large effect? To answer this question, you can look at the set of significant SNPs that also have very high or very low *z*-scores.

	sum(q.pdry<0.01 & abs(z.pdry)>2)
	
Explore the results in this way for this and the other climate variables. Finally, you might want to combine all the *z* and *Q*-values into a single table and then save for use in other software.

	lfmm.results <- cbind(z.pdry, q.pdry, z.pseas, q.pseas, z.tmin, q.tmin, z.tseas, q.tseas)
	write.table(lfmm.results, "lfmm.results", sep="\t", quote=F, row.names=F)
	
The rows of the results files (and R objects) are ordered in the same order as the input data file, so they can easily be combined with other files you have generated (see `~/Worskhop/GBS_Data/Stacks_Output/snp.012.pos`) in or out of R. Outside of R, you can do it with the `paste` command in the Linux command line or in Excel. Exit R (`q()`) and then try making a single file with all the SNP locus names, *z*-scores, and *Q*-values yourself. Once you have succeeded, choose a SNP that is significant and use its name to find the fragment it came from in the FASTA file that was generated in the Stacks (see `~/Worskhop/GBS_Data/Stacks_Output/populations.loci.fa`). Then, you could [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch) your locus to see if any functional information is available. Because of the low sample size in this exercise, the results are unlikely to be meaningful, but the approach can be readily applied your research projects.

