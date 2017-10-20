## EXERCISE 5: Environmental association analysis - outliers

Exercise 5 will focus on detecting specific loci that have strong associations with environment after adjusting for the overall association we observed in the previous exercise. Population genetic variation is influenced by demography, mating patterns, and natural selection, each of which is shaped by the environment. Theoretically, demographic changes and gene flow ("neutral" processes) affect genetic variation throughout the genome, whereas natural selection affects a relatively small number of genes. To disentangle which parts of the genome are under natural selection by the environment, we can use genome-wide data sets to identify SNPs/alleles that are extremely associated with environmental variables after accounting for the background association of genetic variation with the environment due to demographic and mating patterns. These "outliers" are candidate loci that could be under the influence of natural selection for local adaptation along the environmental gradient. We will use latent factor mixed models (LFMM), which is a multivariate method that can handle individual-based data (rather than "population"-based) and are considerd to be a powerful method for detecting loci under selection. For your own data, you could also consider linear mixed models such as EMMA or GEMMA when you have individual-based data ([example here](https://github.com/pgugger/LandscapeGenomics/blob/master/2015/Exercise4.mdown)), or [BayEnv2](http://gcbias.org/bayenv/) and [BayeScan](http://cmpg.unibe.ch/software/BayeScan/) when you have population-based data.

### Preparing the input

LFMM requires two input files: SNP data in 012 format with missing data coded as 9, and a table with environmental variables of interest in which rows are samples and column are variables. Neither can have headers or row names (the order will be preserved in the output).

Let's start by generating the SNP input file using the 012-formatted file that we generated in Exercise 4 with `vcftools` (in `~/Workshop/GBS_Data/Corrected_Output/`). Replace the default missing data symbol of -1 with 9 and remove the first column that has the sample numbers (you won't actually have any missing data if you followed my code exactly).

	sed 's/-1/9/g' snp.012 | cut -f2- > snp.lfmm

For the climate data, we can cut the relevant columns and remove the header from `clim.points`, which we saved in Exercise 4. View the first view lines of the file to see the column headings (`head`), and then `cut` the relevant columns and remove the header row with `tail`:

	head clim.points
	cut -f4-7 clim.points | tail -n+2 > clim.env
	
Copy the new SNP and climate data files to `~/Workshop/LFMM`. Note that instead of using climate variables directly, we could use principal components of the climate variables to generate a set of uncorrelated, synthetic climate variables. When dealing with correlated climate variables (as is typically the case), this strategy may be preferrable, but for today's exercise we will simply use the climate variables themselves as we did in previous analyses.
	
### Assessing population structure

LFMM accounts for overall or background associations of genetic variation with environmental variation using *latent factors* to model unobserved variation. A key step in LFMM is determining the number of latent factors to include, as this can effect the power of the test. It is advised to start by using the number of population clusters (*K*) inferred from a program like Structure or Admixture as the initial value, and then consider values slightly higher and lower. We will estimate *K* using a built-in function called `snmf`, which is conveniently part of the same LEA package for R that implements LFMM. If you want to see an example using Admixture, you can look at [last year's version of this exercise](https://github.com/pgugger/LandscapeGenomics/blob/master/2016/Exercise5.md).

Enter R, set up the work environment, and run `snmf` to estimate *K*.

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
	project = lfmm("snp.lfmm", "clim.env", K = 2, repetitions = 3, CPU = 16, iterations = 1000, burnin = 500, project = "new")

The command specifies the file name with the SNP data, the "best" number of clusters that we inferred above, the number of repetitions of the model to run, the number of processors to use, and the number of iterations and burn-in. Each can be adjusted accordingly, but you would likely want to run 5-10 repetitions and increase the iterations and burnin (perhaps 10-fold).

When the analysis finishes, we need to combine the data from the three repetitions and compute new calibrated *P*-values. To do that, first extract the *z*-scores for all repetitions for a given climate variable, then take the median. This can be done using the LEA function `z.scores` and the base function `apply` to take the median for each locus. Here is an example for the association tests with Pdry:

	z.pdry = z.scores(project, K = 2, d = 1)
	z.pdry <- apply(z.pdry, 1, median)

Next, we need to calculate lambda (the "genomic inflation factor"), which is commonly used for calibration of *P*-values. However, it is often considered too conservative, so some suggest using a value lower than lambda for the calibration. Lambda is calculated from the median of the median *z*-scores (from above) and a chi-squared distribution for each set of associations:
	
	lambda.pdry = median(z.pdry^2)/qchisq(0.5, df = 1)
	lambda.pdry

The calibrated or "adjusted" *P*-values are then calculated as follows:

	p.pdry.adj = pchisq(z.pdry^2/lambda.pdry, df = 1, lower = FALSE)
	
Now, repeat this correction procedure with the other three climate variables.
	
	z.pseas = z.scores(project, K = 2, d = 2)
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
	
The rows of the results files (and R objects) are ordered in the same order as the input data file, so they can easily be combined with other files you have generated (see `~/Worskhop/GBS_Data/Corrected_Output/snp.012.pos`) in or out of R. Outside of R, you can do it with the `paste` command in the Linux command line or in Excel. Exit R (`q()`) and then try making a single file with all the SNP locus names, *z*-scores, and *Q*-values yourself. Once you have succeeded, choose a SNP that is significant and use its name to find the fragment it came from in the FASTA file that was generated in the Stacks (see `~/Worskhop/GBS_Data/Corrected_Output/batch_1.samples-raw.fa`). Then, you could BLAST your locus to see if any functional information is available. Because of the low sample size in this exercise, the results are unlikely to be meaningful, but the approach can be readily applied your research projects.
