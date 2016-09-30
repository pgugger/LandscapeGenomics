## EXERCISE 5: Environmental association analysis - outliers

Exercise 5 will focus on detecting specific loci that have strong associations with environment after adjusting for the overall association we observed in the previous exercise. Population genetic variation is influenced by demography, mating patterns, and natural selection, each of which is shaped by the environment. Theoretically, demographic changes and gene flow ("neutral" processes) affect genetic variation throughout the genome, whereas natural selection affects a relatively small number of genes. To disentangle which parts of the genome are under natural selection by the environment, we can use genome-wide data sets to identify SNPs/alleles that are extremely associated with environmental variables after accounting for the background association of genetic variation with the environment due to demographic and mating patterns. These "outliers" are candidate loci that could be under the influence of natural selection for local adaptation along the environmental gradient. We will use latent factor mixed models, which is a multivariate method that can handle individual-based data (rather than "population"-based) and are considerd to be a powerful method for detecting loci under selection. For your own data, you could also consider linear mixed models such as EMMA or GEMMA when you have individual-based data ([example here](https://github.com/pgugger/LandscapeGenomicsWorkshop_Morelia/blob/master/Exercise4.mdown)), or [BayEnv2](http://gcbias.org/bayenv/) and [BayeScan](http://cmpg.unibe.ch/software/BayeScan/) when you have population-based data.

### Preparing the input

LFMM requires two input files: SNP data in 012 format with missing data coded as -9, and a table with environmental variables of interest in which rows are samples and column are variables. Neither can have headers or row names (the order will be preserved in the output).

Let's start by generating the SNP input file using the 012-formatted file that we generated in Exercise 4 with `vcftools` (in `~/Workshop/GBS_Data/Corrected_Output/`). Replace the default missing data symbol of -1 with -9 and remove the first column that has the sample numbers (you won't actually have any missing data if you followed my code exactly).

	sed 's/-1/-9/g' snp.012 | cut -f2- > snp.lfmm

For the climate data, we can cut the relevant columns and remove the header from `clim.points`, which we saved from Exercise 4. View the first view lines of the file to see the column headings (`head`), and then `cut` the relevant columns and remove the header row with `tail`:

	head clim.points
	cut -f4-7 clim.points | tail -n+2 > clim.env
	
Copy the new SNP and climate data files to `~/Workshop/LFMM`.
	
### Assessing population structure

LFMM accounts for overall or background associations of genetic variation with environmental variation using *latent factors* to model unobserved variation. A key step in LFMM is determining the number of latent factors to include, as this can effect the power of the test. It is advised to start by using the number of population clusters (*K*) inferred from a program like Structure as the initial value, and then consider values slightly higher and lower. We will start by estimating *K* using the program [Admixture](https://www.genetics.ucla.edu/software/admixture/download.html).

First, we need to get the data in the right format. Admixture uses PLINK PED/MAP format in 12 format. It is similar to the 012 format we generated with `vcftools`, but has additional columns and is composed of rearranged files. The easiest way to generate the proper format is to start with the PLINK-formatted data that we generated directly from Stack in Exercise 3 (in `~/Workshop/GBS_Data/Corrected_Output/`). View `plink.ped` and `plink.map` to understand their structure. All we are going to do is run them through PLINK to change the ACGT format to 12 format.

	plink --file batch_1.plink --allow-extra-chr --recode 12 --out snp.plink.12

View the output to see the difference. Now we can run Admixture. As in Structure, one must consider a number of possible values of *K* and then choose among them. To run the basic command for a single value of *K* (let's say 5, for example), you could simply type

	admixture snp.plink.12.ped 5

However, it is easiest to run Admixture for all relevant values of *K* at once. Any thoughts how you might do that?  ... a `for` loop, of course! We also want to add `--cv` so that it will output the cross-validation error, which we will use to choose the "optimal" number of clusters, and `-j16` to enable parallel processing on 16 processors. Here is an example to run *K* = {1, 2, 3, 4} written as one line of code rather than a script:

	for K in 1 2 3 4; do admixture --cv snp.plink.12.ped $K -j16 | tee log${K}.out; done

After running Admixture for a given *K*, it sends the ouput that was shown on the screen to a log file using `tee`, each named accroding to the *K* value. View one of the log files to understand what it contains. We can quickly parse all of the log files for cross-validation error (CV), and choose the *K* with the lowest CV.

	grep -h CV log*.out

Which *K* has the lowest CV?

R code is available in the manual to quickly make the individual ancestry plots that you commonly see with similar programs, such as Structure. For now, we are most interested in estimating the number of clusters, so let's move on.


### Latent factor mixed modeling (LFMM)

We can run LFMM simultaneously for all the climate variables contained in our climate data table (default) or one at a time. The basic command would 

	LFMM -x snp.lfmm -v clim.env -K 1 -o results

Run the command and view the output files after it finishes one of the variables (you can stop it after that using `CTRL c`). You will see output files for each climate variable that contain *z*-scores and *P*-values for each SNP. The file name indicates which climate variable `s1` for the variable in the first column of `clim.env` and the *K* (`.1.`). See the output logged on the screen (or the manual) to understand the contents of the files.
	
However, we want to specify a few options, such as the number of processors to use in parallel (`-p 16`). You should also consider increasing the the burn-in length and the total number of iterations but we will use defaults today. If we have missing data, we would need to add `-m`. In addition, the developers suggest that LFMM be run 5-10 times and that the results then be combined. Again, we can use a loop to run the program repeatedly. First, make a folder called `results` and then run it in a loop 3 times, including the additional arguments I mentioned: 
	
	mkdir results
	for i in 1 2 3; do LFMM -x snp.lfmm -v clim.env -K 1 -p 16 -o results/res$i; done

While it is running read the relevant parts of the [manual](http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/note.pdf) to learn more, or help your neighbor. When it finishes, you could view all the outputs separately, but the LFMM developers provide R code to combine them and recalculate *z* and *P* values, adjusting for multiple testing using a false discovery rate (FDR) method.

Open R and set the working directory.

	R
	setwd("~/Workshop/LFMM/results")

Now we can copy and paste the developers code into the R command line:	
	
	z.table = NULL
	for (i in 1:3){
		file.name = paste("res",i, "_s1.1.zscore", sep="")
		z.table = cbind(z.table, read.table(file.name)[,1])
	}

	z.score = apply(z.table, MARGIN = 1, median)  #combines z-scores
	lambda = median(z.score^2)/0.456
	adjusted.p.values = pchisq(z.score^2/lambda, df = 1, lower = F)  #re-adjust p-values

Notice that a `for` loop is used to combine the *z*-scores by taking the median (the notation is different in R). Adjusted *P*-values are then calculated based on these and lambda, which is the "inflation factor" and summarizes the rate of false positives.

We can save the new *z* and *P* values to file. These represent the final summary of all results.

	write.table(z.score, file="LFMM_zscore_bio15.txt", sep="\t", quote=F, row.names=F)
	write.table(adjusted.p.values, file="LFMM_Pvalues_bio15.txt", sep="\t", quote=F, row.names=F)

Then, we can plot a histogram of the *P*-values to visually assess whether the value of *K* we chose leads to conservative or liberal tests. Use histograms and the estimate of lambda to decide whether the results are coservative, a good balance, or excess false positives.
	
	pdf("LFMM_Histogram_Pvalues_bio15.pdf")
	hist(adjusted.p.values, col = 3)
	dev.off()

	lambda  

Ideally, the histograms would show an even distribution except at very low *P* ("candidate SNPs") and lambda would equal 1. Should we consider changing *K* to better control false positives?

Finally, combine the SNP names (see `~/Worskhop/GBS_Data/Corrected_Output/snp.012.pos`), *z* and *P* values into a single table. This is most easily done by copying and pasting files in Excel, but you can also do it with the command line. Exit R (`q()`) and then try it yourself. You can count how many are significant and then use the SNP IDs to refer back to their sequence context in the FASTA output from Stacks. Then, you could BLAST these loci to see if any functional information is available. Because of the low sample size in this exercise, the results are unlikely to be meaningful, but the approach can be readily applied your research projects.