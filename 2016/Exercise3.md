## EXERCISE 3: Advanced Stacks & optimizing parameters

This exercise builds on the basic pipeline to add steps that ensure high quality results. You will then learn how to optimize parameters.

### The full pipeline

In Exercise 2, you learned the basic pipeline, but more steps are recommended to enhance the quality of the Stacks results. This full pipeline is as follows:

1. `process_radtags`
2. `ustacks` or `pstacks`
3. `cstacks`
4. `sstacks`
5. `rxstacks`
6. `cstacks`
7. `sstacks`
8. `populations`
9. `load_radtags.pl`
10. `index_radtags.pl`
11. Examine output files and MySQL database, filter, decide about blacklisting, decide about changing parameters or re-running pipeline, and consider `export_sql.pl` and  rerunning `populations`

You know Steps 1-4. Steps 9-11 above relate to using a MySQL database to visualize and further filter results. We will not go over Steps 9-11 in this workshop, but I have provided you with some example commands in `MySQL.md`. 

So, let's start with Step 5. `rxstacks` will correct the haplotype/genotype calls of individuals using information from the whole population. It especially will help filter out calls due to sequencing errors and low coverage, as well as cases when multiple loci in an individual map to multiple catalog loci or there are an excessive number of haplotypes at a locus. 

To avoid overwriting previous results, make a new output folder: `mkdir ~/Workshop/GBS_Data/Corrected_Output`. Then, try running this command on your results from yesterday:

	rxstacks -b 1 -P ~/Workshop/GBS_Data/Stacks_Output/ -o ~/Workshop/GBS_Data/Corrected_Output/ --conf_filter --conf_lim 0.25 --prune_haplo --model_type bounded --bound_high 0.05 --lnl_filter --lnl_lim -50.0 -t 16

While it runs, look at the terminal to see what it is filtering and correcting. In the future, if you want it to save this information in a log file, instead of simply watching it flash across the screen, you can add the following to the end of the command next time you run it: `&>> ~/Workshop/GBS_Data/Corrected_Output/Log`. You can use this type of notation on other steps too, if you wish.  As `rxstacks` continues to run, read the [webpage](http://creskolab.uoregon.edu/stacks/comp/rxstacks.php) to learn more about the parameters. Some `rxstacks` parameters can have a large effect on the number of loci that remain after filtering. For `--lnl_lim`, -10.0 is a recommneded starting point but I have found it to be too conservative for our data. The Stacks web manual for `rxstacks` provides code to make a table that can be used to make a histogram of the log-likelihoods, which might help you choose a threshold.

When `rxstacks` finishes, notice that it has written new files for each sample to the output folder. Some look similar to those output by `ustacks` in the first step, but others are new, such as a list of blacklisted loci that were filtered. 

Now, we need to (re)run `cstacks`, `sstacks`, and `populations` on the new filtered data files (Steps 6-8 above). We will use the same settings as we did prior to `rxstacks`, except that the input data are now in `~/Workshop/GBS_Data/Corrected_Output` not `~/Workshop/GBS_Data` or `~/Workshop/GBS_Data/Stacks_Output`, and we will send the output to Corrected_Output. Modify your loop scripts from yesterday for `cstacks` and `sstacks` to adjust the input and output and include all 16 samples. Run these scripts and then run `populations` on the resulting output (total runtime could be more than 1 hour). We will use the resulting filtered SNPs for the next exercise. 

While those are running modify your "master" script to include `rxstacks`, plus your additional `cstacks` and `sstacks` commands.

When `populations` finishes, look at the summary statistics file. How do the results differ from the unfiltered run we did yesterday? Another way to check if the results make sense is to run an ordination (MDS using Hamming distance or PCA) and check that individuals from the same population cluster with each other. We can do that seamlessly by using [PLINK](https://www.cog-genomics.org/plink2), which is another genomics software package separately installed. Here are the commands you can run MDS (a.k.a PCoA) for the output pre-filtering and post-filtering. 

	cd ~/Workshop/GBS_Data/Stacks_Output/
	plink --file batch_1.plink --cluster --mds-plot 3 --allow-extra-chr
	
	cd ~/Workshop/GBS_Data/Corrected_Output/
	plink --file batch_1.plink --cluster --mds-plot 3 --allow-extra-chr
	
Then, you can use your favorite software to plot axes C1 versus C2 from `plink.mds`. Do the samples cluster in any interesting way? Does it change before versus after filtering?

Finally, you can add `plink` to your "master" script. I often make one script from `ustacks` to `populations` (Steps 2-8 above) plus `plink`, but it is very common to encounter a problem or error along the way which can waste time and cause confusion, so an alternative is to have a folder with a separate script for each step (as we started doing yesterday). For this workshop, though, it will be satisfying to see everything you learned in one script. So, make any remaining adjustments to your script and add `plink`.  

Your amended "master" script should now look something like this (there is no need to run it, since we have already run all the steps separately):

	#!/bin/bash
	
	#MASTER SCRIPT FOR STACKS. Assumes process_radtags was already run.
	
	samples="AM2
	BOL3
	CA7
	CR9
	DIN1
	ERO2
	Htwo4
	MA3
	MM6
	SC6
	SNR6
	SR3
	TAT9
	TER1
	TLA1
	ZA7"
	
	#2. RUN USTACKS
	i=1
	for sample in $samples
	do
	   ustacks -t gzfastq -f ~/Workshop/GBS_Data/${sample}.fq.gz -o ~/Workshop/GBS_Data/Stacks_Output/ -i $i -m 3 -M 2 -p 16 -r -d --max_locus_stacks 3 --model_type bounded --bound_high 0.05
	   let "i+=1";
	done
	
	#3. RUN CSTACKS
	samp="" 
	for sample in $samples 
	do 
		samp+="-s ~/Workshop/GBS_Data/Stacks_Output/${sample} "; 
	done

	cstacks -b 1 -o ~/Workshop/GBS_Data/Stacks_Output -p 16 -n 2 $samp
	
	#4. RUN SSTACKS
	for sample in $samples
	do
	   sstacks -b 1 -c ~/Workshop/GBS_Data/Stacks_Output/batch_1 -s ~/Workshop/GBS_Data/Stacks_Output/${sample} -o ~/Workshop/GBS_Data/Stacks_Output -p 16 &>> ~/Workshop/GBS_Data/Stacks_Output/Log
	done
	
	#5. RUN RXSTACKS
	rxstacks -b 1 -P ~/Workshop/GBS_Data/Stacks_Output/ -o ~/Workshop/GBS_Data/Corrected_Output/ --conf_filter --conf_lim 0.25 --prune_haplo --model_type bounded --bound_high 0.05 --lnl_filter --lnl_lim -50.0 -t 16 &>> ~/Workshop/GBS_Data/Corrected_Output/Log
	
	#6. RUN CSTACKS
	samp="" 
	for sample in $samples 
	do 
		samp+="-s ~/Workshop/GBS_Data/Corrected_Output/${sample} "; 
	done

	cstacks -b 1 -o ~/Workshop/GBS_Data/Corrected_Output -p 16 -n 2 $samp
	
	#7. RUN SSTACKS
	for sample in $samples
	do
	   sstacks -b 1 -c ~/Workshop/GBS_Data/Corrected_Output/batch_1 -s ~/Workshop/GBS_Data/Corrected_Output/${sample} -o ~/Workshop/GBS_Data/Corrected_Output -p 16 &>> ~/Workshop/GBS_Data/Corrected_Output/Log
	done
	
	#8. RUN POPULATIONS
	populations -P ~/Workshop/GBS_Data/Corrected_Output/ -M ~/Workshop/GBS_Data/popmap -t 16 -r 1 -p 2 -b 1 -s --lnl_lim -50 --write_single_snp --fasta --plink --structure --vcf
	
	#8A. RUN PLINK
	cd ~/Workshop/GBS_Data/Corrected_Output/
	plink --file batch_1.plink --cluster --mds-plot 3 --allow-extra-chr

That is the whole pipeline! Now let's optimize it.


### Optimizing parameters

There is no easy way to optimize parameter values other than running the program over and over with different combinations of reasonable values. The most important ones to vary are `-m`, `-M`, and `-n` (see [here](http://creskolab.uoregon.edu/stacks/param_tut.php)) and you might also try different values for `--bound_high` and `--max_locus_stacks`. What values are considered reasonable will depend on the biology of your organism(s) and the depth of sequencing. For typical population genetic data, I have tried `-m` 2-5, `-M` 1-3, `-n` 1-3, `--bound_high` 0.01, 0.03, 0.05 and 0.10, and `--max_locus_stacks` 2-4. In my experience, the latter two do not have much impact on the results from my data sets and I set them to 0.05 and 3, respectively. Because it would take far too long to run every combination of every parameter, I usually optimize one at a time holding the other constant. In addition, I set `-M` = `-n` because they are conceptually similar, except that one compares stacks within samples and the other between samples.

To enable running the script repeatedly with different values, we can modify the master script by inserting variables for each parameter value we want to vary. Specifically, we can utilize a particular kind of variable that is automatically generated when arguments are typed into the command line.

For example, we can modify yesterday's example script that we called `ustacks_loop.sh` like this:

	#!/bin/bash
	
	samples="AM2
	SC6
	SR3"
	
	i=1
	for sample in $samples
	do
	   ustacks -t gzfastq -f ~/Workshop/GBS_Data/${sample}.fq.gz -o ~/Workshop/GBS_Data/Stacks_Output -i $i -m $1 -M $2 -p 2
	   let "i+=1";
	done

Do you see the difference? I inserted `$1` for the value of `-m` and `$2` for the value of `-M`. To run the script with `-m 3` and `M 2`, we would now type it in the command line like this (do not run now):

	./ustacks_loop.sh 3 2

So, `$1` is the first command line argument (after the script name) and `$2` is the second, etc. Try modifying the master script to enable command line arguments for `-b` (the batch number should change every time you run Stacks), `-m`, `-M`, `-n`, `--bound_high`, and `--max_locus_stacks`.

If you run the same script multiple times with different parameter values, you will want it to save the results in dfferent folders, so you can also have it create folders with informative names (*e.g.* add `mkdir` to your script).

Try it on your own before reading my version below. 

	#!/bin/bash
	
	#MASTER SCRIPT FOR STACKS. Assumes process_radtags was already run.
	#Order of command line arguments -b -m -M -n --max_locus_stacks --bound_high
	
	#Create output directories and define variables for output folders (so the code is cleaner below)
	mkdir ~/Workshop/GBS_Data/Stacks_Output/b$1_m$2_M$3_n$4_mls$5_bh$6
	OUT=~/Workshop/GBS_Data/Stacks_Output/b$1_m$2_M$3_n$4_mls$5_bh$6
	mkdir ~/Workshop/GBS_Data/Corrected_Output/b$1_m$2_M$3_n$4_mls$5_bh$6
	COR_OUT=~/Workshop/GBS_Data/Corrected_Output/b$1_m$2_M$3_n$4_mls$5_bh$6
	
	#2. RUN USTACKS
	samples="AM2
	BOL3
	CA7
	CR9
	DIN1
	ERO2
	Htwo4
	MA3
	MM6
	SC6
	SNR6
	SR3
	TAT9
	TER1
	TLA1
	ZA7"
	
	i=1
	for sample in $samples
	do
	   ustacks -t gzfastq -f ~/Workshop/GBS_Data/${sample}.fq.gz -o ${OUT}/ -i $i -m $2 -M $3 -p 16 -r -d --max_locus_stacks $5 --model_type bounded --bound_high $6
	   let "i+=1";
	done
	
	#3. RUN CSTACKS
	samp="" 
	for sample in $samples 
	do 
		samp+="-s ${OUT}/${sample} "; 
	done

	cstacks -b $1 -o ${OUT} -p 16 -n $4 $samp
	
	#4. RUN SSTACKS
	for sample in $samples
	do
	   sstacks -b $1 -c ${OUT}/batch_$1 -s ${OUT}/${sample} -o ${OUT} -p 16 &>> ${OUT}/Log
	done
	
	#5. RUN RXSTACKS
	rxstacks -b $1 -P ${OUT}/ -o ${COR_OUT}/ --conf_filter --conf_lim 0.25 --prune_haplo --model_type bounded --bound_high $6 --lnl_filter --lnl_lim -50.0 -t 16 &>> ${COR_OUT}/Log
	
	#6. RUN CSTACKS
	samp="" 
	for sample in $samples 
	do 
		samp+="-s ${COR_OUT}/${sample} "; 
	done

	cstacks -b $1 -o ${COR_OUT} -p 16 -n $4 $samp
	
	#7. RUN SSTACKS
	for sample in $samples
	do
	   sstacks -b $1 -c ${COR_OUT}/batch_$1 -s ${COR_OUT}/${sample} -o ${COR_OUT} -p 16 &>> ${COR_OUT}/Log
	done
	
	#8. RUN POPULATIONS
	populations -P ${COR_OUT}/ -M ~/Workshop/GBS_Data/popmap -t 16 -r 1 -b $1 -s --lnl_lim -50 --write_single_snp --fasta --plink --structure --vcf
	
	#8A. RUN PLINK (MDS and PCA)
	cd ${COR_OUT}/
	plink --file batch_$1.plink --cluster --mds-plot 4 --allow-extra-chr
	plink --file batch_$1.plink --cluster --pca 4 --allow-extra-chr
	
Then, you can run this master script with many sets of command line arguments are compare the outputs (do not run now... it will take too long!). For example, to run it with all the same values we used earlier and yesterday:

	./master_script.sh 1 3 2 2 3 0.05

But if we run this many times, how do we compare all the data!?	Mastretta-Yanes *et al.* (2014) have an excellent way that involves a series of R scripts found [here](https://github.com/AliciaMstt/RAD-error-rates). However, it relies on having replicates of the same individuals in multiple lanes. It is in fact a good idea to do this, but what if you didn't or couldn't?

I use a simpler approach that is in the same spirit of their approach. Instead I choose a small subset of samples from a few sites known (or strongly suspected) to differ genetically (perhaps 3 sample from each of 3 sites), and I run the master script on those with a variety of parameter values. I then choose the set of parameter values that leads to reasonable summary statistics and maximizes SNPs, but most importantly, the one that recovers the expected population structure in a PCA and/or MDS from PLINK. Essentially, we want to minimize the distance between points on the MDS/PCA of individuals from the same sample site while maximizing the difference among sites. If you had identical replicates of individuals their points should be on top of each other (except for error). This logic can be extended to samples from a population under the assumption that they should be more closely related and thus more genetically similar than those for distant populations. 

It will take too long to run an optimization procedure within this workshop, but now you have the script and the tools to do it. Please try it on your own some time and feel free to contact me if you have any questions!
