## EXERCISE 2: Calling SNPs with the Stacks pipeline

[Stacks](http://catchenlab.life.illinois.edu/stacks/) is a sophisticated pipeline for identifying SNPs from RAD-Seq or GBS data. It can be used with or without reference genome, although it is most popular when a reference is not available as in most non-model organisms. There is also a popular [Google Group](https://groups.google.com/forum/#!forum/stacks-users) for asking questions and troubleshooting. For this exercise, we will use the basic *de novo* SNP calling pipeline. We will start with the data we converted to FASTQ in Exercise 1 to demonstrate the quality filtering portion of the pipeline, but then we will switch to some real samples from a Mexican oak tree species, *Quercus rugosa*, for the rest ([Martins *et al.* 2018 *Evol. Appl.*](https://onlinelibrary.wiley.com/doi/full/10.1111/eva.12684)).

### Demultiplexing and removing low quality reads

Stacks provides a convenient tool called [process_radtags](http://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php) for simultaneously demultiplexing and removing low quality reads and adapter contamination. To demultiplex (separate the reads into new files by sample), it requires a file that lists all the barcodes, which I included here as barcodes.txt. First, navigate to `~/Workshop/Test_Data/` and `mkdir Processed_Radtags`.

Now run the `process_radtags` command as follows (all one line).

	process_radtags -p ~/Workshop/Test_Data/ -o ~/Workshop/Test_Data/Processed_Radtags/ -b barcodes.txt -e apeKI --adapter_1 GATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 2 -r -c -q -s 10

Every command in Stacks has many options, so it's best to construct each with the help on the online manual. Let's dissect this command while it's running. With `-p` we indicate the directory with the sequence data, `-o` the output directory, `-b` the barcode file to use for demultiplexing, `-e` the restriction enzyme so it knows what cut site sequence to expect after the barcode, `--adapter_1` to indicate the adapter sequence to identify and remove contaminants, `--adapter_mm` the number of mismatches allowed between the read and real adapter sequence, `-r` to "rescue" (keep) a sequence if there is one error in the barcode, `-c` to remove reads with uncalled bases (N), and `-q` to discard reads with low quality scores which are defined by a sliding window if the average score is less than `-s`. Take a minute to review the `process_radtags` [webpage](http://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php) and see if you understand each part or find additional arguments that could be useful.

We will now switch to using a real data set comprised of GBS samples that were already processed using `process_radtags` from multiple lanes of data.


### The core Stacks pipeline

The core of the Stacks 2 pipeline is
1. `ustacks`
2. `cstacks`
3. `sstacks` 
4. `tsv2bam`
5. `gstacks`
6. `populations`

[ustacks](http://catchenlab.life.illinois.edu/stacks/comp/ustacks.php) is for *de novo* SNP calling. `ustacks` finds sets of exactly matching sequences (*i.e.* "stacks") within samples. [cstacks](http://catchenlab.life.illinois.edu/stacks/comp/cstacks.php) merges similar stacks (alleles) into a set of consensus loci that they call the "catalog", which is a sort or "reference" library of fragments (one per locus). [sstacks](http://catchenlab.life.illinois.edu/stacks/comp/sstacks.php) matches each sample to the catalog. [tsv2bam](http://catchenlab.life.illinois.edu/stacks/comp/tsv2bam.php) rearranges the data for [gstacks](http://catchenlab.life.illinois.edu/stacks/comp/gstacks.php), which calls SNPs and aligns paired-end reads when appropriate. [populations](http://catchenlab.life.illinois.edu/stacks/comp/populations.php) filters SNP data, rearranges it into common file formats, and calculates a variety of population genetic summary statistics. Let's run each step of the pipeline using a few samples of genotyping-by-sequencing (GBS) [Elshire *et al.* 2011 PLoS ONE](https://doi.org/10.1371/journal.pone.0019379) data from *Quercus rugosa* ([Martins *et al.* 2018 *Evol. Appl.*](https://onlinelibrary.wiley.com/doi/full/10.1111/eva.12684)).

We will not go through the reference-based pipeline, but its implementation is simpler and only requires `gstacks` (after directly aligning reads to a reference genome with [bwa](https://sourceforge.net/projects/bio-bwa/files/) or similar to create BAM files) and then `populations`. Let's get started on the *de novo* pipeline.  


### *De novo* SNP calling pipeline
1\. Navigate to `~/Workshop/GBS_Data/` and list the directory contents. You will see 16 compressed FASTQ files (.fq.gz) named by sample ID, the barcodes.txt file you saw earlier, and a popmap file which we will discuss later. Now, create an output folder for the Stacks results: `mkdir Stacks_Output`. Then, we are ready to run `ustacks` for each sample. Try running only the first of the following commands.

	ustacks -t gzfastq -f ~/Workshop/GBS_Data/AM1.fq.gz -o ~/Workshop/GBS_Data/Stacks_Output/ -i 1 -m 3 -M 2 -p 16
	ustacks -t gzfastq -f ~/Workshop/GBS_Data/SC6.fq.gz -o ~/Workshop/GBS_Data/Stacks_Output/ -i 2 -m 3 -M 2 -p 16
	ustacks -t gzfastq -f ~/Workshop/GBS_Data/TAT9.fq.gz -o ~/Workshop/GBS_Data/Stacks_Output/ -i 3 -m 3 -M 2 -p 16
	...

Its progress will be printed to the screen. While its running, look up what each argument does on the [website](http://catchenlab.life.illinois.edu/stacks/comp/ustacks.php) (`-m` and `-M` are especially [important](http://catchenlab.life.illinois.edu/stacks/param_tut.php).). Notice that the only things that change (the only variables) from one line to the next are the input file (*i.e.* the sample) and `-i`. `-i` is an sample number that should be unique for each sample.

Even with only a few samples, you might be noticing how tedious it would to type each command in sequence and then wait for it to finish. Instead, we can write a script to run them all at once using a `for` loop to run the same command for all specified values of any variables. This kind of loop is useful beyond Stacks, in many cases when you would like to execute the same commands on many files. Copy the following example into an empty text file and save it as `ustacks_loop.sh`.

	#!/bin/bash

	samples="AM1
	SC6
	TAT9"
	
	i=1
	for sample in $samples
	do
	   ustacks -t gzfastq -f ~/Workshop/GBS_Data/${sample}.fq.gz -o ~/Workshop/GBS_Data/Stacks_Output/ -i $i -m 3 -M 2 -p 16
	   let "i+=1";
	done

The first line indicates that it is a certain kind of script that runs in the `bash` shell (don't worry about this). Then, it defines the list of samples and assigns them to the variable called `samples`.  I have only listed 3 of 16 samples for this test script. Then, its sets the variable i = 1. Once a variable is defined, you can recall it later with the notation `$variable`, so if you typed i=1 into your Terminal and then typed `echo $i`, it will return the value 1. 

Then, we get to the `for` loop, which says for each `sample` in the variable `samples`, `do` some command (*i.e.* `ustacks ...`). Notice how `$sample` and `$i` recall the variables defined earlier as `sample` and `i`. When run, the first sample name will be inserted in place of `${sample}` and the value for `i` in place of `$i`. Do not be confused by the fact that `samples` and `sample` look similar to us. They are in fact different variables: `samples` is a list of sample names, and `sample` is an item (individual sample) from the list of `samples`. Thus, the `for` loop begins by saying for each individual `sample` in the list of `samples` (defined above the loop), `do`...

After it runs the command inside the loop, it adds 1 to `i` and then goes back to the beginning of the loop. So when the loop runs again, it uses the second sample and `$i` is now 2. So, both `i` and `sample` are changing in parallel. It does this until all the samples in the list have been processed.

To run the script, we first need to make it executable.

	chmod +x ustacks_loop.sh

Then, run it (as is, with only 3 samples). It should take about 15 minutes.

	./ustacks_loop.sh

You just wrote your first shell script and your first loop! This skill will be useful later for Stacks and for bioinformatics in general. When it finishes take a look at the output files using `zcat` and `less`. These files are not very useful at this stage, but it is good to know what each part of the program generates.

2\. For `cstacks`, we run one command that includes all the samples we want to use to build the catalog (~reference). Let's just use the same 3 samples to build the catalog.

	cstacks -P ~/Workshop/GBS_Data/Stacks_Output/ -M ~/Workshop/GBS_Data/popmap.subset -n 2 -p 16

Note that I have indicated the specific samples to be used in constructing the catalog using a file called popmap.subset. Take a look at this file and compare to popmap (use `less -S` to view). These files simply contain the sample name followed by the population that the sample belongs to; popmap has all the samples and popmap.subset has the subset of samples we are experimenting with right now. You can run `cstacks` on any number of samples but it is common to focus on a subset of samples that capture the overall genetic variation (*e.g.* one sample per population) to save computational time. While the command is running [read](http://catchenlab.life.illinois.edu/stacks/comp/cstacks.php) what each argument does (`-n` is especially [important](http://catchenlab.life.illinois.edu/stacks/param_tut.php)). When `cstacks` finishes, look at the output files.
	
3\. `sstacks` runs each sample processed with `ustacks` against the catalog generated with `cstacks`:

	sstacks -P ~/Workshop/GBS_Data/Stacks_Output/ -M ~/Workshop/GBS_Data/popmap.subset -p 16
	
Note that I have specified popmap.subset to to run only the three samples we are testing. After the command finishes, take a look at the output files. So far none of them are easy to read or use in common population genetics software packages.

4\. After `sstacks`, we run `tsv2bam` to convert the Stacks TSV format to the more versatile BAM alignment format, which orients data by locus rather than by sample: 
	
	tsv2bam -P ~/Workshop/GBS_Data/Stacks_Output/ -M ~/Workshop/GBS_Data/popmap.subset -t 16

5\. Final alignment and SNP calling happen with `gstacks`, which you can read about [here](http://catchenlab.life.illinois.edu/stacks/comp/gstacks.php).
	
	gstacks -P ~/Workshop/GBS_Data/Stacks_Output -M ~/Workshop/GBS_Data/popmap.subset -t 16

6\. To filter the outputs of the previous steps, as well as reformat the data for use in other programs, we will run `populations`. It will also calculate summary statistics per population, and to do so requires an additional file that defines the populations (popmap). I have included a file called popmap, which defines all the samples as part of either an eastern or western population. 

	populations -P ~/Workshop/GBS_Data/Stacks_Output/ -M ~/Workshop/GBS_Data/popmap.subset -t 16 -r 0.8 --fasta_loci --plink --genepop --structure --vcf

The output files include a FASTA file of all the loci, a file defining the haplotypes, files in the [PLINK](https://www.cog-genomics.org/plink2) , [STRUCTURE](http://pritchardlab.stanford.edu/structure.html) and [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) (variant call) formats that we requested, and finally two files with summary statistics. Go through each and try to understand what they contain. More details are on the Stacks website.

How many SNPs did you get? Do you think the summary statistics, such as nucletide diversity (pi) and heterozygosity, look reasonable for a tree?


### Additional arguments

The previous steps contain mostly the minimum arguments for each step to run the pipeline. Now, let's explore some of the other arguments that we can add. After this, we will build the beginning of a "master" script that incorporates everything so far.

Open the [webpage](http://catchenlab.life.illinois.edu/stacks/comp/ustacks.php) for `ustacks` and browse through the options. Here are some that I like to include when running `ustacks`:

* `-d` to resolve overmerged stacks
* `--max_locus_stacks` to define the maximum number of stacks that can be merged into a locus
* `--model_type bounded` and `--bound_high` to use a more reasonable model of sequencing error with a maximum rate around 0.01 to 0.10.

Add these terms to your scripts from above (do not run anything yet). For example, the relevant line from `ustacks_loop.sh` might now look something like this:

	ustacks -t gzfastq -f ~/Workshop/GBS_Data/${sample}.fq.gz -o ~/Workshop/GBS_Data/Stacks_Output/ -i $i -m 3 -M 2 -p 16 -d --max_locus_stacks 3 --model_type bounded --bound_high 0.05
	
What effect do you think each argument will have?

Browse the webpages for `cstacks` and `sstacks`. Do you see anything you might add? Personally, I do not use any of the additional options on these parts of Stacks.

Finally, let's look more carefully at the `populations` [manual](http://catchenlab.life.illinois.edu/stacks/comp/populations.php). Here are some potentially important arguments to consider:

* `--write_single_snp` to output only one SNP per fragment. This could be a good idea to minimize the number of closely linked loci that might be in linkage disequilibrium.
* `-r` indicates the proportion of samples for which a a locus must have a genotype call and was included above already. A value of 0.8 is recommended but this could depend on many factors (see [Paris *et al.* 2017](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12775)).
* `-p` to keep only loci found in a at least a certain number of populations. 
* `-B` to specify a file with loci that you do not want to include in the output, for example, if you looked at them and think they are not accurate or if it is a locus that doesn't interest you for some reason.
* `--lnl_lim` to filter low-quality or unlikely loci. Higher log-likelihoods (less negative) are better. Be careful with this setting; I have noticed that setting filters out many good SNPs.

Then there are a few options about how to report summary statistics that I do not typically use. Instead, I use other population genetics software. Finally, you can choose which output file formats you would like, so you can do anaylses in other programs (VCF is useful).

Edit your `populations` command from above by adding `--write_single_snp`, `-p`, `--lnl_lim`, but do not add `-B` today. Do not run yet, but here is how the new `populations` command might look:

	populations -P ~/Workshop/GBS_Data/Stacks_Output/ -M ~/Workshop/GBS_Data/popmap -t 8 -r 0.8 -p 2 --lnl_lim 50 --write_single_snp --fasta_loci --plink --structure --vcf


### The whole *de novo* pipeline in one script
Once you have tested all the parts of the pipeline and added other arguments discussed above, we can build a "master" script with everything together. Try putting yours together before looking at my example below. Name it `master_script.sh`.

	#!/bin/bash
	
	samples="AM1
	BOL3
	CA7
	CR9
	DIN1
	ERO2
	Htwo2
	MA3
	MM6
	SC6
	SNR6
	SR3
	TAT9
	TER1
	TLA1
	ZA7"
	
	#1. RUN USTACKS
		
	i=1
	for sample in $samples
	do
	   ustacks -t gzfastq -f ~/Workshop/GBS_Data/${sample}.fq.gz -o ~/Workshop/GBS_Data/Stacks_Output/ -i $i -m 3 -M 2 -p 16 -d --max_locus_stacks 3 --model_type bounded --bound_high 0.05
	   let "i+=1";
	done
	
	#2. RUN CSTACKS
	cstacks -P ~/Workshop/GBS_Data/Stacks_Output/ -M ~/Workshop/GBS_Data/popmap -n 2 -p 16
	
	#3. RUN SSTACKS
	sstacks -P ~/Workshop/GBS_Data/Stacks_Output/ -M ~/Workshop/GBS_Data/popmap -p 16
	
	#4. RUN TSV2BAM
	tsv2bam -P ~/Workshop/GBS_Data/Stacks_Output/ -M ~/Workshop/GBS_Data/popmap -t 16
	
	#5. RUN GSTACKS
	gstacks -P ~/Workshop/GBS_Data/Stacks_Output -M ~/Workshop/GBS_Data/popmap -t 16
	
	#6. RUN POPULATIONS
	populations -P ~/Workshop/GBS_Data/Stacks_Output/ -M ~/Workshop/GBS_Data/popmap -t 16 -r 0.8 -p 2 --lnl_lim 50 --write_single_snp --fasta_loci --plink --genepop --structure --vcf


**IMPORTANT:** Run your master script with all 16 samples, and make sure it finishes successfully. Before running clear your Stacks_Output folder of the previous outputs (`rm ~/Workshop/GBS_Data/Stacks_Output/*`). It is important to complete this step today because we will use the results tomorrow. However, this script will take a long time (~2 hours), so I recommend letting it run beyond the end of today's workshop, if necessary. To do that, first open a new screen within your terminal by typing `screen`. This screen will remain active even if you disconnect from the server or close your laptop. To get out of ("detach") the screen but not close it, use `CTRL a d`; to return to the same screen, type `screen -r`; to close a screen completely, type `exit` while in the screen. For now, you can just open a screen and run your script in it and worry about the rest tomorrow.

### Optimizing parameters

There is no easy way to optimize parameter values other than running the program over and over with different combinations of reasonable values. The most important ones to vary are `-m`, `-M`, and `-n` (see [here](http://catchenlab.life.illinois.edu/stacks/param_tut.php)) and you might also try different values for `--bound_high` and `--max_locus_stacks`. What values are considered reasonable will depend on the biology of your organism(s) and the depth of sequencing. For typical population genetic data, I have tried `-m` 2-5, `-M` 1-3, `-n` 1-3, `--bound_high` 0.01, 0.03, 0.05 and 0.10, and `--max_locus_stacks` 2-4. In my experience, the latter two do not have much impact on the results from my data sets and I set them to 0.05 and 3, respectively. Because it would take far too long to run every combination of every parameter, I usually optimize one at a time holding the other constant. In addition, I set `-M` = `-n` because they are conceptually similar, except that one compares stacks within samples and the other between samples.

One can run the pipeline over and over manually as we did above, or using the built-in `denovo_map.pl` wrapper script, or by modifying our manual script to facilitate looping. We don't have time to run an optimization procedure within this workshop but if you finish the exercises above early you can try modifying your master script to allow for looping. 

To enable running the script repeatedly with different values, we can modify the master script by inserting variables for each parameter value we want to vary. Specifically, we can utilize a particular kind of variable that is automatically generated when arguments are typed into the command line.

For example, we can modify today's example script that we called `ustacks_loop.sh` like this:

	#!/bin/bash
	
	samples="AM1
	SC6
	SR3"
	
	i=1
	for sample in $samples
	do
	   ustacks -t gzfastq -f ~/Workshop/GBS_Data/${sample}.fq.gz -o ~/Workshop/GBS_Data/Stacks_Output -i $i -m $1 -M $2 -p 16 -d --max_locus_stacks 3 --model_type bounded --bound_high 0.05 --gapped
	   let "i+=1";
	done

Do you see the difference? I inserted `$1` for the value of `-m` and `$2` for the value of `-M`. To run the script with `-m 3` and `M 2`, we would now type it in the command line like this (do not run now):

	./ustacks_loop.sh 3 2

So, `$1` is the first command line argument (after the script name) and `$2` is the second, etc. Try modifying the master script to enable command line arguments for `-m`, `-M`, `-n`, `--bound_high`, and `--max_locus_stacks`.

If you run the same script multiple times with different parameter values, you will want it to save the results in dfferent folders, so you can also have it create folders with informative names (*e.g.* add `mkdir` to your script).

Try it on your own before reading my version below. 

	#!/bin/bash
	
	#MASTER SCRIPT FOR STACKS. Assumes process_radtags was already run.
	#Order of command line arguments -m -M -n --max_locus_stacks --bound_high
	
	#Create output directories and define variables for output folders (so the code is cleaner below)
	mkdir ~/Workshop/GBS_Data/Stacks_Output/m$1_M$2_n$3_mls$4_bh$5
	OUT=~/Workshop/GBS_Data/Stacks_Output/m$1_M$2_n$3_mls$4_bh$5
	
	#1. RUN USTACKS
	samples="AM1
	BOL3
	CA7
	CR9
	DIN1
	ERO2
	Htwo2
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
	   ustacks -t gzfastq -f ~/Workshop/GBS_Data/${sample}.fq.gz -o ${OUT}/ -i $i -m $1 -M $2 -p 16 -d --max_locus_stacks $4 --model_type bounded --bound_high $5
	   let "i+=1";
	done
	
	#2. RUN CSTACKS
	cstacks -P ${OUT}/ -M ~/Workshop/GBS_Data/popmap -n $3 -p 16
	
	#3. RUN SSTACKS
	sstacks -P ${OUT}/ -M ~/Workshop/GBS_Data/popmap -p 16
	
	#4. RUN TSV2BAM
	tsv2bam -P ${OUT}/ -M ~/Workshop/GBS_Data/popmap -t 16
	
	#5. RUN GSTACKS
	gstacks -P ${OUT}/ -M ~/Workshop/GBS_Data/popmap -t 16
	
	#6. RUN POPULATIONS
	populations -P ${OUT}/ -M ~/Workshop/GBS_Data/popmap -t 16 -r 0.8 -p 2 --lnl_lim 50 --write_single_snp --fasta_loci --plink --genepop --structure --vcf
	
	#7. RUN PLINK (MDS and PCA)
	cd ${OUT}/
	plink --file populations.plink --cluster --mds-plot 3 --allow-extra-chr
	plink --file populations.plink --cluster --pca 3 --allow-extra-chr
	
Then, you can run this master script with many sets of command line arguments are compare the outputs (do not run now... it will take too long!). For example, to run it with all the same values we used earlier:

	./master_script.sh 3 2 2 3 0.05

But if we run this many times, how do we compare all the data!? 

[Mastretta-Yanes *et al.* (2014)](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12291) have an excellent way that involves a series of R scripts found [here](https://github.com/AliciaMstt/RAD-error-rates). However, it relies on having replicates of the same individuals in multiple lanes. It is in fact a good idea to do this, but what if you didn't or couldn't?

I use a simpler approach that is in the same spirit of their approach ([Gugger et al. 2018](https://onlinelibrary.wiley.com/doi/full/10.1111/eva.12534)). Instead, I choose a small subset of samples from a few sites known (or strongly suspected) to differ genetically (perhaps 3 samples from each of 3 sites), and I run the master script on those with a variety of parameter values. I then choose the set of parameter values that leads to reasonable summary statistics and maximizes SNPs, but most importantly, the one that recovers the expected population structure in an ordination of the resulting SNP data (*e.g.* PCA and/or MDS run in PLINK). Essentially, we want to minimize the distance between points on the MDS/PCA of individuals from the same sample site while maximizing the difference among sites. If you had identical replicates of individuals their points should be on top of each other (except for error). This logic can be extended to samples from a population under the assumption that they should be more closely related and thus more genetically similar than those for distant populations. 

It will take too long to run an optimization procedure within this workshop, but now you have the script and the tools to do it. Please try it on your own some time and feel free to contact me if you have any questions!

### Other considerations
For the basic pipelines there are a "wrapper" scripts available called [denovo_map.pl](http://catchenlab.life.illinois.edu/stacks/comp/denovo_map.php) and [ref_map.pl](http://catchenlab.life.illinois.edu/stacks/comp/ref_map.php) that allow you to do all the above steps at once. However, I prefer the manual script for more customization and better understanding of the process. In addition, if you have a reference genome sequence you might consider the [GATK](https://www.broadinstitute.org/gatk/) SNP calling pipeline instead.

In Exercise 3, we will analyze the output SNP data for population structure and environmental associations consistent with natural selection.
