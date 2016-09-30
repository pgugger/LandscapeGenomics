## EXERCISE 2: Running the basic Stacks pipeline to call SNPs

[Stacks](http://creskolab.uoregon.edu/stacks/) is a fairly sophisticated and rather complicated pipeline for identifying SNPs from RAD-Seq or GBS data. It can be used with or without reference genome, although it is most popular when a reference is not available as in most non-model organisms. There is also a popular [Google Group](https://groups.google.com/forum/#!forum/stacks-users) for asking questions and troubleshooting. For this exercise, we will start with the basic *de novo* SNP calling pipeline and then you will have an opportunity to try the reference based pipeline. We will start with the data we converted to FASTQ in Exercise 1 to demonstrate the quality filtering portion of the pipeline, but then we will switch to some real samples from *Quercus rugosa* for the rest.

### Demultiplexing and removing low quality reads

Stacks provides a convenient tool called [process_radtags](http://creskolab.uoregon.edu/stacks/comp/process_radtags.php) for simultaneously demultiplexing and removing low quality reads. To demultiplex (separate the reads into new files by sample), it requires a file that lists all the barcodes, which I included here as barcodes.txt. First, navigate to `~/Workshop/QSEQ_Data/` and `mkdir Processed_Radtags`.

Now run the `process_radtags` command as follows (all one line).

	process_radtags -p ~/Workshop/QSEQ_Data/FASTQ/ -i gzfastq -o ~/Workshop/QSEQ_Data/Processed_Radtags/ -y gzfastq -b barcodes.txt -e apeKI --adapter_1 AGATCGGAAGAGCGGTTCAGGAATGCCGAG --adapter_mm 2 -r -c -q -s 10 -t 92

Every command in Stacks has many options, so it's best to construct each with the help on the online manual. Let's dissect this command while it's running. With `-p` we indicate the directory with the sequence data, `-i` the type of sequence file, `-o` the output directory, `-y` the output format, `-b` the barcode file to use for demultiplexing, `-e` the restriction enzyme so it knows what cut site sequence to expect after the barcode, `--adapter_1` to indicate the adapter sequence to identify and remove contaminants, `--adapter_mm` the number of mismatches allowed between the read and real adapter sequence, `-r` to "rescue" (keep) a sequence if there is one error in the barcode, `-c` to remove reads with uncalled bases (N), `-q` to discard reads with low quality scores which are defined by a sliding window if the average score is less than `-s`, and finally `-t` to trim all the reads to the same length of 92 bases. Take a minute to review the `process_radtags` [webpage](http://creskolab.uoregon.edu/stacks/comp/process_radtags.php) and see if you understand each part or find additional arguments that could be useful.

We will now switch to using a real data set comprised of GBS samples that were already processed using the `process-radtags` from multiple lanes of data.


### The core Stacks pipeline

The core of the Stacks pipeline is
1. `ustacks` or `pstacks`
2. `cstacks`
3. `sstacks` 
4. `populations`

[ustacks](http://creskolab.uoregon.edu/stacks/comp/ustacks.php) is for *de novo* SNP calling and [pstacks](http://creskolab.uoregon.edu/stacks/comp/pstacks.php) for reference guided SNP calling. `ustacks` finds sets of exactly matching sequences (*i.e.* "stacks"). [cstacks](http://creskolab.uoregon.edu/stacks/comp/cstacks.php) merges similar stacks (alleles) into a set of consensus loci that they call the "catalog", which is a sort or "reference" library of fragments (one per locus). [sstacks](http://creskolab.uoregon.edu/stacks/comp/sstacks.php) calls SNPs in each sample by matching them to the catalog. [populations](http://creskolab.uoregon.edu/stacks/comp/populations.php) will filter SNP data, rearrange it into common file formats, and calculate a variety of population genetic summary statistics. Let's run each step of the pipeline using a few samples of [genotyping-by-sequencing](http://www.biotech.cornell.edu/brc/genomic-diversity-facility) (GBS) (Elshire *et al.* 2011 PLoS ONE) data from *Quercus rugosa*. 

We will not go through the reference-based pipeline (with `pstacks`), but its implementation would be very similar to the *de novo* pipeline, except that you would first align your samples to the genome sequence (*e.g.*, with bowtie2) and then you would substitute `pstacks` and its arguments where `ustacks` is. Let's get started on the *de novo* pipeline.  


### *De novo* SNP calling pipeline
Navigate to `~/Workshop/GBS_Data/` and list the directory contents. You will see 16 compressed FASTQ files (.fq.gz) named by sample ID, the barcodes.txt file you saw earlier, and a popmap file which we will discuss later. Now, create an output folder for the Stacks results: `mkdir Stacks_Output`. Then, we are ready to run `ustacks` for each sample. Try running the first of the following commands.

	ustacks -t gzfastq -f ~/Workshop/GBS_Data/AM2.fq.gz -o ~/Workshop/GBS_Data/Stacks_Output/ -i 1 -m 3 -M 2 -p 16
	ustacks -t gzfastq -f ~/Workshop/GBS_Data/SC6.fq.gz -o ~/Workshop/GBS_Data/Stacks_Output/ -i 2 -m 3 -M 2 -p 16
	ustacks -t gzfastq -f ~/Workshop/GBS_Data/SR3.fq.gz -o ~/Workshop/GBS_Data/Stacks_Output/ -i 3 -m 3 -M 2 -p 16
	...

Its progress will be printed to the screen. While its running, look up what each argument does on the [website](http://creskolab.uoregon.edu/stacks/comp/ustacks.php) (`-m` and `-M` are especially [important](http://creskolab.uoregon.edu/stacks/param_tut.php).). Notice that the only things that change from one line to the next are the input file (*i.e.* the sample) and `-i`. `-i` is an sample number that should be unique for each sample.

Even with only a few samples, you might be noticing how tedious it would to type each command in sequence and then wait for it to finish. Instead, we can write a script to run them all at once using a `for` loop. Copy the following example into an empty text file and save it as `ustacks_loop.sh`.

	#!/bin/bash

	samples="AM2
	SC6
	SR3"
	
	i=1
	for sample in $samples
	do
	   ustacks -t gzfastq -f ~/Workshop/GBS_Data/${sample}.fq.gz -o ~/Workshop/GBS_Data/Stacks_Output/ -i $i -m 3 -M 2 -p 16
	   let "i+=1";
	done

The first line indicates that it is a certain kind of script that runs in the `bash` shell (don't worry about this). Then, it defines the list of samples and assigns them to the variable called `samples`.  I have only listed 3 of 16 samples for this test script. Then, its sets the variable i = 1. Once a variable is defined, you can recall it later with the notation `$variable`, so if you typed i=1 into your Terminal and then typed `echo $i`, it will return the value 1. 

Then, we get to the `for` loop, which says for each `sample` in the variable `samples`, `do` some command (*i.e.* `ustacks ...`). Notice how `$sample` and `$i` recall the variables defined earlier as `samples` and `i`. When run the first sample name will be inserted in place of `${sample}` and the value for i in place of `$i`. 

After it runs that command it adds 1 to i and then goes back to the beginning of the loop. So when the loop runs again, it uses the second sample and `$i` is now 2. It does this until all the samples in the list have been processed.

To run the script, we first need to make it executable.

	chmod +x ustacks_loop.sh

Then, run it (as is, with only 3 samples). It should take about 15 minutes.

	./ustacks_loop.sh

You just wrote your first shell script and your first loop! This skill will be useful in other parts of Stacks and for bioinformatics in general. When it finishes take a look at the output files using `zcat` and `less`. These files are not very useful at this stage, but it is good to know what each part of the program generates.

For `cstacks`, we run one command that includes all the samples we want to use to build the catalog (~reference). Let's just use the same 3 samples to build the catalog.

	cstacks -b 1 -o ~/Workshop/GBS_Data/Stacks_Output -p 16 -n 2 -s ~/Workshop/GBS_Data/Stacks_Output/AM2 -s ~/Workshop/GBS_Data/Stacks_Output/SC6 -s ~/Workshop/GBS_Data/Stacks_Output/SR3

While the command is running [read](http://creskolab.uoregon.edu/stacks/comp/cstacks.php) what each argument does (`-n` is especially [important](http://creskolab.uoregon.edu/stacks/param_tut.php)). When `cstacks` finishes, look at the output files.

When there are more than just a few samples (as is typically the case), it is helpful to generate the `cstacks` input command itself using a `for` loop. Create a file called `cstacks_loop.sh` with the following code. You can see an example [here](http://creskolab.uoregon.edu/stacks/manual/#phand), which I have edited for our purposes: 
	
	#!/bin/bash

	samples="AM2
	SC6
	SR3"
	
	samp="" 
	for sample in $samples 
	do 
		samp+="-s ~/Workshop/GBS_Data/Stacks_Output/${sample} "; 
	done
	
	cstacks -b 1 -o ~/Workshop/GBS_Data/Stacks_Output -p 16 -n 2 $samp 

The first few lines are the same as our `ustacks` loop and indicate the samples. Then, there is a `for` loop that generates a string of text (in "") for each sample. The string of text is amended (`+=`) to the variable `samp` until the loop finishes all the samples. The string of text is then called with `$samp` in the `cstacks` command. Do you see how this code will essentially generate the exact same command that we ran a moment ago? Although not essential to understanding Stacks, this sort of code is very handy for dealing with large data sets.
	
`sstacks` is the next step and, like `ustacks`, it needs to be run on each sample to call SNPs for each individual relative to the catalog. Write a script with a `for` loop that will run the following example command on the 3 samples we ran with `ustacks` and `cstacks`. Call the script `sstacks_loop.sh` (remember to make it executable). Consult the Stacks manual for help [here](http://creskolab.uoregon.edu/stacks/manual/#phand).

	sstacks -b 1 -c ~/Workshop/GBS_Data/Stacks_Output/batch_1 -s ~/Workshop/GBS_Data/Stacks_Output/AM2 -o ~/Workshop/GBS_Data/Stacks_Output -p 16

When you have created your script and run it, take a look at the output files. So far none of them are easy to read or use in common population genetics software packages.

To filter the outputs of the previous steps, as well as reformat the data for use in other programs, we will run `populations`. It will also calculate summary statistics per population, and to do so requires an additional file that defines the populations. I have included a file called popmap, which defines all the samples as part of either an eastern or western population, but our three test samples are all from the west. 

Take a look at popmap. Then run

	populations -P ~/Workshop/GBS_Data/Stacks_Output/ -M ~/Workshop/GBS_Data/popmap -t 16 -r 1 -b 1 --fasta --plink --structure --vcf

The output files include a FASTA file of all the loci, a file defining the haplotypes, files in the [PLINK](https://www.cog-genomics.org/plink2) , [STRUCTURE](http://pritchardlab.stanford.edu/structure.html) and [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) (variant call) formats that we requested, and finally two files with summary statistics. Go through each and try to understand what they contain. More details are on the Stacks website.

Do you think the summary statistics, such as nucletide diversity (pi) and heterozygosity, look reasonable for a tree?


### Additional arguments

The previous steps contain mostly the minimum arguments for each step to run the pipeline. Now, let's explore some of the other arguments that we can add. After this we will build the beginning of a "master" script that incorporates everything so far.

Open the [webpage](http://creskolab.uoregon.edu/stacks/comp/ustacks.php) for `ustacks` and browse through the options. Here are some that I like to include when running `ustacks`:

* `--gapped` to enable gapped alignments
* `-r` to remove highly repetetive stacks
* `-d` to resolve overmerged stacks
* `--max_locus_stacks` to define the maximum number of stacks that can be merged into a locus
* `--model_type bounded` and `--bound_high` to use a more reasonable model of sequencing error with a maximum rate around 0.01 to 0.10.

Add these terms to your scripts from above (do not run anything yet). For example, the relevant line from `ustacks_loop.sh` might now look something like this:

	ustacks -t gzfastq -f ~/Workshop/GBS_Data/${sample}.fq.gz -o ~/Workshop/GBS_Data/Stacks_Output/ -i $i -m 3 -M 2 -p 16 -r -d --max_locus_stacks 3 --model_type bounded --bound_high 0.05 --gapped
	
What effect do you think each argument will have?

Browse the webpages for `cstacks` and `sstacks`. Do you see anything you might add? Personally, I don't use any of the additional options on these parts of Stacks.

Finally, let's look more carefully at the `populations` [manual](http://creskolab.uoregon.edu/stacks/comp/populations.php). Here are some potentially important arguments to consider:

* `--write_single_snp` to output only one SNP per fragment. This could be a good idea to minimize the number of loci in linkage disequilibrium.
* `-s` to print a file for uploading to a database with a graphical interface. We will not explore this during this workshop, but it can be useful while exploring initial results.
* `-B` to specify a file with loci that you do not want to include in the output, for example, if you looked at them and think they are not accurate or if it is a locus that doesn't interest you for some reason.
* `-p` to keep only loci found in a at least a certain number of populations. It is recommended to use `-p` and `-r` in combination instead of `-m`. 
* `--lnl_lim` to filter low-quality or unlikely loci. Higher log-likelihoods (less negative) are better.

Then there are a few options about how to report summary statistics that I do not typically use. Instead I use other population genetics software. Finally, you can choose which output file formats you would like, so you can do anaylses in other programs (VCF is useful).

Edit your `populations` from above command by adding `--write_single_snp`, `-p`, `--lnl_lim`, but do not add `-s` or `-B` today. Do not run yet, but here is how the new `populations` command might look:

	populations -P ~/Workshop/GBS_Data/Stacks_Output/ -M ~/Workshop/GBS_Data/popmap -t 16 -r 1 -p 2 -b 1 --lnl_lim -50 --write_single_snp --fasta --plink --structure --vcf


### The whole basic *de novo* pipeline in one script
Once you have tested all the parts of the pipeline and added other arguments discussed above, we can build a "master" script with everything together. Try putting yours together before looking at my example below. Name it `master_script.sh`.

	#!/bin/bash
	
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
	
	#1. RUN USTACKS
		
	i=1
	for sample in $samples
	do
	   ustacks -t gzfastq -f ~/Workshop/GBS_Data/${sample}.fq.gz -o ~/Workshop/GBS_Data/Stacks_Output/ -i $i -m 3 -M 2 -p 16 -r -d --max_locus_stacks 3 --model_type bounded --bound_high 0.05 --gapped
	   let "i+=1";
	done
	
	#2. RUN CSTACKS
	samp="" 
	for sample in $samples 
	do 
		samp+="-s ~/Workshop/GBS_Data/Stacks_Output/${sample} ";
	done
	
	cstacks -b 1 -o ~/Workshop/GBS_Data/Stacks_Output -p 16 -n 2 $samp
	
	#3. RUN SSTACKS
	for sample in $samples
	do
	   sstacks -b 1 -c ~/Workshop/GBS_Data/Stacks_Output/batch_1 -s ~/Workshop/GBS_Data/Stacks_Output/${sample} -o ~/Workshop/GBS_Data/Stacks_Output -p 16
	done
	
	#4. RUN POPULATIONS
	populations -P ~/Workshop/GBS_Data/Stacks_Output/ -M ~/Workshop/GBS_Data/popmap -t 16 -r 1 -p 2 -b 1 --lnl_lim -50 --write_single_snp --fasta --plink --structure --vcf


**IMPORTANT:** Run your script with all 16 samples, and make sure it finishes successfully. It is important to complete this step today because we will use the results tomorrow. However, this script will take a long time (~2 hours), so I recommend letting it run beyond the end of today's workshop, if necessary. To do that, first open a new screen within your terminal by typing `screen`. This screen will remain activeeven if you disconnect from the server or close your laptop. To get out of ("detach") the screen but not close it, use `CTRL a d`; to return to the same screen, type `screen -r`; to close a screen completely, type `exit` while in the screen. For now, you can just open a screen and run your script in it and worry about the rest tomorrow.


### Other considerations
For the basic pipelines there are a "wrapper" scripts available called [denovo_map.pl](http://creskolab.uoregon.edu/stacks/comp/denovo_map.php) and [ref_map.pl](http://creskolab.uoregon.edu/stacks/comp/ref_map.php) that allow you to do all the above steps at once. However, I do not recommend them because there are several other steps worth including and customizing that you will learn tomorrow. In addition, if you have a reference genome sequence you might be better off using the [GATK](https://www.broadinstitute.org/gatk/) SNP calling pipeline, which is a bit more rigorous.

In Exercise 3, we will build on the core script to add steps that enhance the quality of the final SNPs.