## EXERCISE 1: Preparing sequence data files

### Viewing sequence files

When you get your Illumina data, it will either come in QSEQ or FASTQ files (the latter is becoming more common). Let's assume you just got your results in QSEQ format.

Log in to the cluster with `ssh` using the Terminal on a Mac or PuTTY on Windows. Type the following and enter your password

	ssh username@132.247.186.44
	
Now copy all the files for the workshop by typing

	cp -r ./ /home/pgugger/Workshop
	
Navigate to the folder on the Desktop called "QSEQ_Data", and list the files in the folder. Remember `cd` and `ls`. Typically, Illumina data are given in a series of files, usually in multiples of 48 or 96, that are compressed (.gz). I have given you just two.

We can look at a compressed file without extracting it:

	zcat s_3_1_1101_qseq.txt.gz | less -S

The `-S` flag indicates not to wrap lines when longer than the screen is wide. If you scroll down you can see the file is very long. Press `q` to quit.

We can find out how many lines the file has, and therefore how many reads it has, using `wc`, which counts the number of lines, words, and characters by default. 

	zcat s_3_1_1101_qseq.txt.gz | wc -l


### Converting from QSEQ to FASTQ

To convert formats we will use a PERL script I downloaded that requires the files be extracted first. Let's extract both at the same time

	gunzip *.gz

List the files to see how big they are. If you don't know how to show the file sizes, look at the man page for `ls` by typing `man ls` to learn about the options.

	ls -lh

Now make a directory to put the FASTQ files in: `mkdir FASTQ`. Finally, we are ready to convert from QSEQ to FASTQ. Note that for some programs it is useful to combine all the files into one large data file (`cat *.txt > alldata.qseq`), but we don't need to do that today.

	qseq2fastq.pl s_3_1_1101_qseq.txt FASTQ/s_3_1_1101.fq

This will take a while. In the meantime, you can open another Terminal and run the other file too.

When it finishes, navigate to the FASTQ folder and view one of the FASTQ files with `less`.


### Assessing quality with FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc) is a convenient program for summarizing base call quality scores, GC content, possible contamination, and other information. Its easy to run but can be slow.

	fastqc s_3_1_1101.fq

When the program finishes, open the "QSEQ_Data" folder by double-clicking with the mouse. Then open the HTML file by double-clicking. How do the Phred scores look at each base?


### Assessing possible contaminants with FastQ_Screen

Let's say you want to know how much contamination is in your data. You might try to search (`grep`) the FASTQ for part of the adapter sequence and count (`-c`) how many sequences have it.

	grep -c 'AGATCGGAAGAGCGGT' s_3_1_1101.fq

However, what if there are errors, so you don't have an exact match to the adapter, and what if there are other contaminants like your own DNA or another organism? One way to get a rough estimate is to use [FastQ_Screen](http://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqscreen), which will align your reads to specified genomes or sequences).

First, we need to set up the confirguration file. Open fastq_screen.conf with `vi` or using your FTP client (*e.g.*, WinSCP). The top section indicates the path to the [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) aligner for aligning short reads to genomes. The Threads section indicates how many processors your computer has available to run the analysis. The last section, called Database, give the path to different genomes or sequences that you want to align your data to. Note that the genomes/sequences must be indexed first (not FASTA) using Bowtie2 commands, which I have already done for you and placed in `~/Workshop/Genomes/`. Edit the path names to reflect your working directory (*e.g.*, in the Database section, change `/home/pgugger/*` to `/home/your_username/*`). To edit in `vi`, type `i` for "insert" then use the keyboard to edit; press `ESC` to stop editing; and then save and close the file by typing `:x` (or close without saving using `:q!`). You can now run the program.

	fastq_screen --aligner bowtie2 --conf ~/Workshop/fastq_screen.conf s_3_1_1101.fq

When the program finishes, open the "QSEQ_Data" folder by double-clicking with the mouse. Then open the PNG file. What fraction of reads align to the specified genomes? How many don't align anywhere? What could that mean!?

Finally, compress the FASTQ file to save space.

	gzip s_3_1_1101.fq
	
In Exercise 2, you will learn how to process the data for SNP calling in Stacks.
