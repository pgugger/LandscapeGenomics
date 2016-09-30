### The MySQL database

[MySQL](https://www.mysql.com/) is an open-source database system. I only know the basics, just to use it for Stacks. The advantages of using the database are that it can be used to visualize the data in browser window (*e.g.* Firefox), has additional filters in the web interface, and can be used to export a summary file that might be useful for some purposes. It is best to limit how much you run these steps because they are slow and take lots of hard drive space. I usually only run them after the parameters have been optimized.

The computer cluster we have used for this workshop does not currently have the MySQL interface, so I provide these commands only as a guide for when you do have the proper setup. Some of them take a very long time to finish. 

First, we would need to create a database name and structure (it will ask for the password). The name must end in "_radtags".

	mysql -p -e "CREATE DATABASE WorkshopExample_radtags"
	mysql -p WorkshopExample_radtags < /usr/local/share/stacks/sql/stacks.sql

Next, we would load the data into the database. Note that this is only possible if you indicated `-i` for each sample in `sstacks`, `-b` in several steps, and `-s` in populations.

	load_radtags.pl -D WorkshopExample_radtags -p ~/Desktop/GBS_Data/Corrected_Output/ -b 1 -M ~/Desktop/GBS_Data/popmap -c -B -e "Q. rugosa example data from workshop" -t population
	index_radtags.pl -D WorkshopExample_radtags -c -t

Then, you would be able to accessing the database using a web browser by typing into the address bar (or actual host name instead of localhost):

	localhost/stacks
	
You can browse by sample or catalog locus, look at the coverage and evidence for each SNP, and filter the data. You can also export the filtered data from here into an Excel format. If you are satisfied you can generate a file with all the data in a compact form using `export_sql.pl` in the command line. On the [webpage](http://creskolab.uoregon.edu/stacks/comp/export_sql.php) you will see many filters you can apply, but here is the basic command with one filter shown. I have also already run this for you and the results are in `~/Desktop/GBS_Data/MySQL_Export`.

	export_sql.pl -D WorkshopExample_radtags -b 1 -f ~/Desktop/GBS_Data/Corrected_Output/WorkshopExample_alleles.tsv -o tsv -F snps_l=1 -F snps_u=3
	
Finally, note that the databases are big and you need to create a new one for each run that you want to see in the browser. When you are finished with one and no longer want it you can remove the database with the following command.

	mysql -p -e "DROP DATABASE WorkshopExample_radtags"