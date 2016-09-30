## EXERCISE 4: Environmental association analysis

In this exercise, we will learn the basics of gradient forest (GF) and generalized dissimilarity modeling (GDM), which allow us to quantify, describe, and map the association spatial and environmental factors with genetic variation. Refer to Fitzpatrick & Keller (2015 *Ecology Letters*) for more detail about concepts and application, and refer to the R package manuals and vignettes for details and examples of syntax and usage (*e.g.*, [GDM](https://cran.r-project.org/web/packages/gdm/vignettes/gdmVignette.pdf), [GF](http://gradientforest.r-forge.r-project.org/biodiversity-survey.pdf)). Both are multivariate approaches that consider potential non-linear patterns in the data. GF is a machine learning approach, whereas GDM is a distance-based approach. Other similar approaches that you might consider are [redundancy analysis (RDA)](https://github.com/pgugger/LandscapeGenomicsWorkshop_Morelia/blob/master/Exercise4.mdown), which is a constrained ordination method that assumes linear relationships among variables, or the R package [BEDASSLE](http://cran.r-project.org/web/packages/BEDASSLE/index.html), which is distance-based. 

Typically, GF and GDM are used with allele frequency data based on sample sites ("populations"), but we will use individual-based data and treat each of the 16 samples as data points. Otherwise, we wouldn't have enough data points to run the models.

### Preparing the data
#### SNP data
First, we need to take the SNP output from Stacks and convert it to a convenient format for R. First, we will use `vcftools`, which is generally very useful for manipulating and filtering VCFs, to output the SNPs as genotypes encoded 0, 1, or 2 for homozygous, heterozygous, and other homozygous, respectively. 

	cd ~/Workshop/GBS_Data/Corrected_Output
	vcftools --vcf batch_1.vcf --012 --out snp

Take a quick look at the output (file names with "012") to understand what each one represents. Then, the easiest thing to do would be to open the files and copy, paste, and edit them the way you want in a text editor or Excel. To stay in the command line environment and acheive the same goal, we can instead run the following commands: 

	cut -f2- snp.012 | sed 's/-1/NA/g' >snp.temp
	tr -d '\t' <snp.012.pos | tr '\n' '\t' | sed 's/[[:space:]]*$//' >header
	paste <(echo "ID" | cat - snp.012.indv) <(echo "" | cat header - snp.temp) > snp.forR
	rm header snp.temp

Don't worry about the details of these commands, but someday you may find it useful to learn them. To briefly summarize: we `cut` to keep from the second to last column of snp.012, then replaced (`sed`) the missing data encoded as -1 to encoded as NA (assuming we allowed missing data); then rearranged the list of SNP positions into a single tab-separated row to create a header; then `paste`d the SNP IDs with the header and SNP data. View the resulting file with `less -S`. When finished copy it to the `GF_GDM` folder with `cp snp.forR ~/Workshop/GF_GDM/`.

Now let's move to R by simply typing `R` in the command line. You are now in the command line environment of R, not the Linux itself. Some of you may be able to connect with the graphical interface of RStudio, but in either case the commands will be the same.
	
Set the directory you will be working in.

	setwd("~/Workshop/GF_GDM")

Load the SNP data

	snp <- read.table("snp.forR", header = T, row.names = 1)

#### Climate data
We will use WorldClim climate data, specifically temperature seasonality (bio4), minimum temperature of coldest month (bio6), precipitation seasonality (bio15), and precipitation of driest quarter (bio17). I put these climate data layers in the folder `bio_22_tif` as GeoTiff format. First, we will load the required packages and the data.
	
	library(raster)
	library(rgdal)
	clim.list <- dir("./bio_22_tif/", full.names=T, pattern='.tif')  #makes list of file paths for each layer
	clim.layer <-  stack(clim.list)  #stacks the layers into a single object

Now, we will crop the climate data layers to just the area of interest. I defined an extent with min and max longitude and min and max latitude, allowing us to then use the `crop` function.

	extent <- c(-104, -96, 18, 22) 
	clim.layer.crop <- crop(clim.layer, extent)

We can compare before and after by simply plotting the climate data.

	pdf("clim.layer.BandA.pdf")
	plot(clim.layer)
	plot(clim.layer.crop)
	dev.off()

These commands save the plots as a PDF. If you are logged into the computer cluster, you will need to retrieve the PDF with `scp` (*e.g.*, WinSCP) to view it on your laptop. If you are running R locally on your laptop, then you can just run the `plot` commands without the `pdf` and `dev.off` commands, which will cause the plots to appear in a window. You can do this for all future plots in this exercise, but you will have to manually save them if you want to keep them.

Finally, we need to load our sample coordinates and extract climate data at each of the points.

	sample.coord <-read.table("sample.coord.txt", header=T, stringsAsFactors=F)
	sample.coord
	crs.wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"  #defines the spatial projection system that the points are in (usually WGS84)
	sample.coord.sp <- SpatialPointsDataFrame(sample.coord[,c('Longitude','Latitude')], proj4string=CRS(crs.wgs), data=sample.coord)

	clim.points <- extract(clim.layer.crop, sample.coord.sp)  #extracts the data for each point (projection of climate layer and coordinates must match)
	clim.points <- cbind(sample.coord, clim.points)  #combines the sample coordinates with the climate data points
	write.table(clim.points, "clim.points", sep="\t", quote=F, row.names=F)  #save the table for later use
	clim.points


### Gradient forest (GF)

We will model the associations of spatial and climate variables with allele frequencies (genotypes) of individuals. You have already seen how we got the climate and genetic data. For the spatial variables, one could use latitude and longitude, but a more sophisticated approach might be to use PCNMs or MEMs (principal coordinates of neighbor matrices or Moran's eigenvector maps). These approaches generate a set of uncorrelated spatial variables. We do not have time to go into these methods in detail, but here is the code to generate the PCNM spatial variables.

	library(vegan)
	coord <- clim.points[,c("Longitude","Latitude")]
	pcnm <- pcnm(dist(coord))  #this generates the PCNMs, could stop here if you want all of them
	keep <- round(length(which(pcnm$value > 0))/2)
	pcnm.keep <- scores(pcnm)[,1:keep]  #keep half of positive ones as suggested by some authors
	pcnm.keep

Now we are ready to begin GF. Load the library and then create a file that contains only the climate and spatial variables (no lat/lon).

	library(gradientForest)
	env.gf <- cbind(clim.points[,c("bio15_22", "bio17_22", "bio4_22", "bio6_22")], pcnm.keep)

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

Do you see any interesting patterns? It appears genetic variation changes abruptly for bio15 values of 60 and then acheives high cumulative importance. Several other plots and tables can be output, including results for each SNP. I recommend generating these with your own data sets, but we will skip them today.

One of the useful features of these models is that the results can be readily projected spatially across a landscape that was not necessarily completely sampled. The mapping can be a little tricky, but first we must extract the climate data for **all** the cells in the climate data layer and remove missing data.

	clim.land <- extract(clim.layer.crop, 1:ncell(clim.layer.crop), df = TRUE)
	clim.land <- na.omit(clim.trns)

Then, we can use the `predict` function to "predict" the allelic turnover across the whole landscape using our fittend model (`gf`) and the climate values on the landscape (`clim.land`).

	pred <- predict(gf, clim.land[,-1])  #note the removal of the cell ID column with [,-1])

These predictions then need to be converted to a color scale for mapping. One way is to use principal components analysis (PCA) on the predictions and use the first three axes to define red, green, blue color scales, respectively. After the values are defined in color space, they can be stacked and mapped. 

	PCs <- prcomp(pred, center=T, scale.=F)
	r <- PCs$x[, 1]
	g <- PCs$x[, 2]
	b <- PCs$x[, 3]
	r <- (r - min(r))/(max(r) - min(r)) * 255
	g <- (g - min(g))/(max(g) - min(g)) * 255
	b <- (b - min(b))/(max(b) - min(b)) * 255
	mask<-clim.layer.crop$bio15_22
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

The colors represent genetic variation as predicted based on the modeled relationships. Similar colors are more similar genetically.

### Generalized dissimilarity modeling (GDM)

GDM is a distance-based method, so we will first need to generate a genetic distance measure among individuals. We will simply use Euclidean distance, but if you had populations as sample points you could use pairwise *F*st or other widely used measures in genetics. We first estimate genetic (Euclidean) distance, then rescale to values between 0 and 1, and finally manipulate the result so it has the proper column with sample IDs.

	snp.dist <- dist(snp, diag = T, upper = T)  #generate Euclidean distance matrix
	snp.dist.1 <- snp.dist/(max(snp.dist))  #rescale by dividing by max value
	snp.dist.1 <- as.matrix(snp.dist.1)
	snp.dist.1 <- cbind(ID = rownames(snp.dist.1), snp.dist.1)  #make the row names an actual column named "ID"
	rownames(snp.dist.1) <- NULL  #remove prior row names
	write.table(snp.dist.1, "snp.dist.1", sep="\t", quote=F, row.names=F)  #save a copy if you want
	snp.dist.1  #view the final distance matrix

Now, we are ready for GDM. Load the library and format the input.

	library(gdm)
	gdm.input <- formatsitepair(bioData=snp.dist.1, bioFormat=3, predData=clim.points, siteColumn="ID", XColumn="Longitude", YColumn="Latitude")

There are several ways to tell GDM the structure of your input data using the `formatsitepair` function, and I refer you to the manual for details. Ours are format "Type 3", our predictor data (climate data points, etc.) are in `clim.points`, and in that object, the samples are given in the column named "ID" (`siteColumn`), *x* (`XColumn`) is in the column "Longitude", and *y* (`YColumn`) is in the column "Latitude". Latitude and longitude will be used to calculate spatial distances.

Once the input is formatted, we can run GDM indicating the input (`gdm.input`) and that we want it to consider spatial distance (`geo = T`).
	
	gdm <- gdm(gdm.input, geo = T, splines = NULL, knots = NULL)

To summarize the results and view an estimate of the deviance explained, type

	summary.gdm(gdm)
	gdm$explained
	
How much deviance is explained by the model? What does that mean?

To better undertand the importance of each variable in the overall model, we can use the `gdm.varImp` function to assess variable importance. Let's try that and then plot the results.

	gdm.importance <- gdm.varImp(gdm.input, geo=T, nPerm=100, parallel=TRUE, cores=16)
	
	pdf("GDM_VariableImportance.pdf")
	barplot(sort(gdm.importance[[2]][,1], decreasing=T))
	dev.off()

Similar to GF, we can get "turnover functions" showing the relationship of genetic composition with geographic and climate gradients.

	pdf("GDM_TurnoverFunctions.pdf")
	plot(gdm, plot.layout = c(3, 3))
	dev.off()

What do these plots suggest? How do the variable importance and turnver functions compare to results from GF?

Finally, we can map the results by transforming the climate data layers (`clim.layer.crop`) based on the modeled results (`gdm`), and using PCA to covert to a color scale.

	clim.trans <- gdm.transform(gdm, clim.layer.crop)
	clim.rast <- na.omit(getValues(clim.trans))
	
	pca <- prcomp(clim.rast)
	pca.rast <- predict(clim.trans, pca, index=1:3)
	pca.rast[[1]] <- (pca.rast[[1]]-pca.rast[[1]]@data@min) / (pca.rast[[1]]@data@max-pca.rast[[1]]@data@min)*255
	pca.rast[[2]] <- (pca.rast[[2]]-pca.rast[[2]]@data@min) / (pca.rast[[2]]@data@max-pca.rast[[2]]@data@min)*255
	pca.rast[[3]] <- (pca.rast[[3]]-pca.rast[[3]]@data@min) / (pca.rast[[3]]@data@max-pca.rast[[3]]@data@min)*255

	pdf("GDM_Map.pdf")
	plotRGB(pca.rast, r=1, g=2, b=3, bgalpha=0)
	points(clim.points$Longitude,clim.points$Latitude)
	dev.off()
	
Compare to the GF results. Overall, the methods appear to suggest that different predictors are important. However, we should be careful interpreting too much due to low sample size. There are additional features that can be used to generate confidence intervals and further explore the data. We will stop here, and in the next exercise we will see one way to detect specific SNPs that have exceptionally strong associations with the environment.