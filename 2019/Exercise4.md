## EXERCISE 4: Landscape genomics - quantifying and mapping patterns with gradient forest

In this exercise, we will learn the basics of [gradient forest (GF)](https://doi.org/10.1890/11-0252.1) to quantify, describe, and map the associations of spatial and environmental variables with genetic variation. Refer to [Fitzpatrick & Keller (2015)](http://dx.doi.org/10.1111/ele.12376) for more detail about concepts and application of GF, and refer to the R package manuals and [vignettes](http://gradientforest.r-forge.r-project.org/biodiversity-survey.pdf) for details and examples of syntax and usage. 

GF is a multivariate, machine learning approach that considers potential non-linear patterns in the data. GF partitions the allele frequency data at split values along the environmental gradients. Split importance, a measure of the amount of variation explained, is high in positions along the gradient where allelic change is large. Moving along the gradient, the split importance values are summed cumulatively to produce a step-like function for allele frequency change along the environmental gradient, thus quantifying and describing the shape of the relationship between genotypic and environmental data. When run on a large "random" SNP data set, we can infer the contribution of isolation by environment or how environment shapes neutral evolutionary processes such as drift and gene flow. When run on adaptive SNP variation, we can learn more about how selection shapes variation across natural landscapes. 

Other approaches that you might consider are generalized dissimilarity modeling (GDM), which is distance-based multivariate approach that also explores nonlinear relationships, or redundancy analysis (RDA), which is a constrained ordination method that assumes linear relationships among variables.

In this exercise, we will run two GF models and compare the results. The first will contain all the SNP data, whereas the second will only contain the putatively adaptive SNPs that we identified yesterday with LFMM. We will quantify the importance of different environmental and spatial predictors and plot how each class of genetic variation changes along the environmental gradients and landscape. 

Typically, GF is used with allele frequency data based on sample sites ("populations"), but we will use individual-based data and treat each of the 16 samples as data points. Otherwise, we wouldn't have enough data points to run the models.

### Preparing the GF input

#### SNP data

GF can read SNP data in a similar 012 format to the one we generated for LFMM. However, we need to change the missing data symbol from "9" to "NA", and it would be nice if we also included sample IDs as row names and locus/SNP IDs as column names to improve clarity. Here is one way to do that in the Linux command line, starting from the original 012-formatted files that we generated with `vcftools`:

	cd ~/Workshop/GBS_Data/Stacks_Output
	cut -f2- snp.012 | sed 's/-1/NA/g' >snp.temp
	tr '\t' '_' <snp.012.pos | tr '\n' '\t' | sed 's/[[:space:]]*$//' >header
	paste <(echo "ID" | cat - snp.012.indv) <(echo "" | cat header - snp.temp) > snp.forR
	rm header snp.temp
	less -S snp.forR
	mv snp.forR ~/Workshop/GF

Don't worry about the details of these commands, but someday you may find it useful to learn them. To briefly summarize: we `cut` to keep from the second column to last column of snp.012, then replaced (`sed`) the missing data encoded as -1 to be encoded as NA; then rearranged the list of SNP positions into a single tab-separated row to create a header of SNP IDs as locus_position; then `paste`d the SNP IDs with the header and SNP data. 

Now that we have nice table of SNP data, we can enter `R`, load the SNP table, and then extract the subset of adaptive loci that we want to model separately and save it as another object.

Set the directory you will be working in.

	setwd("~/Workshop/GF")

Load the SNP data and the LFMM results from yesterday.

	snp <- read.table("snp.forR", header = T, row.names = 1)
	lfmm.results <- read.table("~/Workshop/LFMM/lfmm.results", header = T)

First, we can see examples of which rows represent significant LFMM associations:

	#Just Pdry
	which(lfmm.results$q.pdry < 0.01 & abs(lfmm.results$z.pdry) > 2, )  
	
	#Any climate variable
	which( (lfmm.results$q.pdry < 0.01 & abs(lfmm.results$z.pdry) > 2 ) | (lfmm.results$q.pseas < 0.01 & abs(lfmm.results$z.pseas) > 2) | (lfmm.results$q.tmin < 0.01 & abs(lfmm.results$z.tmin) > 2) | (lfmm.results$q.tseas < 0.01 & abs(lfmm.results$z.tseas) > 2), )  

Now that we know which rows, we can select the corresponding columns in the `snp` data table. Recall that the order of the rows in `lfmm.results` is the same as the order of the columns in `snp.forR`. I will proceed with all the putatively significant SNPs, starting by saving the list of row numbers as a new object:

	lfmm.sig <- which( (lfmm.results$q.pdry < 0.01 & abs(lfmm.results$z.pdry) > 2 ) | (lfmm.results$q.pseas < 0.01 & abs(lfmm.results$z.pseas) > 2) | (lfmm.results$q.tmin < 0.01 & abs(lfmm.results$z.tmin) > 2) | (lfmm.results$q.tseas < 0.01 & abs(lfmm.results$z.tseas) > 2), ) 
	length(lfmm.sig)
	
	snp.adaptive <- snp[ , lfmm.sig]
	dim(snp.adaptive)
	
Now, we have a table of adaptive SNPs (`snp.adaptive`) in the same format as the table of all SNPs (`snp`). You can write the adaptive SNPs to file if you like.

#### Climate data

We already prepared the climate point data yesterday as `clim.points`, which should already be in your GF folder (refer to Exercise 3).

	clim.points <-read.table("clim.points", header=T)

GF will also utilize the climate layers for projecting onto maps, so recall that these files are in `~/Workshop/Climate_Data`. However, we should crop these to focus only on our study area and save them as a new object. First, load and stack the climate layer data as we did yesterday:

	library(raster)
	library(rgdal)
	clim.list <- dir("~/Workshop/Climate_Data/", full.names=T, pattern='.tif')  
	clim.layer <-  stack(clim.list)  #stacks the layers into a single object

To crop the climate data layers to just the area of interest, I defined an extent with minimum and maximum longitude and minimum and maximum latitude and then use the `crop` function.

	extent <- c(-104, -96, 18, 22) 
	clim.layer.crop <- crop(clim.layer, extent)

We can plot maps of the cropped layers, which you can compare to the complete climate layer maps from yesterday. 

	pdf("clim.layer.crop.pdf")
	plot(clim.layer.crop)
	dev.off()
	
#### Spatial data

For the spatial variables, one could use latitude and longitude, but a more sophisticated approach might be to use PCNMs or MEMs (principal coordinates of neighbor matrices or Moran's eigenvector maps). These approaches generate a set of uncorrelated spatial variables. We do not have time to go into these methods in detail, but here is the code to generate the PCNM spatial variables.

	library(vegan)
	coord <- clim.points[,c("Longitude","Latitude")]
	pcnm <- pcnm(dist(coord))  #this generates the PCNMs, you could stop here if you want all of them
	keep <- round(length(which(pcnm$value > 0))/2)
	pcnm.keep <- scores(pcnm)[,1:keep]  #keep half of positive ones as suggested by some authors
	pcnm.keep


### Running the GF model with all SNPs

Now we are ready to begin GF to model the associations of spatial and climate variables with allele frequencies (genotypes) of individuals. Load the library and then create an object that contains only the climate and PCNM spatial variables (no lat/lon).

	library(gradientForest)
	env.gf <- cbind(clim.points[ , c("Pdry", "Pseas", "Tmin", "Tseas") ], pcnm.keep)

In GF, a maximum number of splits to evaluate can be defined following the developers suggestion

	maxLevel <- log2(0.368*nrow(env.gf)/2)

Running GF involves only one command:

	gf <- gradientForest(cbind(env.gf, snp), predictor.vars=colnames(env.gf), response.vars=colnames(snp), ntree=500, maxLevel=maxLevel, trace=T, corr.threshold=0.50)

The input is the combined climate, spatial, and SNP data as input (`cbind(env.gf, snp)`), and the subsequent parts of the command define which variables are predictors and response variables, as well as a number of other parameters that I have left as suggested. GF will take ~30 minutes to run. In the meantime, read about the various settings or how the model works. When it finishes, there will be warnings about having less than five values for response variables, which is because we have only three: 0, 1, or 2. You can ignore them.

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

Do you see any interesting patterns? It appears that genetic variation changes abruptly for Pseas values of 60 and then acheives high cumulative importance.

We can also make plots of turnover functions for individual loci:

	pdf("GF_TurnoverFunctions_bySNP.pdf")
	plot(gf, plot.type = "C", imp.vars = by.importance, show.overall = F, legend = T, leg.posn = "topleft", leg.nspecies = 5, cex.lab = 0.7, cex.legend = 0.4, cex.axis = 0.6, ylim = c(0, 0.5), line.ylab = 0.9, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0, 0.3, 0, 0)))
	dev.off()

Each line within each panel represents allelic change at a single SNP. Notice that in each panel some SNPs show very steep changes along the environmental gradient. One might consider these SNPs as candidates for involvement in local adaptation along the gradient. This approach to "outlier" detection is still being tested, but Fitzpatrick & Keller (2015) show a promising example. If interested, you can check if any of these highly associated SNPs (*e.g.*, those listed in the plot legend) were also significant in your LFMM analyses.

Several other plots and tables can be output from `gradientForest`. I recommend generating these with your own data sets, but we will skip them today.

One of the useful features of these models is that the results can be readily projected spatially across a landscape that was not necessarily completely sampled. The mapping can be a little tricky, but first we must extract the climate data for *all* the cells in the climate data layer and remove missing data.

	clim.land <- extract(clim.layer.crop, 1:ncell(clim.layer.crop), df = TRUE)
	clim.land <- na.omit(clim.land)

Then, we can use the `predict` function to "predict" the allelic turnover across the whole landscape using our fitted model (`gf`) and the climate values on the landscape (`clim.land`).

	pred <- predict(gf, clim.land[,-1])  #note the removal of the cell ID column with [,-1]

These predictions then need to be converted to a color scale for mapping. One way is to use principal components analysis (PCA) on the predictions and use the first three axes to define red-green-blue color scales, respectively. After the values are defined in color space, they can be stacked and mapped. 

	pca <- prcomp(pred, center=T, scale.=F)
	
	#Assign PCs to colors
	r <- pca$x[, 1]
	g <- pca$x[, 2]
	b <- pca$x[, 3]
	
	#Scale colors
	r <- (r - min(r))/(max(r) - min(r)) * 255
	g <- (g - min(g))/(max(g) - min(g)) * 255
	b <- (b - min(b))/(max(b) - min(b)) * 255
	
	#Define raster properties with existing one
	mask <- clim.layer.crop$Pseas
	
	#Assign color to raster
	rastR <- rastG <- rastB <- mask
	rastR[clim.land$ID] <- r
	rastG[clim.land$ID] <- g
	rastB[clim.land$ID] <- b
	
	#Stack color rasters
	rgb.rast <- stack(rastR, rastG, rastB)

Finally, here is how to make the map, along with the original data points added for reference:

	pdf("GF_Map.pdf")
	plotRGB(rgb.rast, bgalpha=0)
	points(clim.points$Longitude, clim.points$Latitude)
	dev.off()

The colors represent genetic variation (allelic composition) as predicted based on the modeled relationships with environmental and spatial variables. Similar colors are predicted to be more similar genetically.

### Running the GF model with adaptive SNPs

Using the above commands as a template, run the model on only the adaptive SNPs that we saved in the object called `snp.adaptive`. Save your results from the `gradientForest` function as `gf.adaptive`, your results from the `predict` function as `pred.adaptive`, and your results from the `prcomp` function as `pca.adaptive`.

### Comparing overall and adaptive genetic variation

One way to compare the results is to visually review the output plots side by side. Which variables explain the most variation in each SNP data set? Do climate variables tend to explain more variation in the putatively climate-adaptive SNPs? What insights might you draw at this stage?

More careful comparison warrants new plots, especially for turnover functions and maps. For example, the allelic turnover functions of the adaptive SNPs and the complete SNP data set can be displayed on the same plot. We would expect climate-adaptive SNPs to be more strongly associated with environmental gradients than the full SNP set. First, we need to extract the cumulative importance data for predictors of interest and then plot. For now, let's just look at one spatial variable (PCNM1) and one environmental variable (Pdry), but you can feel free to plot more on your own.

	#Extract cumulative importance data from GF results. Note that these commands retrieve all the predictors, despite specifying only one. Also note that you could set type="Species" if you want to extract data for individual loci.
	gf.cumimp <- cumimp(gf, predictor="PCNM1", type="Overall")
	gf.adaptive.cumimp <- cumimp(gf.adaptive, predictor="PCNM1", type="Overall" )
	
	pdf("Comparison_PCNM1.pdf")
	plot(cumimp(gf.adaptive, "PCNM1", standardize=FALSE), main="", ylab="Cumulative importance", xlab='PCNM1', type='l', lwd=2, lty="dashed", col="red")
	lines(cumimp(gf, "PCNM1"),type='l',col="black", lwd=2)
	legend("topleft", box.lty=0, legend=c("All", "Adaptive"), lwd=c(2,2), col=c("black", "red"), lty=c("solid", "dashed"))
	dev.off()

	pdf("Comparison_Pdry.pdf")
	plot(cumimp(gf.adaptive, "Pdry", standardize=FALSE), main="", ylab="Cumulative importance", xlab='Precipitation of dryest quarter', type='l', lwd=2, lty="dashed", col="red")
	lines(cumimp(gf, "Pdry"),type='l',col="black", lwd=2)
	legend("topleft",box.lty=0, legend=c("All", "Adaptive"), lwd=c(2,2), col=c("black","red"), lty=c("solid","dashed"))
	dev.off()
	
How does adaptive variation compare to the background? Is it what you expected, and how would you interpret these results?

A second potentially interesting comparison is to map where adaptive genetic variation deviates most from landscape patterns of background genetic variation. The approach will be to use Procrustes rotation to compare the PCAs that you generated from predictions based on the all-SNP and adaptive-SNP models. Recall the following commands that you should have already run above:

	#Do not run (assuming you successfully ran above)
	pred <- predict(gf, clim.land[,-1])  
	pca <- prcomp(pred.gf.all, center=T, scale.=F)
	
	pred.adaptive <- predict(gf.adaptive, clim.land[,-1])  
	pca.adaptive <- prcomp(pred.adaptive, center=T, scale.=F)

To run the Procrustes rotation, we simply provide the two PCAs for comparison. The method rotates the PCAs to orient the data sets similarly. We can then estimate the differences between them as the residuals

	diff.procrustes = procrustes(pca.adaptive, pca, scale=TRUE, symmetrical=FALSE)
	resids = residuals(diff.procrustes)
	
Finally, we can generate a raster layer based on the residuals, choosing any color scale, and then plot as a map.

	#Specify raster properties based on an existing one
	rastProc <- clim.layer.crop$Pseas  
	
	#Assign residuals to raster cells
	rastProc[clim.land$ID] <- resids  
	
	#Map the results
	pdf("ProcrustesMap.pdf")
	plot(rastProc, col = rev( rainbow( 99, start = 0, end = 0.2 ) ))
	points(clim.points$Longitude, clim.points$Latitude)
	dev.off()

Redder areas (in this case) are predicted to have larger differences in adaptive variation than you would expect based on overall genetic variation. How might you interpret these patterns? Perhaps they provide some indication of the strength of local adaptation in certain regions?

### Mapping future predictions 

Another interesting analysis is to use our GF model results to predict what genetic variation would look like in the future under a future climate scenario, and then calculate the difference between the future and present-day projections. This difference has been called the "genetic offset" and can be an indicator of population vulnerability to climate change when analyzing adaptive genetic variation. The basic premise is that regions with large differences are areas where the climate will be most mismatched with current patterns of local adaptation and thus those populations may be maladapted to future conditions at the site. Therefore, we will use the GF results from the adaptive SNPs. 

To start, we need to load future climate data, which I have already prepared for you as a stack of rasters with the same resolution and extent as those we used for the present day. These data were also downloaded from WorldClim, from an arbitrarily choisen model/scenario of future climate. 

	clim.future <- stack("future.tif")
	names(clim.future) <- c("Pdry", "Pseas", "Tmin", "Tseas")

Next, we need to extract the climate data for all the cells in the climate data layer and remove missing data, as we did prior to making predictions above (here as one command).

	clim.land.future <- na.omit(extract(clim.future, 1:ncell(clim.future), df = TRUE))

Again as above, now we can predict allelic turnover across the whole landscape in the *future* (`clim.land.future`) using our fitted model (`gf.adaptive`) that was based on *present-day* climate.
	
	pred.adaptive.future <- predict(gf.adaptive, clim.land.future[,-1])

To calculate the genetic offset between future and present genotype-environment associations, we can use a simple distance metric, such as Euclidean distance:
	
	genetic.offset.adaptive <- sqrt((pred.adaptive.future[,1]-pred.adaptive[,1])^2 + (pred.adaptive.future[,2]-pred.adaptive[,2])^2 + (pred.adaptive.future[,3]-pred.adaptive[,3])^2 + (pred.adaptive.future[,4]-pred.adaptive[,4])^2)

Finally, we are ready to prepare the raster of differences and map.

	#Define raster properties
	rast.offset <- clim.future$Pseas 
	
	#Assign genetic offset values (difference between future and present predictions) to raster                
	rast.offset[clim.land.future$ID] <- genetic.offset.adaptive
	
	#Plot with color scale of choice
	pdf("GeneticOffset.pdf")
	plot(rast.offset, col=rev( rainbow( 99, start=0, end=0.2 )))
	points(clim.points$Longitude, clim.points$Latitude)
	dev.off()
	
Which areas might be most vulnerable to climate change? What are some limitations and caveats?

### Other approaches

If you finish early or are interested in comparing the results of other approaches, you can try generalized dissimilarity modeling (GDM) or redundancy analysis (RDA). See my [2017 tutorial for details on GDM](https://github.com/pgugger/LandscapeGenomics/blob/master/2017/Exercise4.md#generalized-dissimilarity-modeling-gdm) and my [2015 tutorial for details on RDA](https://github.com/pgugger/LandscapeGenomics/blob/master/2015/Exercise4.md).
