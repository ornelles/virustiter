## Synopsis
This is a suite of tools in R to determine viral titers from fluorescent micrograph sets where one image is a DNA image and another a fluorescent image representing the viral signal. The code requires the `EBImage`, `lattice`, `latticeExtra` and `genefilter` packages.

## Overview
This tools in this package are largely structured to process fluorescent virus titers performed in multi-well plates. Typically, pairs of images are collected at different multiplicities of infection or moi. The moi can be expressed as virions (VP) per cell *or* infectious units (IU) per cell *or* a volume (ml, ul, nl) per cell. Although the default order expects the nuclear (DAPI) image file to come before the corresponding viral antigen image file, different orders can be accommodated.

The sets of images associated with each moi can occur either as individual image files in a single directory where each directory is named for the well such as A1, A2, and so on and the files within a directory are identified sequentially as file001.tif, file002.tif, etc. Alternatively, the pairs of images can be part of a multi-layered TIFF file for each moi where each set of images includes the DNA and viral antigen images.

Additional information about the experiment must be provided in a "phenotype" data frame that describes the conditions of the experiment and includes the moi and unit of measure (as VP,  ml, ul, nl, etc). These data are merged with the image data for further analysis.

Individual cells are identified by the nuclear stain which is used to generate a nuclear mask. This mask is applied to the viral antigen image file and the mean fluorescence intensity is measured for each cell defined by the nuclear mask. An option is provided to expand the size of the nuclear mask to include more of the associated cytoplasm. See the help function for `parseImages()` for more details and additional options to optimize detection. 

## Installation
This is provided as an early release of a package that can be installed from github. A few steps are necessary to install it and related packages before use.

First, the supporting packages `EBImage` and `genefilter` need to be installed from the Bioconductor using the latest version of `biocLite.R`. Be sure to have the latest version of R installed before using `biocLite`.
```
source("https://bioconductor.org/biocLite.R")
biocLite("EBImage")
biocLite("genefilter")
```
Second, the `devtools` and `latticeExtra` packages need to be installed.
```
install.packages("devtools")
install.packages("latticeExtra")
```
Finally, the `virustiter` package can be installed from GitHub.
```
library(devtools)
install_github("ornelles/virustiter")
```

## Working notes
Phenotype date should be a data frame with the following variables:
```
  moi   A numeric value, which can be named 'x', indicating the multiplicity
  unit  A character string identifying the units per cell as "VP "IU "ul or "ml"
```
and must include *either* `well` *or* `file` *but not both*:
```
  well  A character string indicated the well such as "A1" or "a01" or "a0001"
  file	The filename as a character string of the multi-layered TIFF file
```
An example with images in individual files in folders is shown here. `parseImages(path)` will examine the folder specified by `path` in order to determine if it contains multi-layered tiff files or additional folders with individual image files and process the files accordingly. 
```
  path <- system.file("extdata", "by_folder", package = "virustiter")
  fpd <- system.file("extdata", "by_folder/phenoData.csv", package = "virustiter")
  
  df <- parseImages(path)
  pd <- read.csv(fpd)
  df <- mergePdata(pd, df)
  cut <- getCut(df)  # this is less than optimal, see the analysis below
  df <- score(df, cut)
  res <- tally(df)
  fm <- getFit(res)
  plotFit(fm)
```
An example with stacked images in a single folder is shown here. Repeat the above sample code using the new files.
```
  path <- system.file("extdata", "by_stack", package = "virustiter")
  fpd <- system.file("extdata", "by_stack/phenoData.csv", package = "virustiter")
```
Typical workflow:
```
   df <- parseImages()   # read paired images with EBImage
   ...or...
   df <- readIJResults() # read data from Fluorescent Cell Count (ImageJ)

   pd <- data.frame(well = well.info(levels(df$well)), moi = moi, unit = unit)
   ...or...
   pd <- data.frame(file = levels(df$file), moi = moi, unit = unit)

   df  <- mergePdata(pd, df) # merge with phenotype data in 'pd'
   cut <- getCut(df)     # determine cutoff by control (or well, row, or column)
   df  <- score(df, cut) # assign positive values from cutoff
   res <- tally(df)      # tally positives and negatives and return data.frame
   fm  <- getFit(res)    # get model fit(s) from either res or scored data frame
   cf  <- getTiter(fm)   # get value in units required for MOI of 1 and 95% CI
```
Supporting functions include these as well as others:
```
   checkImages(path)   # check (and optionally display) paired images
   plotDens(df)        # show cutoff values with densityplot 
   plotHist(df)        # show cutoff values with histogram
   plotPlate(df)       # plot entire plate showing positives
   plotWell(df, well)  # plot a schematic of each image in a given well(s) or file(s)
   plotFit(fm)         # plot fit(s) with calculated values using base graphics
   plotOneFit(fm)      # plot fit with more options to adjust colors
   addOneFit(fm)       # add another best-fit line and points to an existing plot
   getAIC(df, cut, by) # evaluate fitted model(s) from df at cut values
   nucMask(dapi)       # extract nuclear mask from dapi image(s) or file(s)
   trimMask(mask)      # remove objects based on size from image mask
   cellMask(mask)      # expand a nuclear mask into a cell image mask
   bnormalize(img)     # normalize images to a common background value
   p2p()               # interactively measure point-to-point distances
   pnpoly(p, v)        # test if points in p are within polygon (v)
```
Often the cutoff value needs to be optimized with parameters provided to `getCut()` as well as those initially used such as `width` in `parseImages()`. Use the plotting tools `plotDens()` and `plotHist()` to evaluate the choice of cutoff values.

The sample data provided here yields a less than ideal cutoff using default settings. The control values (moi of 0) are so tight that the default value of `mult = 5` for the 'mad' multiplier is too generous.

The following code demonstrates one method of exploring values near the optimal cutoff value with `getAIC()`. The AIC values point to two possible cutoffs but the results from `plotHist()` show that the value with `mult` = 3 is more sensible.
```
  mm <- seq(2, 6, 0.25)
  cuts <- getCut(df, mult = mm)
  aic <- getAIC(df, cuts)
  plot(mm, aic)	# by AIC, the best mult value is 3 or 5
  plotHist(df, cuts[mm == 3], main = sprintf("Cutoff = %0.4f", cuts[mm==3]))
  plotHist(df, cuts[mm == 5], main = sprintf("Cutoff = %0.4f", cuts[mm==5]))
```  
## License
GPL-3
