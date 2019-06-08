## Synopsis
This is a suite of imaging tools developed to determine viral titers from fluorescent micrographs. In typical use, one image is a DNA image and another a fluorescent image representing the viral signal. The code requires the `EBImage`, `lattice`, and `latticeExtra` packages.

## Overview
The tools in this package have primarily been developed to perform fluorescent viral titers in multi-well culture dishes. Typically, pairs of images are collected at different multiplicities of infection or **moi**. The moi can be expressed as virions (VP) per cell *or* infectious units (IU) per cell *or* a volume (ml, ul, nl) per cell. Although the default order has the nuclear (DAPI) image file before the corresponding viral antigen image file, different orders can be accommodated.

The sets of images associated with each moi can occur either as individual image files in a single directory where each directory is named for the well such as A1, A2, and so on and the files within a directory are identified sequentially as file001.tif, file002.tif, etc. Alternatively, the pairs of images can be part of a multi-layered TIFF file for each moi where each set of images includes the DNA and viral antigen images.

Additional information about the experiment must be provided in a "phenotype" data frame that describes the conditions of the experiment and includes the moi and unit of measure (as VP,  ml, ul, nl, etc). These data are merged with the image data for further analysis.

Individual cells are identified by a DNA stain which is used to generate a nuclear mask. This mask is applied to the viral antigen image file and the mean fluorescence intensity is measured for each cell defined by the nuclear mask. An option is provided to expand or contract the size of the nuclear mask in order to include more or less of the associated cytoplasm. See the help function for `parseImages()` and `trimMask()` for more details and additional options to optimize detection. 

## Installation
This is revision 3 of the second "release" of a package that can be installed from github. Functions previously embedded in `parseImages()` are now split between `getImages()` and `parseImages()`. A few steps are necessary to install it and related packages before use.

First, the supporting package `EBImage` must be installed from the Bioconductor using the latest version of `biocLite.R`. Be sure to have the latest version of R installed before using `biocLite`.
```
source("https://bioconductor.org/biocLite.R")
biocLite("EBImage")
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
  unit  A character string identifying the units per cell such as "VP "IU "ul or "ml"
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
  
  img <- getImages(path)
  df <- parseImages(img)
  pd <- read.csv(fpd)
  df <- mergePdata(pd, df)
  bg <- getBgnd(df)  # this is less than optimal, see the analysis below
  df <- score(df, bg)
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
   img <- getImages()     # read paired images with EBImage
   df <- parseImages(img) # extract nuclear information and target mfi

   pd <- data.frame(well.info(levels(df$well)), moi = moi, unit = unit)
   ...or...
   pd <- data.frame(file = levels(df$file), moi = moi, unit = unit)

   df  <- mergePdata(pd, df) # merge with phenotype data in 'pd'
   bg <- getBgnd(df)     # determine background level by control
   df  <- score(df, bg)  # assign positive values from bg
   res <- tally(df)      # tally positives and negatives and return data.frame
   fm  <- getFit(res)    # get model fit(s) from either res or scored data frame
   cf  <- getTiter(fm)   # get value in units required for MOI of 1 and 95% CI
```
Supporting functions include these as well as others:
```
   checkImages(path)   # check (and optionally display) paired images
   plotDens(df)        # show default cutoff values with densityplot 
   plotHist(df)        # show default cutoff values with histogram
   plotPlate(df)       # plot entire plate showing positives
   plotWell(well, df)  # plot a schematic of each image in a given well(s) or file(s)
   plotFit(fm)         # plot fit(s) with calculated values using base graphics
   plotOneFit(fm)      # plot fit with more options to adjust colors
   addOneFit(fm)       # add another best-fit line and points to an existing plot
   getAIC(df, bg, by)  # evaluate fitted model(s) from df at various background values
   nucMask(dapi)       # extract nuclear mask from dapi image(s) or file(s)
   trimMask(mask)      # remove objects based on size from image mask
   cellMask(mask)      # expand a nuclear mask into a cell image mask
   bnormalize(img)     # normalize images to a common background value
   getZero(img)        # determine the background (zero value) pixel
   setZero(img, zero)  # normalize images to a common zero value
   p2p()               # interactively measure point-to-point distances
   pnpoly(p, v)        # test if points in p are within polygon v
```
Often the background value needs to be optimized with parameters provided to `getBgnd()` as well as those provided to `parseImages()`. Use the plotting tools `plotDens()` and `plotHist()` to evaluate different choices of background values.

The sample data provided here yields a less than ideal cutoff value for the background when using default settings. The control values (moi of 0) are so tight that the default value of `mult = 2.5` for the 'mad' multiplier is too generous.

The following code demonstrates one method of exploring values near the optimal cutoff value with `getAIC()`. The AIC values point to two possible background values. The results from `plotHist()` show that the values with `mult` = 1.55 or 1.95 may be better than the default value of 2.5.
```
  mm <- seq(1, 2.5, 0.05)
  bg <- sapply(mm, function(m) getBgnd(df, mult = m))
  aic <- getAIC(df, bg)
  plot(mm, aic, type = "b")	# by AIC, the best mult value is 2.65
  plotHist(df, bg[mm == 2.5], main = "Default 'mult' value = 2.5")
  opt.mm <- mm[which.min(aic)]
  plotHist(df, bg[opt.mm], main = sprintf("Optimal 'mult' value = %0.2f", opt.mm))
```  
## License
GPL-3
