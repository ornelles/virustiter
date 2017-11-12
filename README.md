## Synopsis
This is a suite of tools in R to determine viral titers from fluorescent micrograph pairs where the first of each pair is an (ideally overexposed) DNA image and the second a fluorescent image representing the viral signal. The code requires the `EBImage`, `lattice`, `latticeExtra` and  `genefilter` packages.

## Overview
This was developed to process multi-well plates. Most of the tools are structured for this purpose. Pairs of images are collected at different multiplicities of infection or moi expressed as virions (VP) *or* infectious units (IU) *or* volume (ml, ul, nl, etc.) per cell. The nuclear (DAPI) image file must always come before the corresponding viral antigen image file. 

Pairs of images associated with each moi can be individual files where each folder is named for the well such as a1/file001.tif, a2/file002.tif. Alternatively, the pairs of images can be a single multi-layered tif file for each moi where the first image in each pair file is the DNA image. 

## Installation

**Not yet ready for installation.** The `EBImage` and `genefilter` packages will have to be installed with Bioclite. The CRAN package `latticeExtra` will need to be installed as well. After that, the contents can be cloned and "installed" locally from the local directory with   `devtools::load_all()`.

## Working notes
Phenotype date should be a data frame with the following variables:
```
  moi   numeric value indicating the multiplicity (units per cell)
  unit  character string indicating the unit per cell as "VP "IU "ul or "ml"
```
and also must include either `well` or `file`:
```
  well  character string indicated the well such as "A1" or "a01"
  file	character string identifying the file holding the layered TIF
```

An example with images in individual files in folders is shown here. Note that any file in the subfolder can be used to point `parseImages()` to the collection of files. This function will try to determine if this file is a multilayered tiff file or a single image in a collection of image files and process the files accordingly. 
```
  fimg <- system.file("extdata", "by_folder/b2/file002.tif", package = "virustiter")
  fpd <- system.file("extdata", "by_folder/phenoData.csv", package = "virustiter")

  df <- parseImages(fimg)
  pd <- read.csv(fpd)
  df <- mergePdata(pd, df)
  cut <- getCut(df)  # this is not optimal, see analysis below
  df <- score(df, cut)
  res <- tally(df)
  fm <- getFit(res)
  plotFit(fm)
```
An example with stacked images in a single folder is shown here. Repeat the above sample code using the new files.
```
  fimg <- system.file("extdata", "by_stack/file001.tif", package = "virustiter")
  fpd <- system.file("extdata", "by_stack/phenoData.csv", package = "virustiter")
```
Typical workflow:
```
  df <- parseImages()   # read with EBImage
     or
  df <- readData()      # read data from Fluorescent Cell Count v6 (ImageJ)

  df  <- mergePdata(pd, df)  # optional merge with phenoData in 'pd'
  cut <- getCut(df)     # determine cutoff by control (or well, row, or column)
  df  <- score(df, cut) # assign positive values from cutoff
  res <- tally(df)      # tally positives and negatives and return data.frame
  fm  <- getFit(obj)    # get model fit(s) from scored data frame (df) or from 'res'
  cf  <- getTiter(fm)   # get value in units required for MOI of 1 and 95% CI
```
Supporting functions:
```
  plotCut(df)    # calculate and show cutoff values with densityplot 
  plotHist(df)   # histogram of each well with optional cutoff values
  plotPlate(df)  # lattice plot plate showing cell position and positives
  plotWell(df, well) # plot each file in a well showing positives and cell size
  plotFit(fm)    # plot fit(s) with calculated values using base graphics
  plotOneFit(fm) # plot fit with options to adjust colors
  addOneFit(fm)  # add best-fit line to existing base graph
  getAIC(df, cut, by)  # evaluate fitted model(s) from df at cut values
  displayPairs(f) # display one of paired images with nuclear mask
  nucMask(dapi)   # extract nuclear mask from dapi image(s) or file(s)
  p2p()           # interactively measure point-to-point distances
  pnpoly(p, v)    # test if points in p are within polygon (v)
```
Wrapper to automatically process results data frame or ImageJ 'Results.txt' file
``` 
 fitAndPlot(res, by)")
```
To optimize the fit, the cutoff value needs to be tuned with parameters handed to `getCut()` as well as those initially used such as `width` in  `parseImages()`. Use the graphing tools `plotCut()` and `plotHist()` to evaluate the choice of cutoff.

The sample data provided here yields a less than ideal cutoff using default settings. The control values (moi of 0) are so tight that the default value of `mult = 5` for the 'mad' multiplier is too generous.

The following shows a better initial selection followed by further exploration using values near the optimal cutoff value with `getAIC()`. The AIC values point to two possible cutoffs but the results from `plotHist` show that the value with `mult` = 3 is more sensible.
```
  mm <- seq(2, 6, 0.25)
  cuts <- getCut(df, mult = mm)
  aic <- getAIC(df, cuts)
  plot(mm, aic)	# by AIC, the best mult value is 3 or 5
  plotHist(df, cuts[mm==3], main = sprintf("Cutoff = %0.4f", cuts[mm==3]))
  plotHist(df, cuts[mm==5], main = sprintf("Cutoff = %0.4f", cuts[mm==5]))
```  
## License
GPL-3
