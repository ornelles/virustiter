## Synopsis
This is a suite of tools in R to determine viral titers from fluorescent micrograph pairs where the first of each pair is an (overexposed) DNA image and the second a fluorescent image representing the viral signal. The code requires the `EBImage` and `lattice` packages.

## Overview
This was developed to process multi-well plates. Most of the tools are structured for this purpose. Pairs of images are collected at different multiplicities of infection (moi): viral particles or infectious units or volume per cell. DNA (DAPI) images must always come before the  fluorescent viral image file. 

Pairs of images associated with each moi can be individual files in a single folder (file001.tif, file002.tif, etc.) or multi-layered tif files for each moi with the first image in each file being the DNA image. 

## Installation

Approximately: use `bioclite` to install `EBImage`. Then use  `devtools::install_github("ornelles/virustiter")`.

## Working notes
Phenotype date should be a data frame with the following variables:
```
  moi   numeric value indicating the multiplicity (units per cell)
  unit  character string indicating the unit per cell as "VP "IU "ul or "ml"
```
and then either `well` or `file`
```
  well  character string indicated the well such as "A1" or "a01"
  file	character string identifying the file holding the layered TIF
```

Example with images as individual files in folders:
```
  fimg <- system.file("extdata", "by_folder/b1/file001.tif", package = "virustiter")
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
Example with stacked images in a single folder (repeat sample code above)
```
  fimg <- system.file("extdata", "by_stack/file001.tif", package = "virustiter")
  fpd <- system.file("extdata", "by_stack/phenoData.csv", package = "virustiter")
```
Working functions:
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
  plotCut(df)     # calculate and show cutoff values with densityplot 
  plotPlate(df)   # plot plate showing positives
  plotWell(df, well) # plot each file in a well showing positives and sizes
  plotHist(df)    # histogram on well-by-well basis with optional cutoff values
  plotFit(fm)     # plot fit(s) with calculated values using base graphics
  plotOneFit(fm)  # plot fit with options to adjust colors
  addOneFit(fm)   # add best-fit line to existing base graph
  getAIC(df, cut, by)  #evaluate fitted model(s) from df at cut value(s)
  displayPairs(f, dna = TRUE) # display image pairs in directory containing 'f'
```
Wrapper to automatically process results data frame or ImageJ 'Results.txt' file
``` 
 fitAndPlot(res, by)")
```
To optimize the fit, the cutoff value needs to be tuned with parameters handed to `getCut()` as well as those initially used such as `width` in  `parseImages()`. Use the graphing tools (`plotCut()` and `plotHist()`) to evaluate the cutoff.

For example, the sample data yields a non-optimal cutoff using default settings. The control values (moi of 0) are so tight that the default value of 5 for the 'mad' multiplier (`mult`) is too generous.

The following shows a better initial selection followed by further exploration using values near the optimal cutoff value with `getAIC()`. The AIC values point to two possible cutoffs but the plot shows that the value with `mult` = 3 is more sensible.
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
