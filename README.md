## Synopsis
This is a suite of tools in R to analyze fluorescent micrographs from a fluorescent focus titration assay to determine viral titers. The code makes use of the `lattice` and `EBImage` package.

## Overview
  `f <- system.file("extdata", "sample/b1/file001.tif", package = "virustiter"`
## Installation

## Typical workflow
```
Phenotype date should be a data frame with at least each of the following:
  unit  character string indicating the unit per cell as "VP "IU "ul or "ml"
  well  character string indicated the well such as "A1" or "a01"
  moi   numeric value indicating the multiplicity (units per cell)
See for example:
  fpd <- system.file("extdata", "phenoData.csv", package = "virustiter")
  pd <- read.csv(fpd)
  fi <- system.file("extdata", "sample/b1/file001.tif", package = "virustiter")
  df <- parseImages(fi)

Data for parseImages() comes as paired TIF files with high contrast DAPI images
paired with fluorescent images for viral marker.

Working functions:
  df <- parseImages()   # read with EBImage
     or
  df <- readData()      # read data from Fluorescent Cell Count v6 (ImageJ)

  df  <- mergePdata(pd, df)  # optional merge with phenoData in 'pd'
  cut <- getCut(df)     # determine cutoff by control (or well, row, or column)
  df  <- score(df, cut) # assign positive values from cutoff
  res <- tally(df)      # tally positives and negatives and return data.frame
  fm  <- getFit(obj)    # get model fit(s) from scored data frame (df) or from 'res'
  cf  <- getTiter(fm)   # get value in units required for MOI of 1 and 95% CI

Supporting functions:
  plotCut(df)     # calculate and show cutoff values with densityplot 
  plotPlate(df)   # plot plate showing positives
  plotWell(df, well) # plot each file in a well showing positives and sizes
  plotHist(df)    # histogram on well-by-well basis with optional cutoff values
  plotFit(fm)     # plot fit(s) with calculated values using base graphics
  plotOneFit(fm)  # plot fit with options to adjust colors
  addOneFit(fm)   # add best-fit line to existing base graph
  getAIC(df, cut, by)  #evaluate fitted model(s) from df at cut values
  displayPairs(f, dna = TRUE) # display image pairs in directory containing 'f'

Wrapper to automatically process results data frmae or ImageJ 'Results.txt' file
  fitAndPlot(res, by)")
```

## License
GPL-3
