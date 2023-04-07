## Synopsis
This is a suite of imaging code to determine viral titers from fluorescent micrographs. In a typical application, paired fluorescent images of the cell nucleus and a viral antigen are collected analyzed. The code requires the `EBImage`, `lattice`, and `latticeExtra` packages and can make use of the (private) `EBImageExtra` package.

## Overview
The tools in this package have primarily been developed to determine viral titers from multi-well culture dishes. Typically, pairs of images are collected at different multiplicities of infection or **moi**. The moi can be expressed as virions (VP) per cell *or* infectious units (IU) per cell *or* a volume (ml, ul, nl) per cell. The expected order of images is the nucleus (typically DAPI) image file before the corresponding viral antigen image file. However, different orders can be accommodated with optional arguments.

The sets of images associated with each moi can occur in two forms. Most typically, they are collected as individual image files in separate directories where each directory is named for the well such as A1, A2, and so on. In this case, the files within a directory are identified sequentially as file001.tif, file002.tif, etc. Alternatively, the paired images for each moi can be a set of multi-layered TIFF file where the member of each set includes the DNA and viral antigen images.

Additional information about the experiment must be provided in a "phenotype" data frame that **must** contain the moi and unit of measure (as VP,  ml, ul, nl, etc). This data frame can contain additional information describing the experiment. The "phenotype" data are merged with the image data for further analysis.

Individual cells are identified by the DNA stain which is used to generate a nuclear mask. This nuclear mask is applied to the viral antigen image file and the mean fluorescence intensity is measured for each cell defined by the nuclear mask. An option is provided to expand or contract the size of the nuclear mask in order to include more or less of the associated cytoplasm. See the help function for `parseImages()`, `trimMask()` and `cellMask()` for more details and additional options to optimize detection. 

## Significant Changes in Version 0.0.6.0
`getCutoff` has been added to automatically identify and optionally plot optimal cutoff values for trimming nuclear masks.

## Significant Changes in Version 0.0.5.2
`mergePdata` now has the `formatString` option to force that the phenodata and image data use the same format for wells (a conflict such as "A1" vs "A01" can be fixed.)

The `parseImages` now should no longer fail if a fluorescent target image has no cells.

Data can now include a prefix that serves as an identifier to the plate. The unique identifier for each well is now called `label`. This means that a well identified as 2A6 is parsed as `label = 2A6`, `prefix = 1`, `well = "A6"`, `row = "A"` and `column = 6`. Not all of the code (`mergePdata` in particular) accommodates this new variable and new organization. Helpful (hopefully) information is provided where needed. Note that using a variable named `label` in another context may create a conflict with this. Not entirely certain.

A `which.images = NULL` option has be introduced for `checkImages` to handle images that are not organized as paired images. 

## Significant Changes in Version > 0.0.3.2
The `setZero` function was incorrect in versions prior to 0.0.3.2 by failing to replicate the background pixel value for each frame and allowing it to be replicated across frames instead. 

Area calculations in `getVal` and `trimMask` now use `tabulate` for an 8-fold greater efficiency. In version 0.0.3.x, `getTiter` has been changed to actually report the titer as infectious units per ml (if the unit attribute is a known volume unit). The previous functionality is preserved in `getEC63`. The `getVal` function has been simplified with logic to assign the `computeFeatures` function according to the parameter being measured. The arguments for `parseImages` have been shortened (back) to `args.nMask` and `args.cMask`

## Installation
A few steps are necessary to install this code from github and related packages before use. First, the supporting package `EBImage` must be installed from the Bioconductor. As of R.3.6, this is done with `BiocManager` as follows. Be sure to have the latest version of R installed before using `BiocManager`.
```
 if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("EBImage")
```
Second, ensure that the `devtools` and `latticeExtra` packages are installed.
```
  install.packages("devtools")
  install.packages("latticeExtra")
```
Finally, the `virustiter` package can be installed from GitHub.
```
  library(devtools)
  install_github("ornelles/virustiter")
```
Optionally, it may be useful to install the tools in `EBImageExtra` from Github although this is still under development.
```
  require(devtools)
  install_github("ornelles/EBImageExtra")
```
Once the necessary packages are installed, the need to be loaded like a typical R package.
```
  library(virustiter)
  library(EBImageExtra) # If desired/needed
```

## Working notes
Phenotype data must be a data frame with the following variables:
```
  moi   A numeric value, which can be named 'x', indicating the multiplicity
  unit  A character string identifying the units per cell such as "VP "IU "ul or "ml"
```
and must include *either* `well` *or* `file` *but not both*:
```
  well  A character string indicated the well such as "A1" or "a01" or "a0001"
  file  The filename as a character string of the multi-layered TIFF file
```
An example with images in individual files in folders is shown here. `parseImages(path)` will examine the folder specified by `path` in order to determine if it contains multi-layered tiff files or additional folders with individual image files and process the files accordingly. `checkImages(path)` performs the same logic without processing the data. This function is useful to ensure that the images are organized appropriately before analysis.
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

   pd <- data.frame(file = levels(df$file), moi = moi, unit = unit),
     Or for a single plate organized by well:,
   pd <- data.frame(well.info(unique(df$well)), moi = moi, unit = unit),
     Or if multiple plates are used:,
   pd <- data.frame(well.info(unique(df$labels)), moi = moi, unit = unit),

   df  <- mergePdata(pd, df) # merge with phenotype data in 'pd'
   bg <- getBgnd(df)     # determine background level by control
   df  <- score(df, bg)  # assign positive values from bg
   res <- tally(df)      # tally positives and negatives and return data.frame
   fm  <- getFit(res)    # get model fit(s) from either res or scored data frame
   cf  <- getTiter(fm)   # get value in units required for MOI of 1 and 95% CI
```
Supporting functions include these as well as others:
```
   usage()             # show this list of supporting functions
   checkImages(path)   # check (and display) paired images
   list.images(path)   # list image files in the given path
   getZero(x)          # find most common non-zero pixel in each image frame
   setZero(x)          # rescale each image frame to the most common zero value
   nucMask(dapi)       # extract nuclear mask from dapi image(s) or file(s)
   trimMask(mask)      # remove objects based on size from mask
   cellMask(mask)      # expand a nuclear mask into a cell mask
   edgeObjects(mask)   # identify objects near the edge of a mask
   findObjects(expr, df) # find objects in data.frame identified by expr
   getVal(mask, ref)   # extract one 'computeFeatures' value using mask and ref
   p2p()               # interactively measure point-to-point distances
   plotHist(df)        # histogram of each well with optional background values
   plotDens(df)        # calculate and show background values with densityplot
   plotPlate(df)       # plot plate showing positives
   plotWell(df, well)  # plot each cell in a well showing positives and sizes
   plotFit(fm)         # plot fit(s) with calculated values using base graphics
   plotOneFit(fm)      # plot fit with options to adjust colors
   addOneFit(fm)       # add best-fit line to existing base graph
   getAIC(df, bg)      # evaluate fitted model(s) from df at bg values
   getTiter(fm)        # get titer as IU per ml (if volume units were used)
   getEC63(fm)         # get volume required for 1 infectious units +/- 95% CI
   getShift(mask, tgt) # get optimal x-y shift to align tgt with (nuc) mask
   translate(x, v)     # apply (optimal) x-y shift in v to image x
```
Often the background value needs to be optimized with parameters provided to `getBgnd()` as well as those provided to `parseImages()`. Use the plotting tools `plotDens()` and `plotHist()` to evaluate different choices of background values.

The sample data provided here yields a less than ideal cutoff value for the background when using default settings. The control values (moi of 0) are so tight that the default value of `mult = 2.5` for the 'mad' multiplier is too generous.

The following code demonstrates one method of exploring values near the optimal cutoff value with `getAIC()`. The AIC values point to two possible background values. The results from `plotHist()` show that the values with `mult` = 1.55 or 1.95 may be better than the default value of 2.5.
```
  mm <- seq(1, 2.5, 0.05)
  bg <- sapply(mm, function(m) getBgnd(df, mult = m))
  aic <- getAIC(df, bg)
  plot(mm, aic, type = "b") # by AIC, the best mult value is 2.65
  plotHist(df, bg[mm == 2.5], main = "Default 'mult' value = 2.5")
  opt.mm <- mm[which.min(aic)]
  plotHist(df, bg[opt.mm], main = sprintf("Optimal 'mult' value = %0.2f", opt.mm))
```  
## Stuff to do
- Develop a good vignette for this and the `EBImageExtra` package
- Develop a workflow for determining the optimal nuclear mask
- Continue editing the text of this to be a useful stand-alone instruction
- Adjust the code to work more robustly with the new `prefix` variable 

## License
GPL-3
