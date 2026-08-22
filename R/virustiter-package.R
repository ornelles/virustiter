#' virustiter: Process fluorescent images for virus titer
#'
#' A suite of code that helps determining virus titers from fluorescent
#' micrographs, typically from paired fluorescent images of the cell nucleus
#' and a viral antigen. 
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import EBImage
#' @import lattice
#' @importFrom MASS dose.p fitdistr
#' @importFrom multimode locmodes modetest
#' @importFrom stats AIC glm qnorm
#' @importFrom latticeExtra xscale.components.log10ticks
#' @importFrom latticeExtra xscale.components.logpower
#' @importFrom graphics Axis abline box lines locator par points title
#' @importFrom methods is
#' @importFrom stats IQR aggregate as.formula binomial density mad median
#' @importFrom stats optim predict quantile runmed setNames
#' @importFrom utils flush.console head setTxtProgressBar tail txtProgressBar 
#' @importFrom utils unzip 
## usethis namespace: end
NULL
