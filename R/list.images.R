#' List Image Files
#'
#' List image files in the given path.
#'
#' @param path A character vector of path names.
#' @param type The type of image as defined by the file extension (tiff, jpeg, jpg or png).
#' @param pattern Argument passed to list.files (grep pattern), default is NULL.
#' @param full.names Argument passed to list.files, default is TRUE.
#' @param recursive Argument passed to list.files, default is TRUE.
#' @param ignore.case Argument passed to list.files, default is TRUE.
#'
#' @return
#'
#' Character vector of selected image filenames.
#'
#' @examples
#' path.by.folder <- system.file("extdata", "by_folder", package = "virustiter")
#' list.images(path.by.folder, pattern = "002")
#'
#' @export
#' 
list.images <- function(path = ".", type = c("tiff", "jpeg", "jpg", "png"),
	pattern = NULL, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
{
# assign file extension pattern from 'type'
	ext <- switch(match.arg(type),
		tiff = "\\.tif{1,2}$",
		jpeg = "\\.jpe{0,1}g$",
		jpg = "\\.jpe{0,1}g$",
		png = "\\.png$")
# collect file list matching pattern and then by type (extension)
	ff <- list.files(path = path, pattern = pattern, full.names = full.names,
		recursive = recursive, ignore.case = ignore.case)
	ff <- ff[grep(ext, ff, ignore.case = TRUE)]
	return(ff)
}
