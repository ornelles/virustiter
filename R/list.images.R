# list image files in enclosing directory or directories
#
#	path		starting file path 
#	enclosing	levels up from path to start recursive search (0, 1, 2, ...)
#	type		image type defined by extension (tiff, jpeg or png)
#	pattern		passed to list.files
#	full.names	passed to list.files
#	recursive	passed to list.files
#	ignore.case	passed to list.files
#
# returns		character vector of matching filenames
#
list.images <- function(path = ".", enclosing = 1,
	type = c("tiff", "jpeg", "jpg", "png"), pattern = NULL, full.names = TRUE,
	recursive = TRUE, ignore.case = TRUE)
{
# set extension pattern based on type of image file
	ext <- switch(match.arg(type),
		tiff = "tif{1,2}$",
		jpeg = "jpe{0,1}g$", jpg = "jpe{0,1}g$",
		png = "png$")
# adjust path to level specified by 'enclosing'
	while (enclosing > 0) {
		path <- dirname(path)
		enclosing <- enclosing - 1
	}
# collect file list matching pattern and then by type (extension)
	ff <- list.files(path = path, pattern = pattern, full.names = full.names,
		recursive = recursive, ignore.case = ignore.case)
	ff <- ff[grep(ext, ff, ignore.case = TRUE)]
	return(ff)
}
