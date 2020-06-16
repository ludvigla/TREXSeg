#' Filter function
#'
#' Takes an image of class 'Image' as input and applies a
#' 2D FFT convolution filter.
#'
#' @param im image of class 'Image'
#' @param brush.size diameter size of brush [default: 9 pixels]
#' @param verbose Print messages
#'
#' @importFrom EBImage makeBrush filter2
#'
#' @export
#'
filter_cells <- function (
  im,
  brush.size = 9,
  verbose = FALSE
) {
  if (verbose) cat("Applying 2D convolution filter to image ... \n")
  f = makeBrush(brush.size, shape = 'disc', step = FALSE)
  f = f/sum(f)
  imfiltered <- filter2(im, filter = f)
  return(imfiltered)
}


#' Correction function
#'
#' Takes an image of class 'Image' and a filtered version of this
#' image, also of class 'Image' as input. The filtered version is subtracted
#' from the original image and the intensity values are then normalized to
#' produce a corrected image.
#'
#' @param im Image of class 'Image'
#' @param imfiltered Image of class 'Image'
#' @param verbose Print messages
#'
#' @importFrom EBImage normalize
#'
#' @export
#'
correct_cells <- function (
  im,
  imfiltered,
  verbose = FALSE
) {
  if (class(im) != "Image") stop(paste0("Invalid input format of im: ", class(im)))
  if (class(imfiltered) != "Image") stop(paste0("Invalid input format of imfiltered: ", class(imfiltered)))
  if (verbose) cat("Correcting image ... \n")
  imcorrected <- normalize(im - imfiltered)
  return(imcorrected)
}

#' Threshold function
#'
#' Takes a "corrected" image of class 'Image' as input and segments the image
#' by thresholding the intensity values based on 'nsd' number of
#' standard deviations from the mean intensity.
#'
#' @param imcorrected Image of class 'Image'
#' @param nsd Number of standard deviations [default: 2]
#' @param verbose Print messages
#'
#' @export
#'
threshold_cells <- function (
  imcorrected,
  nsd = 2,
  verbose = FALSE
) {
  if (class(imcorrected) != "Image") stop(paste0("Invalid input format: ", class(imcorrected)))
  if (verbose) cat("Thresholding cells ... \n")
  thr <- mean(imcorrected) + nsd*sd(imcorrected)
  imthreshold <- imcorrected > thr
  return(imthreshold)
}

#' Cleaning function
#'
#' Takes an image of class 'Image' as with segmented cells as input,
#' labels the cells (sets of connected pixels) and applies a threshold
#' with both a lower and an upper bound to remove cells that are
#' considered too small or too large (based on area).
#'
#' @param imthreshold An image of class 'Image'
#' @param thr thresholds (lower and upper bound)
#' @param verbose Print messages
#'
#' @importFrom EBImage bwlabel rmObjects
#'
#' @export
#'
clean_cells <- function (
  imthreshold,
  thr = c(5, 40),
  verbose = FALSE
) {
  if (class(imthreshold) != "Image") stop(paste0("Invalid input format: ", class(imthreshold)))
  if (verbose) cat("Cleaning up unwanted speckles ... \n")
  imthreshold <- bwlabel(imthreshold)
  areas <- table(imthreshold)[-1]
  inds <- which(areas < thr[1] | areas > thr[2])
  imclean <- rmObjects(x = imthreshold, index = inds)
  return(imclean)
}

#' Watershed function
#'
#' Take an image of class 'Image' as input and applies a watershed
#' algorithm of split merged cells.
#'
#' @param imclean Image of class 'Image'
#' @param tol tolerance used for watershedding [default: 0.1]
#' @param verbose Print messages
#'
#' @importFrom EBImage watershed distmap
#'
#' @export
#'
watershed_cells <- function (
  imclean,
  tol = 0.1,
  verbose = FALSE
) {
  if (class(imclean) != "Image") stop(paste0("Invalid input format: ", class(imclean)))
  if (verbose) cat("Applying watershed ... \n")
  imwatershed <- watershed(x = distmap(imclean), tolerance = tol)
  return(imwatershed)
}

#' Segment cells
#'
#' This function can be used on immunofluorescence images to
#' segment out stained cells. The segmentation workflow is computed in 5 steps:
#' \itemize{
#'   \item{filter}
#'   \item{correction}
#'   \item{threshold}
#'   \item{cleaning}
#'   \item{watershed}
#' }
#'
#' @param impath path to immunofluorescence image
#' @param crop.window vector of length four specifying a window to crop out
#' from the provided image. Can be useful if there are artefacts along the edges
#' of the image.
#' @param return.all Return images from each step as a list
#' @param verbose Print messages
#'
#' @inheritParams filter_cells
#' @inheritParams correct_cells
#' @inheritParams threshold_cells
#' @inheritParams clean_cells
#' @inheritParams watershed_cells
#'
#' @importFrom EBImage readImage normalize
#'
#' @export
#'
SegmentCells <- function (
  impath,
  crop.window = NULL,
  brush.size = 9,
  nsd = 2,
  feature.threshold = c(5, 40),
  tolerance = 0.1,
  return.all = FALSE,
  verbose = FALSE
) {
  if (!file.exists(impath)) stop(paste0("File ", impath, " does not exist \n"))
  cells <- readImage(impath)
  cells <- normalize(cells)

  if (!is.null(crop.window)) {
    if (!length(crop.window) == 4 & class(crop.window) %in% c("numeric", "integer")) stop("Invalid crop window \n")
    cells <- cells[crop.window[1]:crop.window[2], crop.window[3]:crop.window[4]]
  }

  cells_filtered <- filter_cells(cells, brush.size, verbose) # 1. filter
  cells_corrected <- correct_cells(cells, cells_filtered, verbose) # 2. correct
  cells_th <- threshold_cells(cells_corrected, nsd, verbose) # 3. threshold
  cells_clean <- clean_cells(cells_th, feature.threshold, verbose) # 4. clean
  cells_split <- watershed_cells(cells_clean, tolerance, verbose) # 5. watershed

  if (return.all) {
    return(list(cells, cells_filtered, cells_corrected, cells_clean, cells_split))
  } else {
    return(cells_split)
  }
}


#' Overlap estimate
#'
#' Calculates the overlap of two sets x and y
#' and their intersect i.
#'
#' @param x size of x
#' @param y size of y
#' @param i size of intersect of x and y
#'
overlap_fkn <- function (
  x,
  y,
  i
) {
  return(i/pmin(x, y))
}

#' Find overlapping signals
#'
#' Takes two segmented images "ima" and "imb" of class 'Image' as input and
#' tries to find cells with overlapping signals by computing shapes with
#' an overlapping area larger than 'overlap.min'.
#'
#' @param ima Image of class 'Image' with segmented data
#' @param imb Image of class 'Image' with segmented data
#' @param overlap.min Minimum overlap between two cells to consider them the same cell
#' @param return.merged Returns overlapping cells but also cells that are strictly non-overlapping
#' @param return.indices Also returns the indices of the removed cells
#' @param verbose Print messages
#'
#' @importFrom EBImage bwlabel
#'
#' @export
#'
OverlapImages <- function (
  ima,
  imb,
  overlap.min = 0.5,
  return.merged = FALSE,
  return.indices = FALSE,
  verbose = FALSE
) {

  if (length(dim(ima)) == 3 | length(dim(imb)) == 3) {
    stop("Invalid dims for ima or imb")
  }
  if (length(table(ima)) <= 2) stop("ima has not been labeled")
  if (length(table(imb)) <= 2) stop("imb has not been labeled")
  # Summarize overlap across images
  if (verbose) cat(paste0("Calculating intersect between images \n"))
  d <- data.frame(a = as.numeric(ima), b = as.numeric(imb))
  d <- table(d) %>% as.matrix()
  d <- d[-1, -1]

  # Extract area for image a and image b
  area_a <- table(ima)[-1]
  aa <- setNames(as.numeric(area_a), names(area_a))
  area_b <- table(imb)[-1]
  ab <- setNames(as.numeric(area_b), names(area_b))
  if (verbose) cat(paste0("Finished calculating intersect between ", length(area_a), " shapes in image a and ", length(area_b), " shapes in image b \n"))

  # Collect indices for overlaping shapes
  intersect_ab <- which(d > 0, arr.ind = T)
  indices <- setNames(data.frame(apply(do.call(rbind, (lapply(1:nrow(intersect_ab), function(i) {
    inds <- intersect_ab[i, ]
    c(rownames(d)[inds[1]], colnames(d)[inds[2]], d[inds[1], inds[2]])
  }))), 2, function(x) as.character(x)), stringsAsFactors = F), nm = c("ind.a", "ind.b", "intersect"))

  # Collect shape areas, indices and intersect
  df <- data.frame(area.a = aa[indices$ind.a], area.b = ab[indices$ind.b],
                   inda = as.integer(indices$ind.a), indb = as.integer(indices$ind.b),
                   intersect = as.numeric(indices$intersect), stringsAsFactors = F)

  # Calculate overlap
  if (verbose) cat(paste0("Calculating overlap between images using shape intersect \n"))
  df$overlap <- overlap_fkn(x = df$area.a, y = df$area.b, i = df$intersect)
  df$keep <- df$overlap > overlap.min
  df <- subset(df, keep)

  # Clean up image a and image b
  ima_clean <- rmObjects(ima, index = setdiff(as.integer(names(table(ima)[-1])), df$inda))
  imb_clean <- rmObjects(imb, index = setdiff(as.integer(names(table(imb)[-1])), df$indb))

  # Should the clean images be merged?
  if (return.merged) {
    imres <- ima_clean | imb_clean
  } else {
    imres <- ima_clean & imb_clean
  }

  # return extra data
  if (return.indices) {
    return(list(bwlabel(imres), df))
  } else {
    return(bwlabel(imres))
  }
}
