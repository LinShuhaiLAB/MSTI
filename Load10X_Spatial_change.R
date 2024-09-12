
Load10X_Spatial_change <- function(data.dir, filename = "filtered_feature_bc_matrix.h5", 
                                   image.name = "tissue_lowres_image.png",
                                   assay = "Spatial", slice = "slice1", filter.matrix = TRUE, 
                                   to.upper = FALSE, image = NULL, ...) 
{
  if (length(x = data.dir) > 1) {
    warning("'Load10X_Spatial' accepts only one 'data.dir'", 
            immediate. = TRUE)
    data.dir <- data.dir[1]
  }
  if(grepl("h5", filename )){
    data <- Read10X_h5(filename = file.path(data.dir, filename), 
                       ...)
  } else {
    #data = Read10X(filtered_dir)
    data = Read10X(paste(data.dir, filename,sep = "/") )
  }
  if (to.upper) {
    rownames(x = data) <- toupper(x = rownames(x = data))
  }
  object <- CreateSeuratObject(counts = data, assay = assay)
  if (is.null(x = image)) {
    image <- Read10X_Image(image.dir = file.path(data.dir, 
                                                 "spatial"), 
                           image.name = image.name,
                           filter.matrix = filter.matrix)
  }
  else {
    if (!inherits(x = image, what = "VisiumV1")) 
      stop("Image must be an object of class 'VisiumV1'.")
  }
  image <- image[Cells(x = object)]
  DefaultAssay(object = image) <- assay
  object[[slice]] <- image
  if(image.name == "tissue_lowres_image.png") {
    object = object
  }else {
    object@images[[1]]@scale.factors$lowres = object@images[[1]]@scale.factors$hires
    #object@images$slice1@scale.factors$lowres = object@images$slice1@scale.factors$hires
    
  }
  return(object)
}
