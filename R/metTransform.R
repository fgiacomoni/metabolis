#### metTransform (mset) ####

#' Transformation of the dataMatrix
#'
#' Transformation of the dataMatrix
#'
#' @param x An S4 object of class \code{MultiDataSet}
#' @param methodC Character: Factor of interest (name of a column from the
#' pData(x))
#' @param info.txtC Character: File name for the printed results (call to
#' 'sink()'); if NA (default), messages will be printed on the screen; if NULL,
#' no verbose will be generated
#' @return \code{MultiDataSet} including the adjusted p-values in the fData data frames
#' @rdname metTransform
#' @export
#' @examples
#' prometMset <- metRead(system.file("extdata/promet", package="metabolis"))
#' prometMset <- metTransform(prometMset, "sqrt")
#' printStr <- lapply(assayData(prometMset), function(aData) ropls::strF(aData[["exprs"]]))
setMethod("metTransform", signature(x = "MultiDataSet"),
          function(x,
                   methodC = c("log2", "log10", "sqrt")[1],
                   info.txtC = NA) {
            
            if (!is.null(info.txtC) &&
                !is.na(info.txtC))
              sink(info.txtC, append = TRUE)
            
            infTxtC <- info.txtC
            if (!is.null(infTxtC))
              infTxtC <- NA

            for (setC in names(x)) {
              
              if (!is.null(info.txtC))
                cat("\nTransforming the '",
                    setC,
                    "' dataset:\n",
                    sep = "")
              
              ese <- metabolis::metTransform(x[[setC]],
                                              methodC,
                                              infTxtC)
              
              x <- MultiDataSet::add_eset(x,
                                          ese,
                                          dataset.type = setC,
                                          GRanges = NA,
                                          overwrite = TRUE,
                                          warnings = FALSE)
              
            }
            
            if (!is.null(info.txtC) &&
                !is.na(info.txtC))
              sink()
            
            return(invisible(x))

          })


#### metTransform (eset) ####

#' Transformation of the dataMatrix
#'
#' Transformation of the dataMatrix
#'
#' @param x An S4 object of class \code{ExpressionSet}
#' @param methodC Character: Factor of interest (name of a column from the
#' pData(x))
#' @param info.txtC Character: File name for the printed results (call to
#' 'sink()'); if NA (default), messages will be printed on the screen; if NULL,
#' no verbose will be generated
#' @return \code{ExpressionSet} including the adjusted p-values in fData (or the vector of adjusted p-values when returnAdjustOnlyL is TRUE)
#' @rdname metTransform
#' @export
#' @examples
#' sacSet <- metRead(system.file("extdata/sacurine", package = "metabolis"))
#' ropls::strF(t(exprs(sacSet)))
#' sacSet <- metTransform(sacSet, "log2")
#' ropls::strF(t(exprs(sacSet)))
setMethod("metTransform", signature(x = "ExpressionSet"),
          function(x,
                   methodC = c("log2", "log10", "sqrt")[1],
                   info.txtC = NA) {
            
            if (!is.null(info.txtC) &&
                !is.na(info.txtC))
              sink(info.txtC, append = TRUE)
            
            x <- .metTransform(x,
                               methodC)

            validObject(x)
            
            if (!is.null(info.txtC) &&
                !is.na(info.txtC))
              sink()

            return(invisible(x))

          })

.metTransform <- function(eset,
                            methC = c("log2", "log10", "sqrt")[1]) {
  
  ## checking		     
  
  if (length(which(Biobase::exprs(eset) < 0)))
    stop("The 'dataMatrix' contains negative values", call. = FALSE)
  
  ## Number of missing values
  nasN <- length(which(is.na(Biobase::exprs(eset))))
  cat("\nMissing values in the 'dataMatrix': ",
      nasN,
      " (",
      round(nasN / cumprod(dim(Biobase::exprs(eset)))[2] * 100),
      "%)\n",
      sep = "")
  
  ## Number of zero values
  zerN <- length(which(Biobase::exprs(eset) == 0))
  cat("\nZero values in the 'dataMatrix': ",
      zerN,
      " (",
      round(zerN / cumprod(dim(Biobase::exprs(eset)))[2] * 100),
      "%)\n",
      sep = "")
  
  ## metTransform
  
  switch(methC,
         log2 = {
           
           cat("\n'log2' transformation\n", sep = "")
           
           Biobase::exprs(eset) <- log2(1 + Biobase::exprs(eset))
           
         },
         log10 = {
           
           cat("\n'log10' transformation\n", sep = "")
           
           Biobase::exprs(eset) <- log10(1 + Biobase::exprs(eset))
           
         },
         sqrt = {
           
           cat("\n'Square root' transformation\n", sep = "")
           
           Biobase::exprs(eset) <- sqrt(Biobase::exprs(eset))
           
         })
  
  eset
  
}