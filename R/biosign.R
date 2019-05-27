#### biosign (mset) ####

#' Feature selection
#'
#' Feature selection
#'
#' @param x An S4 object of class \code{MultiDataSet}
#' @param y a character indicating the name of the column of the pData
#' to be used, when x is an ExpressionSet object
#' @param seedI integer: optional seed to obtain exactly the same signature when rerunning biosigner
#' @param fig.pdfC Character: Name of the file for the graphics from the significant features;
#' if NA, the graphics are displayed on the screen; if NULL, no plot is generated
#' @param plotTierMaxC Character: Maximum level of tiers to display: Either 'S' and
#' 'A', (for boxplot), or also 'B', 'C', 'D', and 'E' (for tiers) by decreasing
#' number of selections
#' @param info.txtC Character: File name for the printed results (call to 'sink()')
#' @param ... additional parameters for the ropls::opls modeling
#' @return list of 'mset', the \code{MultiDataSet} including the tiers in fData 
#' and 'signLs', the list of 'biosign' outputs 
#' @rdname biosign
#' @export
#' @examples
#' prometMset <- readSet(system.file("extdata/promet", package="metabolis"))
#' ## selecting a seed (for reproducibility) and a low bootstrap number (for speed in this demo)
#'\dontrun{
#' set.seed(123)
#' prometMset <- biosign(prometMset, "gene", bootI = 5)
#' ## viewing selected features:
#' lapply(fData(prometMset), function(fda) {
#'   tierDF <- fda[, grep("biosign", colnames(fda))]
#'   tierDF[rowSums(tierDF == "S") > 0, ]
#' })
#' }

setMethod("biosign", signature(x = "MultiDataSet"),
          function(x,
                   y,
                   seedI = NULL,
                   fig.pdfC = NA,
                   plotTierMaxC = "S",
                   info.txtC = NA,
                   ...) {
            
            if (!is.null(info.txtC) && !is.na(info.txtC))
              sink(info.txtC, append = TRUE)
            
            infTxtC <- info.txtC
            if (!is.null(infTxtC))
              infTxtC <- NA
            
            if (!is.null(fig.pdfC) && !is.na(fig.pdfC))
              pdf(fig.pdfC)
            
            msetSignLs <- vector(mode = "list",
                                 length = length(names(x)))
            names(msetSignLs) <- names(x)
            
            for (setC in names(x)) {
              
              if (!is.null(info.txtC))
                cat("\nSelecting features for the '",
                    setC,
                    "' dataset:\n",
                    sep = "")
              
              if (!is.null(seedI))
                set.seed(seedI)
              
              optWrnN <- options()$warn
              options(warn = -1)
              
              msetSignLs[[setC]] <- biosigner::biosign(x[[setC]],
                                                       y = y,
                                                       fig.pdfC = NULL,
                                                       info.txtC = infTxtC,
                                                       ...)
              
              options(warn = optWrnN)
              
              if (!is.null(seedI))
                set.seed(NULL)
              
              if (sum(apply(msetSignLs[[setC]]@tierMC, 2,
                            function(colVc) sum(colVc == "S"))) < 1) {
                
                if (!is.null(info.txtC))
                  cat("No significant features could be selected for the '",
                      setC,
                      "' dataset.\n",
                      sep = "")
                
              } else {
                
                if (!is.null(fig.pdfC)) {
                  
                  plot(msetSignLs[[setC]],
                       tierMaxC = plotTierMaxC,
                       typeC = "tier")
                  title(setC, line = 1, adj = 0, outer = TRUE)
                  
                  plot(msetSignLs[[setC]],
                       tierMaxC = plotTierMaxC,
                       typeC = "boxplot")
                  title(setC, line = 1, adj = 0, outer = TRUE)
                  
                }
                
                x <- MultiDataSet::add_eset(x,
                                            biosigner::getEset(msetSignLs[[setC]]),
                                            dataset.type = setC,
                                            GRanges = NA,
                                            overwrite = TRUE,
                                            warnings = FALSE)
                
              }
              
            }
            
            if (!is.null(fig.pdfC) && !is.na(fig.pdfC))
              dev.off()
            
            if (!is.null(info.txtC) && !is.na(info.txtC))
              sink()
            
            return(invisible(list(mset = x,
                                  signLs = msetSignLs)))
            
          })

