#### opls (mset) ####

#' PCA and (O)PLS-DA modeling
#'
#' PCA and (O)PLS-DA modeling
#'
#' @param x An S4 object of class \code{MultiDataSet}
#' @param y Character: Name of the fData column to be used for (O)PLS(-DA) modeling;
#' (set to NULL in case of PCA)
#' @param fig.pdfC Character: Name of the file for the graphics from the significant features;
#' if NA, the graphics are displayed on the screen; if NULL, no plot is generated
#' @param info.txtC Character: File name for the printed results (call to 'sink()')
#' @param ... additional parameters for the ropls::opls modeling
#' @return list of 'mset', the \code{MultiDataSet} including the scores, loadings,
#' VIP, etc. in the corresponding pData and fData when a model could be built
#' (>=1 component) and was significant (pQ2 <= 5%), and 'oplsLs', the list
#' of 'opls' outputs
#' @rdname opls
#' @export
#' @examples
#' prometMset <- readSet(system.file("extdata/promet", package="metabolis"))
#' prometMset <- opls(prometMset, "gene", predI = 1, orthoI = 1)
setMethod("opls", signature(x = "MultiDataSet"),
          function(x,
                   y = NULL,
                   fig.pdfC = NA,
                   info.txtC = NULL,
                   ...) {
            
            if (!is.null(info.txtC) && !is.na(info.txtC))
              sink(info.txtC, append = TRUE)
            
            infTxtC <- info.txtC
            if (!is.null(infTxtC))
              infTxtC <- NA
            
            if (!is.null(fig.pdfC) && !is.na(fig.pdfC))
              pdf(fig.pdfC)
            
            msetOplsLs <- vector(mode = "list",
                                 length = length(names(x)))
            names(msetOplsLs) <- names(x)
            
            for (setC in names(x)) {
              
              if (!is.null(info.txtC))
                cat("\nBuilding model for the '",
                    setC,
                    "' dataset:\n",
                    sep = "")
              
              setOpls <- tryCatch(ropls::opls(x[[setC]],
                                              y = y,
                                              info.txtC = infTxtC,
                                              fig.pdfC = NULL,
                                              ...),
                                  error = function(e) NULL)
              
              if (is.null(setOpls)) {
                
                if (!is.null(info.txtC))
                  cat("No model could be built for the '",
                      setC,
                      "' dataset.\n",
                      sep = "")
                
                msetOplsLs[[setC]] <- NULL
                
              } else if (grepl("PLS", setOpls@typeC) &&
                         ropls::getSummaryDF(setOpls)["Total", "pQ2"] > 0.05) {
                
                if (!is.null(info.txtC))                  
                  cat("No model was included for the '",
                      setC,
                      "' dataset because pQ2 was above 5%.\n",
                      sep = "")
                
                msetOplsLs[[setC]] <- NULL
                
              } else {
                
                msetOplsLs[[setC]] <- setOpls
                
                if (!is.null(fig.pdfC))
                  plot(msetOplsLs[[setC]],
                       fig.pdfC = fig.pdfC)
                
                x <- MultiDataSet::add_eset(x,
                                            ropls::getEset(msetOplsLs[[setC]]),
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
                                  oplsLs = msetOplsLs)))
            
          })