#### metTest (mset) ####

#' Univariate hypothesis testing
#'
#' Provides numerical metrics and graphical overview of an ExpressionSet instance
#'
#' @param x An S4 object of class \code{MultiDataSet}
#' @param factorC Character: Factor of interest (name of a column from the
#' pData(x))
#' @param testC Character: One of the 6 available hypothesis tests can be selected
#' (either 'ttest', 'limma', 'wilcoxon', 'anova', 'kruskal', 'pearson', or 'spearman');
#' if set to 'param' (or 'nonparam' [default]), the appropriate (non-)parametric test will be used
#' @param adjustC Character: Name of the method for correction of multiple testing
#' (the p.adjust function is used)
#' @param adjustThreshN Numeric: Threshold for (corrected) p-values
#' @param fig.pdfC Character: Name of the file for the graphics from the significant features;
#' if NA, the graphics are displayed on the screen; if NULL, no plot is generated
#' @param info.txtC Character: File name for the printed results (call to
#' 'sink()'); if NA (default), messages will be printed on the screen; if NULL,
#' no verbose will be generated
#' @return \code{MultiDataSet} including the adjusted p-values in the fData data frames
#' @rdname metTest
#' @export
#' @examples
#' prometMset <- readSet(system.file("extdata/promet", package="metabolis"))
#' prometMset <- metTest(prometMset, "gene")
setMethod("metTest", signature(x = "MultiDataSet"),
          function(x,
                   factorC,
                   testC = c("param", "nonparam", "ttest", "limma", "wilcoxon", "anova", "kruskal", "pearson", "spearman")[2],
                   adjustC = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")[5],
                   adjustThreshN = 0.05,
                   fig.pdfC = NA,
                   info.txtC = NA) {
            
            if (!is.null(info.txtC) && !is.na(info.txtC))
              sink(info.txtC, append = TRUE)
            
            infTxtC <- info.txtC
            if (!is.null(infTxtC))
              infTxtC <- NA
            
            if (!is.null(fig.pdfC) && !is.na(fig.pdfC))
              pdf(fig.pdfC)
            
            figPdfC <- fig.pdfC
            if (!is.null(figPdfC))
              figPdfC <- NA            
            
            for (setC in names(x)) {
              
              if (!is.null(info.txtC))
                cat("\nHypothesis testing for the '",
                    setC,
                    "' dataset:\n",
                    sep = "")
              
              ese <- x[[setC]]
              
              if (!(factorC %in% colnames(Biobase::pData(ese)))) {
                
                if (!is.null(info.txtC))
                  cat("The '",
                      factorC,
                      "' factor was not found in the pData columns from the '",
                      setC,
                      "' dataset.\nNo test will be computed for this dataset.\n",
                      sep = "")
                
              } else {
                
                ese <- metTest(ese,
                               factorC,
                               testC = testC,
                               adjustC = adjustC,
                               adjustThreshN = adjustThreshN,
                               plotMainC = setC,
                               fig.pdfC = figPdfC,
                               info.txtC = infTxtC)
                
                x <- MultiDataSet::add_eset(x,
                                            ese,
                                            dataset.type = setC,
                                            GRanges = NA,
                                            overwrite = TRUE,
                                            warnings = FALSE)
                
              }
              
            }
            
            if (!is.null(fig.pdfC) && !is.na(fig.pdfC))
              dev.off()
            
            if (!is.null(info.txtC) &&
                !is.na(info.txtC))
              sink()
            
            return(invisible(x))
            
          })


#### metTest (eset) ####

#' Univariate hypothesis testing
#'
#' Provides numerical metrics and graphical overview of an ExpressionSet instance
#'
#' @param x An S4 object of class \code{ExpressionSet}
#' @param factorC Character: Factor of interest (name of a column from the
#' pData(x))
#' @param testC Character: One of the 6 available hypothesis tests can be selected
#' (either 'ttest', 'limma', 'wilcoxon', 'anova', 'kruskal', 'pearson', or 'spearman');
#' if set to 'param' (or 'nonparam' [default]), the appropriate (non-)parametric test will be used
#' @param adjustC Character: Name of the method for correction of multiple testing
#' (the p.adjust function is used)
#' @param adjustThreshN Numeric: Threshold for (corrected) p-values
#' @param fig.pdfC Character: Name of the file for the graphics from the significant features;
#' if NA, the graphics are displayed on the screen; if NULL, no plot is generated
#' @param plotMainC character: title to be displayed on the first plot
#' @param info.txtC Character: File name for the printed results (call to
#' 'sink()'); if NA (default), messages will be printed on the screen; if NULL,
#' no verbose will be generated
#' @return \code{ExpressionSet} including the adjusted p-values in fData (or the vector of adjusted p-values when returnAdjustOnlyL is TRUE)
#' @rdname metTest
#' @export
#' @examples
#' proSet <- readSet(system.file("extdata/promet/proteo", package="metabolis"))
#' proSet <- metTest(proSet, "gene", fig.pdfC = NULL)
#' head(fData(proSet))
#'\dontrun{
#' proSet <- metTest(proSet, "gene", fig.pdfC = NA)
#' proSet <- metTest(proSet, "gene", fig.pdfC = "proSet_metTest-gene.pdf")
#'}
setMethod("metTest", signature(x = "ExpressionSet"),
          function(x,
                   factorC,
                   testC = c("param", "nonparam", "ttest", "limma", "wilcoxon", "anova", "kruskal", "pearson", "spearman")[2],
                   adjustC = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")[5],
                   adjustThreshN = 0.05,
                   fig.pdfC = NA,
                   plotMainC = NA,
                   info.txtC = NA) {
            
            if (!is.null(info.txtC) &&
                !is.na(info.txtC))
              sink(info.txtC, append = TRUE)
            
            if (!(factorC %in% colnames(Biobase::pData(x))))
              stop("'", factorC, "' was not found in the column names of pData(x).",
                   call. = FALSE)
            
            if (is.na(plotMainC))
              plotMainC <- Biobase::experimentData(x)@title
            
            fdaDF <- Biobase::fData(x)
            
            ## Getting the response (either a factor or a numeric)
            
            tesLs <- metabolis:::.chooseTest(x, factorC, testC)
            
            facFcVn <- tesLs[["facFcVn"]]
            facLevVc <- tesLs[["facLevVc"]]
            testC <- tesLs[["tesC"]]
            varPfxC <- tesLs[["varPfxC"]]
            
            if (!is.null(info.txtC))
              cat("\nPerforming '", testC, "'\n", sep = "")
            
            ## Two-level or quantitative
            
            if (testC %in% c("ttest", "limma", "wilcoxon",
                             "pearson", "spearman")) {
              
              staC <- ifelse(testC %in% c("ttest", "limma", "wilcoxon"),
                             "dif",
                             "cor")
              
              switch(testC,
                     ttest = {
                       staF <- function(y) diff(tapply(y,
                                                       facFcVn,
                                                       function(x)
                                                         mean(x, na.rm = TRUE)))
                       tesF <- function(y) t.test(y ~ facFcVn)[["p.value"]]
                     },
                     limma = {
                       staF <- function(y) diff(tapply(y,
                                                       facFcVn,
                                                       function(x)
                                                         mean(x, na.rm = TRUE)))
                     },
                     wilcoxon = {
                       staF <- function(y) diff(tapply(y,
                                                       facFcVn,
                                                       function(x)
                                                         median(x, na.rm = TRUE)))
                       tesF <- function(y) wilcox.test(y ~ facFcVn)[["p.value"]]
                     },
                     pearson = {
                       staF <- function(y) cor(facFcVn, y, method = "pearson",
                                               use = "pairwise.complete.obs")
                       tesF <- function(y) cor.test(facFcVn, y, method = "pearson",
                                                    use = "pairwise.complete.obs")[["p.value"]]
                     },
                     spearman = {
                       staF <- function(y) cor(facFcVn, y, method = "spearman",
                                               use = "pairwise.complete.obs")
                       tesF <- function(y) cor.test(facFcVn, y, method = "spearman",
                                                    use = "pairwise.complete.obs")[["p.value"]]
                     })
              
              staVn <- apply(Biobase::exprs(x), 1, staF)
              
              optWrnN <- options()$warn
              options(warn = -1)
              if (testC != "limma") {
                adjVn <- p.adjust(apply(Biobase::exprs(x), 1, tesF),
                                  method = adjustC)
              } else {
                dsgMN <- cbind(lev2 = rep(1, Biobase::dims(x)["Samples", ]),
                               lev1.lev2 = as.numeric(tesLs[["facFcVn"]] == tesLs[["facLevVc"]][1]))
                rownames(dsgMN) <- Biobase::sampleNames(x)
                limFit <- limma::lmFit(x, dsgMN)
                limFit <- limma::eBayes(limFit)
                adjVn <- p.adjust(limFit$p.value[, "lev1.lev2"],
                                  method = adjustC)
              }
              options(warn = optWrnN)
              
              sigVn <- as.numeric(adjVn <= adjustThreshN)
              
              resMN <- cbind(staVn, adjVn, sigVn)
              rownames(resMN) <- Biobase::featureNames(x)
              colnames(resMN) <- paste0(varPfxC, c(staC, adjustC, "sig"))
              
              for (colC in colnames(resMN))
                fdaDF[, colC] <- resMN[, colC]
              
              if (sum(sigVn, na.rm = TRUE) > 0) {
                resSigMN <- resMN[sigVn > 0, 1:2, drop = FALSE]
                resSigMN <- resSigMN[order(adjVn[sigVn > 0]), , drop = FALSE]
              } else {
                resSigMN <- NULL
              }
              
              ## graphic
              
              if (!is.null(fig.pdfC)) {
                
                if (!is.na(fig.pdfC))
                  pdf(fig.pdfC, onefile = TRUE)
                
                .volcano(staVn = staVn,
                         adjVn = adjVn,
                         adjC = adjustC,
                         cexN = 1,
                         colC = "green4",
                         facC = factorC,
                         labVc = Biobase::featureNames(x),
                         levVc = facLevVc,
                         maiC = plotMainC,
                         staC = staC,
                         tesC = testC,
                         thrN = adjustThreshN)
                
                ## if (!is.na(fig.pdfC)) {
                
                varVi <- which(sigVn > 0)
                
                if (testC %in% c("ttest", "limma", "wilcoxon")) {
                  
                  facVc <- as.character(facFcVn)
                  names(facVc) <- Biobase::sampleNames(x)
                  
                  for (varI in varVi) {
                    
                    varC <- Biobase::featureNames(x)[varI]
                    
                    .boxplot(facFcVn,
                             t(Biobase::exprs(x))[, varI],
                             paste0(varC, " (", adjustC, " = ", signif(adjVn[varI], 2), ")"),
                             facVc)
                    
                  }
                  
                } else {## pearson or spearman
                  
                  for (varI in varVi) {
                    
                    varC <- Biobase::featureNames(x)[varI]
                    
                    mod <- lm(t(Biobase::exprs(x))[, varI] ~  facFcVn)
                    
                    plot(facFcVn, t(Biobase::exprs(x))[, varI],
                         xlab = factorC,
                         ylab = "",
                         pch = 18,
                         main = paste0(varC, " (", adjustC, " = ", signif(adjVn[varI], 2), ", R2 = ", signif(summary(mod)$r.squared, 2), ")"))
                    
                    abline(mod, col = "red")
                    
                  }
                  
                }
                
                if (!is.na(fig.pdfC))
                  dev.off()
                
                ## }
                
              }
              
              ## More than two levels
              
            } else if (testC == "anova") {
              
              ## getting the names of the pairwise comparisons 'class1Vclass2'
              prwVc <- rownames(TukeyHSD(aov(t(Biobase::exprs(x))[, 1] ~ facFcVn))[["facFcVn"]])
              
              prwVc <- gsub("-", ".", prwVc, fixed = TRUE) ## 2016-08-05: '-' character in dataframe column names seems not to be converted to "." by write.table on ubuntu R-3.3.1
              
              ## omnibus and post-hoc tests
              
              aovMN <- t(apply(t(Biobase::exprs(x)),
                               2,
                               function(varVn) {
                                 
                                 aovMod <- aov(varVn ~ facFcVn)
                                 pvaN <- summary(aovMod)[[1]][1, "Pr(>F)"]
                                 hsdMN <- TukeyHSD(aovMod)[["facFcVn"]]
                                 c(pvaN, c(hsdMN[, c("diff", "p adj")]))
                                 
                               }))
              
              difVi <- 1:length(prwVc) + 1
              
              ## difference of the means for each pairwise comparison
              
              difMN <- aovMN[, difVi]
              colnames(difMN) <- paste0(varPfxC, prwVc, "_dif")
              
              ## correction for multiple testing
              
              aovMN <- aovMN[, -difVi, drop = FALSE]
              aovMN <- apply(aovMN, 2, function(pvaVn) p.adjust(pvaVn, method = adjustC))
              
              ## significance coding (0 = not significant, 1 = significant)
              
              adjVn <- aovMN[, 1]
              sigVn <-  as.numeric(adjVn < adjustThreshN)
              
              aovMN <- aovMN[, -1, drop = FALSE]
              colnames(aovMN) <- paste0(varPfxC, prwVc, "_", adjustC)
              
              aovSigMN <- aovMN < adjustThreshN
              mode(aovSigMN) <- "numeric"
              colnames(aovSigMN) <- paste0(varPfxC, prwVc, "_sig")
              
              ## final aggregated table
              
              resMN <- cbind(adjVn, sigVn, difMN, aovMN, aovSigMN)
              rownames(resMN) <- Biobase::featureNames(x)
              colnames(resMN)[1:2] <- paste0(varPfxC, c(adjustC, "sig"))
              
              for (colC in colnames(resMN))
                fdaDF[, colC] <- resMN[, colC]
              
              if (sum(sigVn, na.rm = TRUE) > 0) {
                resSigMN <- resMN[sigVn > 0, 1, drop = FALSE]
                resSigMN <- resSigMN[order(adjVn[sigVn > 0]), , drop = FALSE]
              } else {
                resSigMN <- NULL
              }
              
              ## graphic
              
              if (!is.null(fig.pdfC)) {
                
                if (!is.na(fig.pdfC))
                  pdf(fig.pdfC, onefile = TRUE)
                
                for (varI in 1:nrow(Biobase::fData(x))) {
                  
                  if (sum(aovSigMN[varI, ]) > 0) {
                    
                    varC <- Biobase::featureNames(x)[varI]
                    
                    boxplot(t(exprs(x))[, varI] ~ facFcVn,
                            main = paste0(varC, " (", adjustC, " = ", signif(adjVn[varI], 2), ")"))
                    
                    for (prwI in 1:length(prwVc)) {
                      
                      if (aovSigMN[varI, paste0(varPfxC, prwVc[prwI], "_sig")] == 1) {
                        
                        claVc <- unlist(strsplit(prwVc[prwI], ".", fixed = TRUE))
                        aovClaVl <- facFcVn %in% claVc
                        aovFc <- facFcVn[aovClaVl, drop = TRUE]
                        aovVc <- as.character(aovFc)
                        names(aovVc) <- Biobase::sampleNames(x)[aovClaVl]
                        .boxplot(aovFc,
                                 t(Biobase::exprs(x))[aovClaVl, varI],
                                 paste0(varC, " (", adjustC, " = ", signif(aovMN[varI, paste0(varPfxC, prwVc[prwI], "_", adjustC)], 2), ")"),
                                 aovVc)
                        
                      }
                      
                    }
                    
                  }
                  
                }
                
                if (!is.na(fig.pdfC))
                  dev.off()
                
              }
              
            } else if (testC == "kruskal") {
              
              ## getting the names of the pairwise comparisons 'class1.class2'
              
              optWrnN <- options()$warn
              options(warn = -1)
              nemMN <- PMCMRplus::kwAllPairsNemenyiTest(t(Biobase::exprs(x))[, 1], facFcVn, "Tukey")[["p.value"]]
              options(warn = optWrnN)
              
              nemVl <- c(lower.tri(nemMN, diag = TRUE))
              nemClaMC <- cbind(rownames(nemMN)[c(row(nemMN))][nemVl],
                                colnames(nemMN)[c(col(nemMN))][nemVl])
              nemNamVc <- paste0(nemClaMC[, 1], ".", nemClaMC[, 2])
              pfxNemVc <- paste0(varPfxC, nemNamVc)
              
              ## omnibus and post-hoc tests
              
              optWrnN <- options()$warn
              options(warn = -1)
              nemMN <- t(apply(Biobase::exprs(x), 1, function(varVn) {
                
                pvaN <- kruskal.test(varVn ~ facFcVn)[["p.value"]]
                varNemMN <- PMCMRplus::kwAllPairsNemenyiTest(varVn, facFcVn, "Tukey")[["p.value"]]
                c(pvaN, c(varNemMN))
                
              }))
              options(warn = optWrnN)
              
              ## correction for multiple testing
              
              nemMN <- apply(nemMN, 2,
                             function(pvaVn) p.adjust(pvaVn, method = adjustC))
              adjVn <- nemMN[, 1]
              sigVn <- as.numeric(adjVn < adjustThreshN)
              nemMN <- nemMN[, c(FALSE, nemVl)]
              colnames(nemMN) <- paste0(pfxNemVc, "_", adjustC)
              
              ## significance coding (0 = not significant, 1 = significant)
              
              nemSigMN <- nemMN < adjustThreshN
              mode(nemSigMN) <- "numeric"
              colnames(nemSigMN) <- paste0(pfxNemVc, "_sig")
              nemSigMN[is.na(nemSigMN)] <- 0
              
              ## difference of the medians for each pairwise comparison
              
              difMN <- sapply(1:nrow(nemClaMC), function(prwI) {
                prwVc <- nemClaMC[prwI, ]
                prwVi <- which(facFcVn %in% prwVc)
                prwFacFc <- factor(as.character(facFcVn)[prwVi], levels = prwVc)
                apply(Biobase::exprs(x)[, prwVi], 1,
                      function(varVn)
                        -diff(as.numeric(tapply(varVn, prwFacFc,
                                                function(x)
                                                  median(x, na.rm = TRUE)))))
              })
              colnames(difMN) <- gsub("_sig", "_dif", colnames(nemSigMN))
              
              ## final aggregated table
              
              resMN <- cbind(adjVn, sigVn, difMN, nemMN, nemSigMN)
              rownames(resMN) <- Biobase::featureNames(x)
              colnames(resMN)[1:2] <- paste0(varPfxC, c(adjustC, "sig"))
              
              for (colC in colnames(resMN))
                fdaDF[, colC] <- resMN[, colC]
              
              if (sum(sigVn, na.rm = TRUE) > 0) {
                resSigMN <- resMN[sigVn > 0, 1, drop = FALSE]
                resSigMN <- resSigMN[order(adjVn[sigVn > 0]), , drop = FALSE]
              } else {
                resSigMN <- NULL
              }
              
              ## graphic
              
              if (!is.null(fig.pdfC)) {
                
                if (!is.na(fig.pdfC))
                  pdf(fig.pdfC, onefile = TRUE)
                
                for (varI in 1:nrow(Biobase::fData(x))) {
                  
                  if (sum(nemSigMN[varI, ]) > 0) {
                    
                    varC <- Biobase::featureNames(x)[varI]
                    
                    boxplot(t(Biobase::exprs(x))[, varI] ~ facFcVn,
                            main = paste0(varC, " (", adjustC, " = ", signif(adjVn[varI], 2), ")"))
                    
                    for (nemI in 1:length(nemNamVc)) {
                      
                      if (nemSigMN[varI, paste0(varPfxC, nemNamVc[nemI], "_sig")] == 1) {
                        
                        nemClaVc <- nemClaMC[nemI, ]
                        nemClaVl <- facFcVn %in% nemClaVc
                        nemFc <- facFcVn[nemClaVl, drop = TRUE]
                        nemVc <- as.character(nemFc)
                        names(nemVc) <- Biobase::sampleNames(x)[nemClaVl]
                        .boxplot(nemFc,
                                 t(Biobase::exprs(x))[nemClaVl, varI],
                                 paste0(varC, " (", adjustC, " = ", signif(nemMN[varI, paste0(varPfxC, nemNamVc[nemI], "_", adjustC)], 2), ")"),
                                 nemVc)
                        
                      }
                      
                    }
                    
                  }
                  
                }
                
                if (!is.na(fig.pdfC))
                  dev.off()
                
              }
              
            }
            
            ## Finalizing
            
            if (!is.null(info.txtC)) {
              if (!is.null(resSigMN)) {
                
                cat("\n", nrow(resSigMN), " variable",
                    ifelse(nrow(resSigMN) > 1, "s", ""),
                    " (", round(nrow(resSigMN) / length(sigVn) * 100), "%) ",
                    ifelse(nrow(resSigMN) > 1, "were", "was"),
                    " found significant at the ", adjustThreshN,
                    " level (after the '", adjustC, "' correction).\n",
                    ifelse(nrow(resSigMN) == 1,
                           "It is",
                           ifelse(nrow(resSigMN) <= 15,
                                  "They are",
                                  "The first 15 are")),
                    " displayed below",
                    ifelse(nrow(resSigMN) > 1,
                           " (sorted by increasing corrected p-values)",
                           ""),
                    ":\n",
                    sep = "")
                print(resSigMN[1:min(15, nrow(resSigMN)), , drop = FALSE])
                
              } else
                cat("\nNo significant variable found at the selected ", adjustThreshN, " level\n", sep = "")
            }
            
            Biobase::fData(x) <- fdaDF
            
            validObject(x)
            
            if (!is.null(info.txtC) &&
                !is.na(info.txtC))
              sink()
            
            return(invisible(x))
            
          })

.boxplot <- function(xFc,
                     yVn,
                     maiC,
                     xVc) {
  
  boxLs <- boxplot(yVn ~  xFc,
                   main = maiC)
  
  outVn <- boxLs[["out"]]
  
  if (length(outVn)) {
    
    for (outI in 1:length(outVn)) {
      levI <- which(levels(xFc) == xVc[names(outVn)[outI]])
      text(levI,
           outVn[outI],
           labels = names(outVn)[outI],
           pos = ifelse(levI == 2, 2, 4))
    }
    
  }
  
}

.chooseTest <- function(eset, facC, tesC){
  
  if (mode(Biobase::pData(eset)[, facC]) == "character") {
    facFcVn <- factor(Biobase::pData(eset)[, facC])
    facLevVc <- levels(facFcVn)
    facLevN <- nlevels(facFcVn)
    if (tesC %in% c("param", "nonparam")) {
      if (facLevN == 1) {
        stop("'", facC, "' has only one level.", call. = FALSE)
      } else if (facLevN == 2) {
        tesC <- switch(tesC, param = "ttest", nonparam = "wilcoxon")
      } else {
        tesC <- switch(tesC, param = "anova", nonparam = "kruskal")
      }
    }
  } else {
    facFcVn <- Biobase::pData(eset)[, facC]
    facLevVc <- NA
    if (tesC %in% c("param", "nonparam"))
      tesC <- switch(tesC, param = "pearson", nonparam = "spearman")
  }
  
  varPfxC <- paste0(make.names(facC), "_", tesC, "_")
  if (tesC %in% c("ttest", "limma", "wilcoxon"))
    varPfxC <- paste0(varPfxC, paste(facLevVc, collapse = "."), "_")
  
  return(list(facFcVn = facFcVn,
              facLevVc = facLevVc,
              tesC = tesC,
              varPfxC = varPfxC))
  
}
