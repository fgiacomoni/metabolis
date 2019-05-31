#### metView (mset) ####

#' Viewing the MultiDataSet object.
#'
#' Provides numerical metrics and graphical overview of each dataset
#'
#' @param x MultiDataSet: multiple dataset to be explored
#' @param factorC Character: column name of pData to be used for sample color and
#' univariate test
#' @param fig.pdfC Character: Name of the file for the graphics; if NA, the
#' graphics are displayed on the screen;
#' @param poolAsPool1L Logical: should pool be included (as pool1) in the correlation
#' with the dilution factor?
#' @param poolCvN Numeric: threshold for the coefficient of variation of the pools
#' @param univTestC Character: if 'factorC' is specified, which test should be used
#' (see the univariate method)
#' @param univAdjustC Character: if 'factorC' is specified, which correction for
#' multiple testing should be used (see the univariate method)
#' @param univTreshN Character: if 'factorC' is specified, which threshold for
#' the corrected p-values should be used to reject the null hypothesis?
#' @param plotMainC Character: title; note that all 'plot' parameters are used only
#' for display
#' @param plotOverviewL Logical: should an overview of the number of samples and
#' variables in all datasets be barplotted?
#' @param plotSpanN Numeric: span parameter used in the loess trend estimation
#' @param plotSampIntensityC Character: function to be used to display the global
#' sample intensity
#' @param info.txtC Character: File name for the printed results (call to
#' 'sink()'); if NA (default), messages will be printed on the screen; if NULL,
#' no verbose will be generated
#' @return MultiDataSet including the computed sample and variable metrics
#' @rdname metView
#' @examples
#' prometMset <- metRead(system.file("extdata/promet",
#'                                   package = "metabolis"))
#'\dontrun{
#' prometMset <- metView(prometMset)
#' prometMset <- metView(prometMset, factorC = "gene")
#'}
#'
setMethod("metView", signature(x = "MultiDataSet"),
          function(x,
                   factorC = NA,
                   fig.pdfC = NA,
                   poolAsPool1L = FALSE,
                   poolCvN = 0.3,
                   univTestC = "nonparam",
                   univAdjustC = "BH",
                   univThreshN = 0.05,
                   plotMainC = "",
                   plotOverviewL = TRUE,
                   plotSpanN = 1,
                   plotSampIntensityC = c("median", "mean")[1],
                   info.txtC = NA) {
            
            if (!is.null(info.txtC) &&
                !is.na(info.txtC))
              sink(info.txtC, append = TRUE)
            
            infTxtC <- info.txtC
            if (!is.null(infTxtC))
              infTxtC <- NA
            
            if (!is.na(fig.pdfC))
              pdf(fig.pdfC)
            
            if (plotOverviewL)
              .msetBarplot(x, plotMainC)
            
            for (setC in names(x)) {
              
              if (!is.null(info.txtC))
                cat("\n\nViewing the '",
                    setC,
                    "' dataset:\n",
                    sep = "")
              
              ese <- x[[setC]]
              
              ese <- metView(ese,
                             factorC,
                             fig.pdfC = NA,
                             poolAsPool1L,
                             poolCvN,
                             univTestC,
                             univAdjustC,
                             univThreshN,
                             plotMainC = setC,
                             plotSpanN,
                             plotSampIntensityC,
                             info.txtC = infTxtC)
              
              x <- MultiDataSet::add_eset(x, ese,
                                          dataset.type = setC,
                                          GRanges = NA,
                                          overwrite = TRUE,
                                          warnings = FALSE)
              
            }
            
            if (!is.na(fig.pdfC))
              dev.off()
            
            if (!is.null(info.txtC) &&
                !is.na(info.txtC))
              sink()
            
            return(invisible(x))
            
          })

.msetBarplot <- function(x,
                         plotMainC = NA) {
  
  dimMN <- as.matrix(as.data.frame(Biobase::dims(x)))
  
  layout(matrix(2:1, nrow = 1),
         widths = c(0.7, 1.3))
  
  par(mar = c(2.6, 0.6, 2.1, 3.1))
  
  barMN <- barplot(dimMN[, ncol(dimMN):1],
                   beside = TRUE,
                   col = c("black", "grey"),
                   horiz = TRUE,
                   main = "",
                   names.arg = rep("", ncol(dimMN)),
                   log = "x",
                   xlab = "",
                   ylab = "")
  legend(10^(par("usr")[1] + 0.85 * diff(par("usr")[1:2])),
         1.05 * par("usr")[4],
         bty = "n",
         text.col = c("grey", "black"),
         legend =  c("samples", "variables"),
         xpd = TRUE)
  
  text(c(dimMN[, ncol(dimMN):1]),
       c(barMN),
       labels = format(c(dimMN[, ncol(dimMN):1]), big.mark = ","),
       pos = 4,
       adj = c(0, 0.5),
       xpd = TRUE)
  
  # Set names
  
  par(mar = c(2.6, 0.6, 2.1, 0.6))
  
  plot(c(0, 1),
       range(barMN),
       bty = "n",
       xlim = c(0, 5),
       ylim = par("usr")[3:4],
       xlab = "",
       ylab = "",
       type = "n",
       xaxt = "n",
       yaxt = "n",
       xaxs = "i",
       yaxs = "i")
  
  text(0,
       colMeans(barMN),
       cex = 1,
       labels = sapply(names(x)[length(names(x)):1],
                       function(namC) {
                         if (nchar(namC) > 20) {
                           paste0(substr(namC,
                                         1,
                                         20),
                                  ".")
                         } else {
                           namC
                         }
                       }),
       pos = 4)
  
  
  # End
  
  title(plotMainC, outer = TRUE, line = -1.5)
  
}

#### metView (eset) ####

#' Viewing the ExpressionSet dataset
#'
#' Provides numerical metrics and graphical overview of an ExpressionSet instance
#'
#' @param x ExpressionSet: dataset to be explored
#' @param factorC Character: column name of pData to be used for sample color and
#' univariate test
#' @param fig.pdfC Character: name of the pdf file for output
#' @param poolAsPool1L Logical: should pool be included (as pool1) in the correlation
#' with the dilution factor?
#' @param poolCvN Numeric: threshold for the coefficient of variation of the pools
#' @param univTestC Character: if 'factorC' is specified, which test should be used
#' (see the univariate method)
#' @param univAdjustC Character: if 'factorC' is specified, which correction for
#' multiple testing should be used (see the univariate method)
#' @param univTreshN Character: if 'factorC' is specified, which threshold for
#' the corrected p-values should be used to reject the null hypothesis?
#' @param plotMainC Character: title; note that all 'plot' parameters are used only for display
#' @param plotSpanN Numeric: span parameter used in the loess trend estimation
#' @param plotSampIntensityC Character: function to be used to display the global
#' sample intensity
#' @param info.txtC Character: File name for the printed results (call to
#' 'sink()'); if NA (default), messages will be printed on the screen; if NULL,
#' no verbose will be generated
#' @return ExpressionSet including the computed sample and variable metrics
#' @rdname metView
#' @examples
#' sacSet <- metRead(system.file("extdata/sacurine",
#'                               package = "metabolis"))
#' sacSet <- metView(sacSet)
#' sacSet <- metCorrect(sacSet)
#' sacSet <- metView(sacSet)
#' sacSet <- metTransform(sacSet)
#' sacSet <- metView(sacSet)
#' sacSet <- metView(sacSet, factorC = "gender")
setMethod("metView", signature(x = "ExpressionSet"),
          function(x,
                   factorC = NA,
                   fig.pdfC = NA,
                   poolAsPool1L = FALSE,
                   poolCvN = 0.3,
                   univTestC = "nonparam",
                   univAdjustC = "BH",
                   univThreshN = 0.05,
                   plotMainC = NA,
                   plotSpanN = 1,
                   plotSampIntensityC = c("median",
                                          "mean")[1],
                   info.txtC = NA) {
            
            if (!is.na(factorC) && !(factorC %in% colnames(Biobase::pData(x))))
              stop("'", factorC, "' was not found in the column names of pData(x).",
                   call. = FALSE)
            
            if (!is.null(info.txtC) &&
                !is.na(info.txtC))
              sink(info.txtC, append = TRUE)
            
            if (is.na(plotMainC))
              plotMainC <- Biobase::experimentData(x)@title
            
            ## Checks
            
            chkLs <- .checkW4Mformat(t(Biobase::exprs(x)),
                                     Biobase::pData(x),
                                     Biobase::fData(x))
            
            if (!chkLs[["chkL"]]) {
              stop("Sample and/or variable names do not match between your tables.",
                   call. = FALSE)
            } else if (chkLs[["ordL"]]) {
              Biobase::exprs(x) <- t(chkLs[["datMN"]])
            }
            
            ## Description
            
            if (!is.null(info.txtC)) {
              cat("\n\nData description:\n\n", sep = "")
              cat("observations:", ncol(Biobase::exprs(x)), "\n")
              cat("variables:", nrow(Biobase::exprs(x)), "\n")
              cat("missing:", sum(is.na(Biobase::exprs(x))), "\n")
              cat("0 values (%):",
                  sum(abs(Biobase::exprs(x)) < .Machine[["double.eps"]], na.rm = TRUE) / cumprod(dim(Biobase::exprs(x)))[2] * 100, "\n")
              cat("min:", min(Biobase::exprs(x), na.rm = TRUE), "\n")
              cat("mean:", signif(mean(Biobase::exprs(x), na.rm = TRUE), 2), "\n")
              cat("median:", signif(median(Biobase::exprs(x), na.rm = TRUE), 2), "\n")
              cat("max:", signif(max(Biobase::exprs(x), na.rm = TRUE), 2), "\n")
              
              if ("sampleType" %in% colnames(Biobase::pData(x))) {
                cat("\nSample types:\n", sep = "")
                print(table(Biobase::pData(x)[, "sampleType"]))
                cat("\n", sep = "")
              }
            }
            
            ## Sample metrics
            
            sampMetLs <- .sampleMetrics(x)
            
            x <- sampMetLs[["x"]]
            pcaLs <- sampMetLs[["pcaLs"]]
            
            ## Variable metrics
            
            x <- .variableMetrics(x,
                                  factorC,
                                  poolAsPool1L,
                                  univTestC,
                                  univAdjustC,
                                  univThreshN,
                                  info.txtC)
            
            ## Figure
            
            if (!is.na(fig.pdfC)) {
              pdf(fig.pdfC)
            }
            
            .plotMetrics(x,
                         factorC,
                         fig.pdfC,
                         pcaLs,
                         poolCvN,
                         univTestC,
                         univAdjustC,
                         univThreshN,
                         plotMainC,
                         plotSpanN,
                         plotSampIntensityC)
            
            if (!is.na(fig.pdfC))
              dev.off()
            
            ## End
            
            if (!is.null(info.txtC) &&
                !is.na(info.txtC))
              sink()
            
            return(invisible(x))
            
          })


.sampleMetrics <- function(x) {
  
  ## Hotelling: p-value associated to the distance from the center in the first PCA score plane
  
  pcaMod <- ropls::opls(t(Biobase::exprs(x)), predI = 2,
                        crossvalI = min(ncol(Biobase::exprs(x)), 7),
                        fig.pdfC = NULL, info.txtC = NULL)
  
  varRelVn <- ropls::getPcaVarVn(pcaMod) / nrow(Biobase::exprs(x)) ## for plotting
  
  pcaScoreMN <- ropls::getScoreMN(pcaMod)
  
  pdaDF <- Biobase::pData(x)
  
  pdaDF[, "pca_sco1"] <- pcaScoreMN[, 1]
  pdaDF[, "pca_sco2"] <- pcaScoreMN[, 2]
  
  invCovScoMN <- solve(cov(pcaScoreMN))
  
  n <- ncol(Biobase::exprs(x))
  hotN <- 2 * (n - 1) * (n^2 - 1) / (n^2 * (n - 2))
  
  hotPvaVn <- apply(pcaScoreMN,
                    1,
                    function(x)
                      1 - pf(1 / hotN * t(as.matrix(x)) %*% invCovScoMN %*% as.matrix(x), 2, n - 2))
  
  pdaDF[, "hotel_pval"] <- hotPvaVn
  
  pcaLs <- list(hotN = hotN,
                varRelVn = varRelVn)
  
  
  ## p-value associated to number of missing values
  
  missZscoVn <- .zscore(apply(t(Biobase::exprs(x)),
                              1,
                              function(rowVn) {
                                sum(is.na(rowVn))
                              }))
  
  pdaDF[, "miss_pval"] <- sapply(missZscoVn, function(zscoN) 2 * (1 - pnorm(abs(zscoN))))
  
  ## p-value associated to the deciles of the profiles
  
  deciMN <- t(as.matrix(apply(t(Biobase::exprs(x)),
                              1,
                              function(x) quantile(x, 0.1 * 1:9, na.rm = TRUE))))
  
  deciZscoMN <- apply(deciMN, 2, .zscore)
  
  deciZscoMaxVn <- apply(deciZscoMN, 1, function(rowVn) rowVn[which.max(abs(rowVn))])
  
  pdaDF[, "deci_pval"] <- sapply(deciZscoMaxVn, function(zscoN) 2 * (1 - pnorm(abs(zscoN))))
  
  Biobase::pData(x) <- pdaDF
  
  resLs <- list(x = x,
                pcaLs = pcaLs)
  
  return(resLs)
  
}


.variableMetrics <- function(x,
                             factorC, ## for the call to 'univariate'
                             poolAsPool1L,
                             uniTesC,
                             uniAdjC,
                             uniThrN,
                             info.txtC) { ## for the call to 'univariate'
  
  fdaDF <- Biobase::fData(x)
  
  tmpDF <- fdaDF ## some of the intermediate metrics will not be included in fData
  
  ## 'blank' observations
  
  if ("sampleType" %in% colnames(Biobase::pData(x)) && "blank" %in% Biobase::pData(x)[, "sampleType"]) {
    
    blkVl <- Biobase::pData(x)[, "sampleType"] == "blank"
    
    if (sum(blkVl) == 1)
      tmpDF[, "blank_mean"] <- t(Biobase::exprs(x))[blkVl, ]
    else
      tmpDF[, "blank_mean"] <- apply(t(Biobase::exprs(x))[blkVl, , drop = FALSE],
                                     2,
                                     function(varVn) mean(varVn, na.rm = TRUE))
    
    if (sum(blkVl) == 1)
      tmpDF[, "blank_sd"] <- rep(0, nrow(fdaDF))
    else
      tmpDF[, "blank_sd"] <- apply(t(Biobase::exprs(x))[blkVl, , drop = FALSE],
                                   2,
                                   function(varVn) sd(varVn, na.rm = TRUE))
    
    tmpDF[, "blank_CV"] <- tmpDF[, "blank_sd"] / tmpDF[, "blank_mean"]
    
  }
  
  ## 'sample' observations
  
  if ("sampleType" %in% colnames(Biobase::pData(x)) && "sample" %in% Biobase::pData(x)[, "sampleType"]) {
    
    samVl <- Biobase::pData(x)[, "sampleType"] == "sample"
    
    if (sum(samVl) == 1)
      tmpDF[, "sample_mean"] <- t(Biobase::exprs(x))[samVl, ]
    else
      tmpDF[, "sample_mean"] <- apply(t(Biobase::exprs(x))[samVl, , drop = FALSE], 2, function(varVn) mean(varVn, na.rm = TRUE))
    
    if (sum(samVl) == 1)
      tmpDF[, "sample_sd"] <- rep(0, nrow(fdaDF))
    else
      tmpDF[, "sample_sd"] <- apply(t(Biobase::exprs(x))[samVl, , drop = FALSE], 2, function(varVn) sd(varVn, na.rm = TRUE))
    
    tmpDF[, "sample_CV"] <- tmpDF[, "sample_sd"] / tmpDF[, "sample_mean"]
    
  }
  
  ## 'blank' mean / 'sample' mean ratio
  
  if (all(c("blank_mean", "sample_mean") %in% colnames(tmpDF))) {
    
    fdaDF[, "blankMean_over_sampleMean"] <- tmpDF[, "blank_mean"] / tmpDF[, "sample_mean"]
    
  }
  
  ## 'pool' observations
  
  if ("sampleType" %in% colnames(Biobase::pData(x)) && "pool" %in% Biobase::pData(x)[, "sampleType"]) {
    
    pooVl <- Biobase::pData(x)[, "sampleType"] == "pool"
    
    if (sum(pooVl) == 1)
      tmpDF[, "pool_mean"] <- t(Biobase::exprs(x))[pooVl, ]
    else
      tmpDF[, "pool_mean"] <- apply(t(Biobase::exprs(x))[pooVl, , drop = FALSE], 2, function(varVn) mean(varVn, na.rm = TRUE))
    
    if (sum(pooVl) == 1)
      tmpDF[, "pool_sd"] <- rep(0, nrow(fdaDF))
    else
      tmpDF[, "pool_sd"] <- apply(t(Biobase::exprs(x))[pooVl, , drop = FALSE], 2, function(varVn) sd(varVn, na.rm = TRUE))
    
    fdaDF[, "pool_CV"] <- tmpDF[, "pool_sd"] / tmpDF[, "pool_mean"]
    
  }
  
  ## 'pool' CV / 'sample' CV ratio
  
  if (all(c("pool_CV", "sample_CV") %in% colnames(Biobase::fData(x)))) {
    
    fdaDF[, "poolCV_over_sampleCV"] <- Biobase::fData(x)[, "pool_CV"] / tmpDF[, "sample_CV"]
    
  }
  
  ## 'pool' dilutions
  
  if ("sampleType" %in% colnames(Biobase::pData(x)) && any(grepl("pool.+", Biobase::pData(x)[, "sampleType"]))) {
    
    pooVi <- grep("pool.*", Biobase::pData(x)[, "sampleType"]) ## pool, pool2, pool4, poolInter, ...
    
    pooNamVc <- Biobase::pData(x)[pooVi, "sampleType"]
    
    if (poolAsPool1L) {
      
      pooNamVc[pooNamVc == "pool"] <- "pool1" ## 'pool' -> 'pool1'
      
    } else {
      
      pooVl <- pooNamVc == "pool"
      pooVi <- pooVi[!pooVl]
      pooNamVc <- pooNamVc[!pooVl]
      
    }
    
    pooDilVc <- gsub("pool", "", pooNamVc)
    
    pooDilVl <- sapply(pooDilVc, .allDigits)
    
    if (sum(pooDilVl)) {
      
      pooNamVc <- pooNamVc[pooDilVl]
      
      pooVi <- pooVi[pooDilVl]
      
      dilVn <- 1 / as.numeric(pooDilVc[pooDilVl])
      
      fdaDF[, "poolDil_cor"] <- apply(t(Biobase::exprs(x))[pooVi, , drop = FALSE], 2,
                                      function(varVn) cor(dilVn, varVn))
      
      fdaDF[, "poolDil_pval"] <- apply(t(Biobase::exprs(x))[pooVi, , drop = FALSE], 2,
                                       function(varVn) cor.test(dilVn, varVn)[["p.value"]])
      
    }
    
  }
  
  Biobase::fData(x) <- fdaDF
  
  ## Hypothesis testing
  
  if (!is.na(factorC)) {
    optWrnN <- options()$warn
    options(warn = -1)
    x <- metTest(x,
                 factorC,
                 testC = uniTesC,
                 adjustC = uniAdjC,
                 adjustThreshN = uniThrN,
                 fig.pdfC = NULL,
                 info.txtC = info.txtC)
    options(warn = optWrnN)
  }
  
  return(x)
  
}


.plotMetrics <- function(x,
                         factorC,
                         fig.pdfC,
                         pcaLs,
                         poolCvN,
                         uniTesC,
                         uniAdjC,
                         uniThrN,
                         plotMainC,
                         plotSpanN,
                         plotSampIntensityC) {
  
  ## Constants
  
  marLs <- list(tit = c(0.6, 1.1, 1.1, 0.6),
                drivol = c(3.5, 3.6, 4.1, 0.6),
                sca = c(0.6, 3.1, 4.1, 0.6),
                scavol = c(3.5, 3.1, 4.1, 0.6),
                ima = c(0.6, 2.6, 4.1, 0.9),
                imavol = c(3.5, 2.6, 4.1, 0.9),
                dri = c(3.5, 3.6, 1.1, 0.6),
                vol = c(3.5, 3.6, 1.1, 0.9),
                pca = c(3.5, 3.6, 1.1, 0.9))
  
  palHeaVc <- rev(rainbow(ceiling(256 * 1.5))[1:256])
  
  ## Functions
  
  .prettyAxis <- function(valVn,
                          lenN) {
    
    if (NA %in% valVn) {
      warning("NA in valVn")
      valVn <- as.vector(na.omit(valVn))
    }
    
    if (lenN < length(valVn))
      stop("The length of in vector must be inferior to the length of the length parameter.")
    
    if (length(valVn) < lenN)
      valVn <- seq(from = min(valVn), to = max(valVn),
                   length.out = lenN)
    
    preValVn <- pretty(valVn)
    
    preLabVn <- preAtVn <- c()
    
    for (n in 1:length(preValVn))
      if (min(valVn) < preValVn[n] && preValVn[n] < max(valVn)) {
        preLabVn <- c(preLabVn, preValVn[n])
        preAtVn <- c(preAtVn, which(abs(valVn - preValVn[n]) == min(abs(valVn - preValVn[n])))[1])
      }
    
    return(list(atVn = preAtVn,
                labVn = preLabVn))
    
  }
  
  ## Script
  
  facStaC <- NA
  
  if (!is.na(factorC)) {
    
    facTesLs <- .chooseTest(x, factorC, uniTesC)
    facLevVc <- facTesLs[["facLevVc"]]
    facTesC <- facTesLs[["tesC"]]
    facPfxC <- facTesLs[["varPfxC"]]
    
    if (facTesC %in% c("ttest",
                       "limma",
                       "wilcoxon")) {
      facStaC <- "dif"
    } else if (facTesC %in% c("pearson",
                              "spearman")) {
      facStaC <- "cor"
    }
    
    facAdjVn <- Biobase::fData(x)[, paste0(facPfxC, uniAdjC)]
    facTesMN <- cbind(adj = facAdjVn,
                      sig = as.numeric(facAdjVn <= uniThrN))
    if (!is.na(facStaC))
      facTesMN <- cbind(sta = Biobase::fData(x)[, paste0(facPfxC, facStaC)],
                        facTesMN)
    rownames(facTesMN) <- rownames(Biobase::fData(x))
    
  }
  
  opar <- par(font = 2,
              font.axis = 2,
              font.lab = 2,
              pch = 18)
  
  if (!is.na(facStaC)) {
    layVn <- c(4, 4, 2, 3,
               1, 1, 1, 5)
  } else
    layVn <- c(1, 2, 3, 3,
               4, 4, 4, 5)
  
  layout(matrix(layVn,
                byrow = TRUE,
                nrow = 2),
         heights = c(3.5, 3.5),
         widths = c(1.5, 1, 1, 3.5))
  
  # Colors
  
  if (!is.na(factorC)) {
    obsColVc <- ropls:::.plotColorF(as.vector(Biobase::pData(x)[, factorC]))[["colVc"]]
  } else if ("sampleType" %in% colnames(Biobase::pData(x))) {
    obsColVc <- .colorType(Biobase::pData(x)[, "sampleType"])
  } else
    obsColVc <- rep("black", nrow(Biobase::pData(x)))
  
  if (!is.na(facStaC)) {
    
    ## vol: Volcano Plot
    
    par(mar = marLs[["vol"]])
    
    .volcano(staVn = facTesMN[, "sta"],
             adjVn = facTesMN[, "adj"],
             cexN = 0.7,
             colC = "green4",
             adjC = uniAdjC,
             facC = factorC,
             labVc = rownames(facTesMN),
             levVc = facLevVc,
             staC = facStaC,
             tesC = facTesC,
             thrN = uniThrN)
    
  } else {
    
    ## tit: Title
    
    par(mar = marLs[["tit"]])
    plot(0:1, bty = "n", type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    text(1, 0.95, adj = 0, cex = 1.2, labels = plotMainC)
    text(1, 0.75, adj = 0, labels = paste0("NAs: ",
                                           round(length(which(is.na(c(t(Biobase::exprs(x)))))) / cumprod(dim(t(Biobase::exprs(x))))[2] * 100), "%"))
    text(1, 0.68, adj = 0, labels = paste0("0 values: ",
                                           round(sum(abs(t(Biobase::exprs(x))) < .Machine[["double.eps"]], na.rm = TRUE) / cumprod(dim(t(Biobase::exprs(x))))[2] * 100, 2), "%"))
    text(1, 0.61, adj = 0, labels = paste0("min: ", signif(min(t(Biobase::exprs(x)), na.rm = TRUE), 2)))
    text(1, 0.54, adj = 0, labels = paste0("median: ", signif(median(t(Biobase::exprs(x)), na.rm = TRUE), 2)))
    text(1, 0.47, adj = 0, labels = paste0("mean: ", signif(mean(t(Biobase::exprs(x)), na.rm = TRUE), 2)))
    text(1, 0.40, adj = 0, labels = paste0("max: ", signif(max(t(Biobase::exprs(x)), na.rm = TRUE), 2)))
    if ("sampleType" %in% colnames(Biobase::pData(x)) &&
        "pool" %in% Biobase::pData(x)[, "sampleType"])
      text(1,
           0.33,
           adj = 0,
           labels = paste0("CVpool<",
                           round(poolCvN * 100), "%: ",
                           round(sum(Biobase::fData(x)[, "pool_CV"] < poolCvN, na.rm = TRUE) / nrow(Biobase::exprs(x)) * 100),
                           "%"))
    
  }
  
  ## sca: Color scale
  
  if (!is.na(facStaC)) {
    par(mar = marLs[["scavol"]])
  } else
    par(mar = marLs[["sca"]])
  
  ylimVn <- c(0, 256)
  ybottomVn <- 0:255
  ytopVn <- 1:256
  
  plot(x = 0,
       y = 0,
       font.axis = 2,
       font.lab = 2,
       type = "n",
       xlim = c(0, 1),
       ylim = ylimVn,
       xlab = "",
       ylab = "",
       xaxs = "i",
       yaxs = "i",
       xaxt = "n",
       yaxt = "n")
  
  rect(xleft = 0,
       ybottom = ybottomVn,
       xright = 1,
       ytop = ytopVn,
       col = palHeaVc,
       border = NA)
  
  eval(parse(text = paste0("axis(at = .prettyAxis(c(ifelse(min(t(Biobase::exprs(x)), na.rm = TRUE) == -Inf, yes = 0, no = min(t(Biobase::exprs(x)), na.rm = TRUE)) , max(t(Biobase::exprs(x)), na.rm = TRUE)), 256)$atVn,
                           font = 2,
                           font.axis = 2,
                           labels = .prettyAxis(c(ifelse(min(t(Biobase::exprs(x)), na.rm = TRUE) == -Inf, yes = 0, no = min(t(Biobase::exprs(x)), na.rm = TRUE)), max(t(Biobase::exprs(x)), na.rm = TRUE)), 256)$labVn,
                           las = 1,
                           lwd = 2,
                           lwd.ticks = 2,
                           side = 2,
                           xpd = TRUE)")))
  
  arrows(par("usr")[1],
         par("usr")[4],
         par("usr")[1],
         par("usr")[3],
         code = 0,
         lwd = 2,
         xpd = TRUE)
  
  box(lwd = 2)
  
  ## ima: Image
  
  if (!is.na(facStaC)) {
    par(mar = marLs[["imavol"]])
  } else
    par(mar = marLs[["ima"]])
  
  imaMN <- Biobase::exprs(x)[, rev(1:ncol(Biobase::exprs(x))), drop = FALSE]
  
  image(x = 1:nrow(imaMN),
        y = 1:ncol(imaMN),
        z = imaMN,
        col = palHeaVc,
        font.axis = 2,
        font.lab = 2,
        xaxt = "n",
        yaxt = "n",
        xlab = "",
        ylab = "")
  
  if (length(rownames(t(Biobase::exprs(x)))) == 0) {
    rowNamVc <- rep("", times = ncol(Biobase::exprs(x)))
  } else
    rowNamVc <- rownames(t(Biobase::exprs(x)))
  
  if (length(colnames(t(Biobase::exprs(x)))) == 0) {
    colNamVc <- rep("", times = nrow(Biobase::exprs(x)))
  } else
    colNamVc <- colnames(t(Biobase::exprs(x)))
  
  xlaVc <- paste(paste(rep("[", 2),
                       c(1, nrow(imaMN)),
                       rep("] ", 2),
                       sep = ""),
                 rep("\n", times = 2),
                 c(colNamVc[1], tail(colNamVc, 1)),
                 sep = "")
  
  for (k in 1:2)
    axis(side = 3,
         hadj = c(0, 1)[k],
         at = c(1, nrow(imaMN))[k],
         cex = 0.8,
         font = 2,
         labels = xlaVc[k],
         line = -0.5,
         tick = FALSE)
  
  
  ylaVc <- paste(paste(rep("[", times = 2),
                       c(ncol(imaMN), 1),
                       rep("]", times = 2),
                       sep = ""),
                 rep("\n", times = 2),
                 c(tail(rowNamVc, 1), rowNamVc[1]),
                 sep = "")
  
  for (k in 1:2)
    axis(side = 2,
         at = c(1, ncol(imaMN))[k],
         cex = 0.8,
         font = 2,
         hadj = c(0, 1)[k],
         labels = ylaVc[k],
         las = 0,
         line = -0.5,
         lty = "blank",
         tick = FALSE)
  
  box(lwd = 2)
  
  ## dri: Analytical drift
  
  if (!is.na(facStaC)) {
    par(mar = marLs[["drivol"]])
  } else
    par(mar = marLs[["dri"]])
  
  ## ordering
  
  driProMN <- t(Biobase::exprs(x))
  driObsDF <- Biobase::pData(x)
  
  driObsDF[, "ordIniVi"] <- 1:nrow(driProMN)
  
  if ("injectionOrder" %in% colnames(driObsDF)) {
    ordNamC <- "Injection Order"
    if ("batch" %in% colnames(driObsDF)) {
      ordVi <- order(driObsDF[, "batch"],
                     driObsDF[, "injectionOrder"])
    } else
      ordVi <- order(driObsDF[, "injectionOrder"])
  } else {
    ordNamC <- "Samples"
    ordVi <- 1:nrow(driProMN)
  }
  
  driProMN <- driProMN[ordVi, ]
  driObsDF <- driObsDF[ordVi, ]
  driColVc <- obsColVc[ordVi]
  
  if ("batch" %in% colnames(driObsDF))
    batTab <- table(driObsDF[, "batch"])
  
  driSumVn <- eval(parse(text = paste0("apply(driProMN, 1, function(obsVn) ", plotSampIntensityC, "(obsVn, na.rm = TRUE))")))
  
  plot(driSumVn,
       col = driColVc,
       pch = 18,
       ylim = range(driProMN, na.rm = TRUE),
       type = "n",
       xaxs = "i",
       xlab = "",
       ylab = "")
  
  for (samI in 1:nrow(driProMN))
    boxplot(driProMN[samI, ],
            at = samI,
            add = TRUE)
  
  points(driSumVn,
         col = driColVc,
         pch = 18)
  
  mtext(ordNamC,
        cex = 0.7,
        line = 2,
        side = 1)
  
  mtext(paste0(toupper(substr(plotSampIntensityC, 1, 1)), substr(plotSampIntensityC, 2, nchar(plotSampIntensityC)), " of variable intensities"),
        cex = 0.7,
        line = 2,
        side = 2)
  
  if ("batch" %in% colnames(driObsDF)) {
    
    abline(v = cumsum(batTab) + 0.5,
           col = "red")
    
    mtext(names(batTab),
          at = batTab / 2 + c(0, cumsum(batTab[-length(batTab)])),
          cex = 0.7)
    
    for (batC in names(batTab)) {
      
      batSeqVi <- which(driObsDF[, "batch"] == batC)
      
      if ("sampleType" %in% colnames(driObsDF)) {
        batSamVi <- intersect(batSeqVi,
                              grep("sample", driObsDF[, "sampleType"]))
      } else
        batSamVi <- batSeqVi
      
      lines(batSeqVi,
            .loess(driSumVn, batSamVi, batSeqVi, plotSpanN),
            col = .colorType("sample"))
      
      if ("sampleType" %in% colnames(driObsDF) &&
          "pool" %in% driObsDF[, "sampleType"]) {
        
        batPooVi <- intersect(batSeqVi,
                              grep("^pool$", driObsDF[, "sampleType"]))
        
        lines(batSeqVi,
              .loess(driSumVn, batPooVi, batSeqVi, plotSpanN),
              col = .colorType("pool"))
        
      }
      
    }
    
  } else {
    
    batSeqVi <- 1:nrow(driObsDF)
    
    if ("sampleType" %in% colnames(driObsDF)) {
      batSamVi <- intersect(batSeqVi,
                            grep("sample", driObsDF[, "sampleType"]))
    } else
      batSamVi <- batSeqVi
    
    lines(batSeqVi,
          .loess(driSumVn, batSamVi, batSeqVi, plotSpanN),
          col = .colorType("sample"))
    
    if ("sampleType" %in% colnames(driObsDF) &&
        "pool" %in% driObsDF[, "sampleType"]) {
      
      batPooVi <- intersect(batSeqVi,
                            grep("^pool$", driObsDF[, "sampleType"]))
      
      lines(batSeqVi,
            .loess(driSumVn, batPooVi, batSeqVi, plotSpanN),
            col = .colorType("pool"))
      
      
    }
    
  }
  
  if (!is.na(facStaC))
    mtext(plotMainC,
          line = 2,
          side = 3)
  
  ## pca: PCA and Hotelling ellipse
  
  par(mar = marLs[["pca"]])
  
  plot(Biobase::pData(x)[, c("pca_sco1", "pca_sco2")],
       type = "n",
       xlab = "",
       ylab = "",
       xlim = range(Biobase::pData(x)[, "pca_sco1"]) * 1.1)
  mtext(paste("t1 (", round(pcaLs[["varRelVn"]][1] * 100), "%)", sep = ""),
        cex = 0.7,
        line = 2,
        side = 1)
  mtext(paste("t2 (", round(pcaLs[["varRelVn"]][2] * 100), "%)", sep = ""),
        cex = 0.7,
        las = 0,
        line = 2,
        side = 2)
  abline(h = 0, lty = "dashed")
  abline(v = 0, lty = "dashed")
  radVn <- seq(0, 2 * pi, length.out = 100)
  
  hotFisN <- pcaLs[["hotN"]] * qf(1 - 0.05, 2, ncol(Biobase::exprs(x)) - 2)
  
  lines(sqrt(var(Biobase::pData(x)[, "pca_sco1"]) * hotFisN) * cos(radVn),
        sqrt(var(Biobase::pData(x)[, "pca_sco2"]) * hotFisN) * sin(radVn))
  
  text(Biobase::pData(x)[, "pca_sco1"],
       Biobase::pData(x)[, "pca_sco2"],
       cex = 0.7,
       col = obsColVc,
       labels = Biobase::sampleNames(x))
  
  if (!is.na(factorC)) {
    ropls:::.plotLegendF(as.vector(Biobase::pData(x)[, factorC]),
                         as.matrix(Biobase::pData(x)[, c("pca_sco1", "pca_sco2")]))
  } else if ("sampleType" %in% colnames(Biobase::pData(x))) {
    obsColVuc <- obsColVc[sort(unique(names(obsColVc)))]
    legOrdVc <- c("blank", paste0("pool", 8:1), "pool", "other", "sample")
    obsColVuc <- obsColVuc[legOrdVc[legOrdVc %in% names(obsColVuc)]]
    
    text(rep(par("usr")[1], times = length(obsColVuc)),
         par("usr")[3] + (0.97 - length(obsColVuc) * 0.03 + 1:length(obsColVuc) * 0.03) * diff(par("usr")[3:4]),
         col = obsColVuc,
         font = 2,
         labels = names(obsColVuc),
         pos = 4)
  }
  
  par(mfrow = c(1, 1))
  par(opar)
  
}


.allDigits <- function(string) { ## from the Hmisc package (all.digits)
  k <- length(string)
  result <- logical(k)
  for (i in 1:k) {
    st <- string[i]
    ls <- nchar(st)
    ex <- substring(st, 1:ls, 1:ls)
    result[i] <- all(match(ex, c("0", "1", "2", "3", "4",
                                 "5", "6", "7", "8", "9"), nomatch = 0) > 0)
  }
  result
}


.colorType <- function(typVc) {
  
  typColVc <- c(sample = "green4",
                pool = "red",
                pool1 = RColorBrewer::brewer.pal(9, "Reds")[7],
                pool2 = RColorBrewer::brewer.pal(9, "Reds")[5],
                pool4 = RColorBrewer::brewer.pal(9, "Reds")[3],
                pool8 = RColorBrewer::brewer.pal(9, "Reds")[2],
                blank = "black",
                other = "yellow")
  
  typVc[!(typVc %in% setdiff(names(typColVc), "other"))] <- "other"
  
  typColVc[typVc]
  
}


.loess <- function(datVn, qcaVi, preVi, spanN) {
  
  if (length(qcaVi) < 5) {
    
    return(predict(lm(datVn[qcaVi] ~ qcaVi),
                   newdata = data.frame(qcaVi = preVi)))
    
  } else {
    
    return(predict(loess(datVn[qcaVi] ~ qcaVi,
                         control = loess.control(surface = "direct"),
                         span = spanN),
                   newdata = data.frame(qcaVi = preVi)))
    
  }
  
  ## Note:
  ##  the surface = 'direct' argument allows extrapolation
  
}


.volcano <- function(staVn,
                     adjVn,
                     adjC = "BH",
                     cexN = 1,
                     colC = "green4",
                     facC,
                     labVc = NULL,
                     levVc,
                     maiC = NA,
                     staC = c("dif", "cor")[1],
                     tesC = NA,
                     thrN = 0.05) {
  
  sigVl <- adjVn <= thrN
  sigVl[is.na(sigVl)] <- FALSE
  
  maiC <- paste0(ifelse(is.na(maiC),
                        "",
                        paste0(maiC, " ")),
                 tesC,
                 " '",
                 facC,
                 "', ",
                 adjC,
                 " (sig: ",
                 format(sum(sigVl),
                        big.mark = ","),
                 ")")
  
  plot(staVn,
       -log10(adjVn),
       xlab = "",
       ylab = "",
       main = maiC,
       type = "n")
  mtext(paste0(staC, " (", facC, ")"),
        cex = cexN,
        line = ifelse(cexN == 1, 2.5, 2),
        side = 1)
  mtext(paste0("-log10(", adjC, ")"),
        cex = cexN,
        las = 0,
        line = ifelse(cexN == 1, 2.5, 2),
        side = 2)
  
  abline(h = -log10(thrN), col = colC,
         lwd = 1)
  
  abline(v = axTicks(1), col = "gray",
         lwd = 1)
  
  if (staC == "dif") {
    mtext(levVc[1],
          adj = 0.1,
          line = -1.5,
          side = 3)
    mtext(levVc[2],
          adj = 0.9,
          line = -1.5,
          side = 3)
  }
  
  points(staVn,
         -log10(adjVn),
         col = ifelse(sigVl, colC, "black"))
  
  if (sum(sigVl) && !is.null(labVc)) {
    text(staVn[sigVl],
         -log10(adjVn)[sigVl],
         labels = labVc[sigVl],
         pos = ifelse(staVn[sigVl] > 0, 4, 2),
         col = colC)
  }
  
  
}


.zscore <- function(x) {
  sdxN <- sd(x, na.rm = TRUE)
  if (sdxN < .Machine[["double.eps"]]) {
    return(rep(0, length(x)))
  } else
    return((x - mean(x, na.rm = TRUE)) / sdxN)
}

#' Venn Diagram
#'
#' Venn Diagram of a list containing up to 5 vectors
#'
#' @param inputLs list containing up to 5 vectors
#' @param file.tiffC file name for the figure, if set to
#' NULL (default), the figure is displayed on the screen
#' @param lwdN size of the circle lines
#' @param mainC title (by default, the name of the matrixMN variable will be used)
#' @param palVc color palette to be used
#' @param subVc subtitle
#' @return No output.
#' @export
#' @examples
#' print("hello")
vennF <- function(inputLs,
                  file.tiffC = NULL,
                  lwdN = 2,
                  mainC = NA,
                  palVc = RColorBrewer::brewer.pal(9, "Set1")[1:5],
                  subC = "") {
  
  if (is.na(mainC))
    mainC <- gsub("MN", "",
                  gsub("Ls", "",
                       deparse(substitute(inputLs))))
  
  if (class(inputLs) != "list")
    stop("'inputLs' must be a list for Venn plot", call. = FALSE)
  
  if (length(inputLs) > 5)
    stop("'inputLs' list must be of maximum length 5 for Venn plot", call. = FALSE)
  
  ## palVc <- c(c100 = brewer.pal(9, "Reds")[7],
  ##            c110 = brewer.pal(9, "Purples")[4],
  ##            c010 = brewer.pal(9, "Blues")[7],
  ##            c011 = brewer.pal(12, "Set3")[2],
  ##            c001 = brewer.pal(9, "Greens")[7],
  ##            c101 = brewer.pal(9, "Oranges")[4],
  ##            c111 = "white")
  ## palVc <- palVc[c("c100", "c010", "c001",
  ##                  brewer.pal(9, )[7],
  ##                  brewer.pal(9, )[7])]
  
  catN <- length(inputLs)
  
  if (length(palVc) < catN)
    stop("'palVc' must contain at least 'length(inputLs)' colors", call. = FALSE)
  
  futile.logger::flog.threshold(futile.logger::ERROR,
                                name = "VennDiagramLogger")  
  
  if (catN <= 3) { 
    ven = VennDiagram::venn.diagram(x = inputLs,
                                    
                                    alpha = 0.8,
                                    
                                    col = palVc[1:catN],
                                    
                                    cat.cex = 1.7, ## 2.5
                                    cat.col = palVc[1:catN],
                                    
                                    ## additional argument for catN <= 3
                                    cat.dist = c(rep(ifelse(catN == 2,
                                                            0.03, 0.05), 2),
                                                 ifelse(catN == 3, 0.04, 0.05),
                                                 rep(0.04, 2))[1:catN],
                                    
                                    cex = 1.7, ## 3
                                    
                                    fill = palVc[1:catN],
                                    
                                    label.col = "white",
                                    lwd = lwdN, ## 4
                                    
                                    main = mainC,
                                    main.cex = 1.7,
                                    main.pos = c(0.5, 1.05),
                                    
                                    margin = 0.01,
                                    
                                    filename = file.tiffC,
                                    cat.fontfamily = "sans",
                                    cat.fontface = "bold",
                                    fontfamily = "sans",
                                    fontface = "bold",
                                    main.fontfamily = "sans",
                                    main.fontface = "bold",
                                    sub = subC,
                                    sub.cex = 1,
                                    sub.fontfamily = "sans",
                                    sub.pos = c(0.5, 1.0))
  } else {
    ven <- VennDiagram::venn.diagram(x = inputLs,
                                     
                                     alpha = 0.8,
                                     
                                     col = palVc[1:catN],
                                     
                                     cat.cex = 1.7, ## 2.5
                                     cat.col = palVc[1:catN],
                                     
                                     cex = 1.7, ## 3
                                     
                                     fill = palVc[1:catN],
                                     
                                     label.col = "white",
                                     lwd = lwdN, ## 4
                                     
                                     main = mainC,
                                     main.cex = 1.7,
                                     main.pos = c(0.5, 1.05),
                                     
                                     margin = 0.01,
                                     
                                     filename = file.tiffC,
                                     cat.fontfamily = "sans",
                                     cat.fontface = "bold",
                                     fontfamily = "sans",
                                     fontface = "bold",
                                     main.fontfamily = "sans",
                                     main.fontface = "bold",
                                     sub = subC,
                                     sub.cex = 1,
                                     sub.fontfamily = "sans",
                                     sub.pos = c(0.5, 1.0))
  }
  
  if (is.null(file.tiffC)) {
    grid::grid.newpage()
    grid::grid.draw(ven)
  }
  
  ## venMN <- venn(inputVnMNLs, show.plot = FALSE,
  ##               intersections = FALSE)
  ## venMN <- venMN[rownames(venMN) != "000", ]
  ## venVn <- venMN[, "num"]
  ## venMN <- venMN[, colnames(venMN) != "num"]
  ## names(venVn) <- as.character(sapply(names(venVn),
  ##                                     function(codC) {
  ##                                         paste(colnames(venMN)[as.logical(venMN[codC, ])], collapse = ":")
  ##                                     }))
  
  return(invisible(ven))
  
}
