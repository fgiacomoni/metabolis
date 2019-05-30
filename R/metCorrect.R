#### metCorrect (mset) ####

#' Signal drift and batch effect correction
#'
#' Signal drift and batch effect correction
#'
#'
#' @param x An S4 object of class \code{MultiDataSet}
#' @param referenceSampeTypeC Character: sample type to be used as reference for the correction (as indicated in the 'colnameSampleType' column from the
#' pData(x); e.g. 'pool')
#' @param colnameBatchC Character: name of the column from pData(x) containing the batch information (encoded as characters)
#' @param colnameInjectionOrder Character: name of the column from pData(x) containing the injection order information (encoded as numerics)
#' @param colnameSampleTypeC Character:  name of the column from pData(x) containing the sample type information (encoded as characters) 
#' @param spanN Numeric: smoothing parameter for the loess regression; between 0 and 1; (default set to 1)
#' @param fig.pdfC Character: Name of the file for the graphics from the significant features;
#' if NA, the graphics are displayed on the screen;
#' @param info.txtC Character: File name for the printed results (call to
#' 'sink()'); if NA (default), messages will be printed on the screen; if NULL,
#' no verbose will be generated
#' @return \code{MultiDataSet} including the adjusted p-values in fData (or the vector of adjusted p-values when returnAdjustOnlyL is TRUE)
#' @rdname metCorrect
#' @export
#' @examples
#' # to be provided
setMethod("metCorrect", signature(x = "MultiDataSet"),
          function(x,
                   referenceSampleTypeC = c("pool", "sample")[1],
                   colnameBatchC = "batch",
                   colnameInjectionOrderC = "injectionOrder",
                   colnameSampleTypeC = "sampleType",
                   spanN = 1,
                   fig.pdfC = NA,
                   info.txtC = NA) {
            
            if (!is.null(info.txtC) &&
                !is.na(info.txtC))
              sink(info.txtC, append = TRUE)
            
            infTxtC <- info.txtC
            if (!is.null(infTxtC))
              infTxtC <- NA
            
            for (setC in names(x)) {
              
              if (!is.null(info.txtC))
                cat("\nCorrecting the '",
                    setC,
                    "' dataset:\n",
                    sep = "")
              
              ese <- metabolis::metCorrect(x[[setC]],
                                           referenceSampleTypeC = referenceSampleTypeC,
                                           colnameBatchC = colnameBatchC,
                                           colnameInjectionOrderC = colnameInjectionOrderC,
                                           colnameSampleTypeC = colnameSampleTypeC,
                                           spanN = spanN,
                                           fig.pdfC = fig.pdfC,
                                           info.txtC = infTxtC)
              
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


#### metCorrect (eset) ####

#' Signal drift and batch effect correction
#'
#' Signal drift and batch effect correction
#'
#' @param x An S4 object of class \code{ExpressionSet}
#' @param referenceSampeTypeC Character: sample type to be used as reference for the correction (as indicated in the 'colnameSampleType' column from the
#' pData(x); e.g. 'pool')
#' @param colnameBatchC Character: name of the column from pData(x) containing the batch information (encoded as characters)
#' @param colnameInjectionOrder Character: name of the column from pData(x) containing the injection order information (encoded as numerics)
#' @param colnameSampleTypeC Character:  name of the column from pData(x) containing the sample type information (encoded as characters) 
#' @param spanN Numeric: smoothing parameter for the loess regression; between 0 and 1; (default set to 1)
#' @param fig.pdfC Character: Name of the file for the graphics from the significant features;
#' if NA, the graphics are displayed on the screen;
#' @param info.txtC Character: File name for the printed results (call to
#' 'sink()'); if NA (default), messages will be printed on the screen; if NULL,
#' no verbose will be generated
#' @return \code{ExpressionSet} including the adjusted p-values in fData (or the vector of adjusted p-values when returnAdjustOnlyL is TRUE)
#' @rdname metCorrect
#' @export
#' @examples
#' sacSet <- metRead(system.file("extdata/sacurine", package = "metabolis"))
#' sacSet <- metCorrect(sacSet)
setMethod("metCorrect", signature(x = "ExpressionSet"),
          function(x,
                   referenceSampleTypeC = c("pool", "sample")[1],
                   colnameBatchC = "batch",
                   colnameInjectionOrderC = "injectionOrder",
                   colnameSampleTypeC = "sampleType",
                   spanN = 1,
                   fig.pdfC = NA,
                   info.txtC = NA) {
            
            if (!is.null(info.txtC) &&
                !is.na(info.txtC))
              sink(info.txtC, append = TRUE)
            
            x <- .metCorrect(x,
                             refC = referenceSampleTypeC,
                             colBatC = colnameBatchC,
                             colInjC = colnameInjectionOrderC,
                             colSamC = colnameSampleTypeC,
                             spnN = spanN,
                             fig.pdfC = fig.pdfC)
            
            validObject(x)
            
            if (!is.null(info.txtC) &&
                !is.na(info.txtC))
              sink()
            
            return(invisible(x))
            
          })

.metCorrect <- function(eset,
                        refC,
                        colBatC,
                        colInjC,
                        colSamC,
                        spnN,
                        fig.pdfC) {
  
  rawMN <- t(Biobase::exprs(eset))
  samDF <- Biobase::pData(eset)
  varDF <- Biobase::fData(eset)
  
  ## checking

  if (sum(grepl(refC, samDF[, colSamC])) == 0)
    stop("No '", refC, "' reference sample type found in the '", colSamC, "' column of the sampleMetadata.")
  
  refMN <- rawMN[samDF[, colSamC] == refC, ]
  refNasZerVl <- apply(refMN, 2,
                       function(refVn)
                         all(sapply(refVn,
                                    function(refN) {is.na(refN) || refN == 0})))
  
  if (sum(refNasZerVl)) {
    
    cat("The following variables have 'NA' or 0 values in all reference samples; they will be removed from the data:\n", sep = "")
    rawMN <- rawMN[, !refNasZerVl, drop = FALSE]
    varDF <- varDF[!refNasZerVl, , drop = FALSE]
    
  }
  
  ## Computation
  
  ## ordering (batch and injection order)
  
  samDF[, "ordIniVi"] <- 1:nrow(rawMN)
  ordBatInjVi <- order(samDF[, colBatC],
                       samDF[, colInjC])
  rawMN <- rawMN[ordBatInjVi, ]
  samDF <- samDF[ordBatInjVi, ]
  
  ## signal drift and batch-effect correction
  
  nrmMN <- .shiftBatchCorrect(rawMN,
                              samDF,
                              refC,
                              spnN,
                              colBatC,
                              colSamC)
  
  ## figure
  
  if (!is.null(fig.pdfC)) {
    
    if (!is.na(fig.pdfC))
      pdf(fig.pdfC,
          onefile = TRUE,
          width = 11,
          height = 7)
    
    .plotBatch(rawMN, samDF, spnN, colBatC, colSamC)
    .plotBatch(nrmMN, samDF, spnN, colBatC, colSamC)
    
    if (!is.na(fig.pdfC))
      dev.off()
    
  }
  
  ## returning to initial order
  
  ordIniVi <- order(samDF[, "ordIniVi"])
  nrmMN <- nrmMN[ordIniVi, ]
  samDF <- samDF[ordIniVi, ]
  samDF <- samDF[, colnames(samDF) != "ordIniVi", drop = FALSE]
  
  ## ending
  
  Biobase::exprs(eset) <- t(nrmMN)
  Biobase::pData(eset) <- samDF
  Biobase::fData(eset) <- varDF
  
  ## returning
  
  eset
  
}

.loess <- function(datVn, qcaVi, preVi, span) {
  
  if (length(qcaVi) < 5) {
    
    return(predict(lm(datVn[qcaVi] ~ qcaVi),
                   newdata = data.frame(qcaVi = preVi)))
    
  } else {
    
    return(predict(loess(datVn[qcaVi] ~ qcaVi,
                         control = loess.control(surface = "direct"),
                         span = span),
                   newdata = data.frame(qcaVi = preVi)))
    
  }
  
  ## Note:
  ##  the surface = 'direct' argument allows extrapolation
  
}

.plotBatch <- function(datMN,
                       samDF.arg,
                       spnN.arg,
                       colBatC.arg,
                       colSamC.arg) {
  
  maiC <- switch(gsub("MN", "", deparse(substitute(datMN))),
                 raw = "Raw",
                 nrm = "Normalized")
  
  colVc <- c(samp = "green4",
             biol = "green4",
             pool = "red",
             blan = "black",
             other = "yellow")
  
  par(font = 2, font.axis = 2, font.lab = 2, lwd = 2, pch = 18)
  
  layout(matrix(c(1, 1, 2, 3), nrow = 2),
         widths = c(0.7, 0.3))
  
  obsNamVc <- rownames(datMN)
  
  obsColVc <- sapply(substr(samDF.arg[, colSamC.arg], 1, 4),
                     function(typC)
                       ifelse(typC %in% names(colVc), colVc[typC], colVc["other"]))
  
  ## Graphic 1: Sum of intensities for each sample
  
  par(mar = c(3.6, 3.6, 3.1, 0.6))
  
  batTab <- table(samDF.arg[, colBatC.arg])
  
  sumVn <- rowSums(datMN, na.rm = TRUE)
  
  plot(sumVn,
       cex = 1.2,
       col = obsColVc,
       pch = 18,
       xaxs = "i",
       xlab = "",
       ylab = "")
  
  mtext("Injection order",
        line = 2.2,
        side = 1)
  mtext("Sum of variable intensities",
        line = 2.2,
        side = 2)
  
  mtext(maiC, cex = 1.2, line = 1.5, side = 3)
  
  abline(v = cumsum(batTab) + 0.5,
         col = "red")
  
  mtext(names(batTab),
        at = batTab / 2 + c(0, cumsum(batTab[-length(batTab)])))
  
  obsColVuc <- obsColVc[sort(unique(names(obsColVc)))]
  
  text(rep(batTab[1], times = length(obsColVuc)),
       par("usr")[3] + (0.97 - length(obsColVuc) * 0.03 + 1:length(obsColVuc) * 0.03) * diff(par("usr")[3:4]),
       col = obsColVuc,
       font = 2,
       labels = names(obsColVuc),
       pos = 2)
  
  for (batC in names(batTab)) {
    
    batSeqVi <- which(samDF.arg[, colBatC.arg] == batC)
    batPooVi <- intersect(batSeqVi,
                          which(samDF.arg[, colSamC.arg] == "pool"))
    batSamVi <- intersect(batSeqVi,
                          which(samDF.arg[, colSamC.arg] == "sample"))
    if (length(batPooVi))
      lines(batSeqVi,
            .loess(sumVn, batPooVi, batSeqVi, spnN.arg),
            col = colVc["pool"])
    if (length(batSamVi))
      lines(batSeqVi,
            .loess(sumVn, batSamVi, batSeqVi, spnN.arg),
            col = colVc["samp"])
    
  }
  
  ## Graphics 2 and 3 (right): PCA score plots of components 1-4
  
  radVn <- seq(0, 2 * pi, length.out = 100)
  epsN <- .Machine[["double.eps"]] ## [1] 2.22e-16
  
  pcaMN <- datMN
  
  if (any(is.na(pcaMN))) {
    minN <- min(pcaMN, na.rm = TRUE)
    pcaMN[is.na(pcaMN)] <- minN
  }
  
  pcaLs <- ropls::opls(pcaMN, predI = 4, algoC = "svd", fig.pdfC = NULL, info.txtC = NULL)
  tMN <- ropls::getScoreMN(pcaLs)
  vRelVn <- pcaLs@modelDF[, "R2X"]
  
  n <- nrow(tMN)
  hotN <- 2 * (n - 1) * (n^2 - 1) / (n^2 * (n - 2))
  
  hotFisN <- hotN * qf(0.95, 2, n - 2)
  
  pcsLs <- list(c(1, 2), c(3, 4))
  
  par(mar = c(3.6, 3.6, 0.6, 1.1))
  
  for (pcsN in 1:length(pcsLs)) {
    
    pcsVn <- pcsLs[[pcsN]]
    
    tcsMN <- tMN[, pcsVn]
    
    micMN <- solve(cov(tcsMN))
    
    n <- nrow(tMN)
    hotN <- 2 * (n - 1) * (n^2 - 1) / (n^2 * (n - 2))
    
    hotFisN <- hotN * qf(0.95, 2, n - 2)
    
    hotVn <- apply(tcsMN,
                   1,
                   function(x) 1 - pf(1 / hotN * t(as.matrix(x)) %*% micMN %*% as.matrix(x), 2, n - 2))
    
    obsHotVi <- which(hotVn < 0.05)
    
    xLabC <- paste("t",
                   pcsVn[1],
                   "(",
                   round(vRelVn[pcsVn[1]] * 100),
                   "%)",
                   sep = "")
    
    yLabC <- paste("t",
                   pcsVn[2],
                   "(",
                   round(vRelVn[pcsVn[2]] * 100),
                   "%)",
                   sep = "")
    
    xLimVn <- c(-1, 1) * max(sqrt(var(tcsMN[, 1]) * hotFisN), max(abs(tcsMN[, 1])))
    yLimVn <- c(-1, 1) * max(sqrt(var(tcsMN[, 2]) * hotFisN), max(abs(tcsMN[, 2])))
    
    plot(tcsMN,
         main = "",
         type = "n",
         xlab = "",
         ylab = "",
         xlim = xLimVn,
         ylim = yLimVn)
    
    mtext(xLabC,
          line = 2.2,
          side = 1)
    mtext(yLabC,
          line = 2.2,
          side = 2)
    
    par(lwd = 1)
    
    abline(v = axTicks(1),
           col = "grey")
    
    abline(h = axTicks(2),
           col = "grey")
    
    abline(v = 0)
    abline(h = 0)
    
    lines(sqrt(var(tcsMN[, 1]) * hotFisN) * cos(radVn),
          sqrt(var(tcsMN[, 2]) * hotFisN) * sin(radVn))
    
    points(tcsMN,
           col = obsColVc,
           pch = 18)
    
    if (length(obsHotVi))
      text(tcsMN[obsHotVi, 1],
           tcsMN[obsHotVi, 2],
           col = obsColVc[obsHotVi],
           labels = obsNamVc[obsHotVi],
           pos = 3)
    
  } ## for(pcsN in 1:length(pcsLs)) {
  
  return(invisible(list(sumVn = sumVn,
                        tcsMN = tcsMN)))
  
} ## plotBatchF

.shiftBatchCorrect <- function(rawMN.arg,
                               samDF.arg,
                               refC.arg,
                               spnN.arg,
                               colBatC.arg,
                               colSamC.arg) {
  
  cat("\nReference observations are: ", refC.arg, "\n")
  
  ## computing median off all pools (or samples) for each variable
  
  refMeaVn <- apply(rawMN.arg[samDF.arg[, colSamC.arg] == refC.arg, ],
                    2,
                    function(feaRefVn) mean(feaRefVn, na.rm = TRUE))
  
  ## splitting data and sample metadata from each batch
  
  batRawLs <- split(as.data.frame(rawMN.arg),
                    f = samDF.arg[, colBatC.arg])
  batRawLs <- lapply(batRawLs, function(inpDF) as.matrix(inpDF))
  
  batSamLs <- split(as.data.frame(samDF.arg),
                    f = samDF.arg[, colBatC.arg])
  
  ## checking extrapolation: are there pools at the first and last observations of each batch

    pooExtML <- matrix(FALSE, nrow = 2, ncol = length(batRawLs),
                       dimnames = list(c("first", "last"), names(batRawLs)))
    for (batC in names(batSamLs)) {
      batSamTypVc <- batSamLs[[batC]][, colSamC.arg]
      pooExtML["first", batC] <- head(batSamTypVc, 1) == refC.arg
      pooExtML["last", batC] <- tail(batSamTypVc, 1) == refC.arg
    }
    if (!all(c(pooExtML))) {
      cat("\nWarning: Reference samples are missing at the first and/or last position of the following batches:\n")
      pooExtBatVi <- which(!apply(pooExtML, 2, all))
      for (i in 1:length(pooExtBatVi))
        cat(names(pooExtBatVi)[i], ": ",
            paste(rownames(pooExtML)[!pooExtML[, pooExtBatVi[i]]], collapse = ", "), "\n", sep = "")
      cat("Extrapolating loess fits for these batches may result in inaccurate modeling!\n")
    }
  
  ## normalizing
  
  nrmMN <- NULL ## normalized data matrix to be computed
  
  cat("\nProcessing batch:")
  
  for (batC in names(batRawLs)) { ## processing each batch individually
    
    cat("\n", batC)
    
    batRawMN <- batRawLs[[batC]]
    batSamDF <- batSamLs[[batC]]
    
    batAllVi <- 1:nrow(batRawMN)
    
    batRefVi <- which(batSamDF[, colSamC.arg] == refC.arg)
    
    if (length(batRefVi) < 5)
      cat("\nWarning: less than 5 '", refC.arg,
          "'; linear regression will be performed instead of loess regression for this batch\n",
          sep = "")
    
    ## prediction of the loess fit
    
    batLoeMN <- apply(batRawMN,
                      2,
                      function(rawVn) .loess(rawVn,
                                             batRefVi,
                                             batAllVi,
                                             spnN.arg))
    
    ## normalization
    
    batLoeMN[batLoeMN <= 0] <- NA
    
    batNrmMN <- batRawMN / batLoeMN
    
    nrmMN <- rbind(nrmMN,
                   batNrmMN)
    
  }
  
  cat("\n")
  
  nrmMN <- sweep(nrmMN, MARGIN = 2, STATS = refMeaVn, FUN = "*")
  
  return(nrmMN)
  
} ## shiftBatchCorrectF
