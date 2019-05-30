#### metHeatmap (mset) ####

#' Heatmap
#'
#' Heatmap of a MultiDataSet instance
#'
#' @param x an S4 object of class \code{MultiDataSet}
#' @param fig.pdfC Character: Name of the file for the graphics from the significant features;
#' if NA, the graphics are displayed on the screen;
#' @param info.txtC Character: File name for the printed results (call to
#' 'sink()'); if NA (default), messages will be printed on the screen; if NULL,
#' no verbose will be generated
#' @param ... additional parameters (see the documentation)
#' @return \code{MultiDataSet} including the column with the
#' cluster numbers in pData and fData 
#' @rdname metHeatmap
#' @export
#' @examples
#' prometMset <- metRead(system.file("extdata/promet/",
#' package="metabolis"))
#' prometMset <- metHeatmap(prometMset)
setMethod("metHeatmap", signature(x = "MultiDataSet"),
          function(x,
                   fig.pdfC = NA,
                   info.txtC = NA,
                   ...) {

            if (!is.null(info.txtC) &&
                !is.na(info.txtC))
              sink(info.txtC, append = TRUE)
            
            infTxtC <- info.txtC
            if (!is.null(infTxtC))
              infTxtC <- NA
                        
            if (!is.na(fig.pdfC))
              pdf(fig.pdfC)
            
            for (setC in names(x)) {
              
              if (!is.null(info.txtC))
                cat("\nHeatmap of the '",
                    setC,
                    "' dataset:\n",
                    sep = "")
              
              ese <- x[[setC]]
              
              ese <- metHeatmap(ese,
                                 fig.pdfC = NA,
                                 info.txtC = infTxtC,
                                 ...)
              
              title(setC,
                    line = -1,
                    adj = 1,
                    outer = TRUE)
              
              x <- MultiDataSet::add_eset(x,
                                          ese,
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


#### metHeatmap (eset) ####

#' Heatmap
#'
#' Heatmap of an ExpressionSet instance
#'
#' @param x An S4 object of class \code{ExpressionSet}
#' @param dissymC Character: [hclust]
#' @param correlC Character: correlation coefficient (in case
#' '1-cor' or '1-abs(cor)' are selected as dissymilarity) 
#' @param aggloC charcter: agglomeration method
#' @param clustVi tupple of integers: number of sample and variable clusters, respectively
#' @param fig.pdfC Character: Name of the file for the graphics from the significant features;
#' if NA, the graphics are displayed on the screen;
#' @param plotCexVn vector of numerics:
#' @param plotPaletteC character: color palette
#' @param plotScaleL logical: scaling (mean-centering and unit
#' variance scaling) to enhance contrast (for plotting only)
#' @param info.txtC Character: File name for the printed results (call to
#' 'sink()'); if NA (default), messages will be printed on the screen; if NULL,
#' no verbose will be generated
#' @return \code{ExpressionSet} including the adjusted p-values in fData (or the vector of adjusted p-values when returnAdjustOnlyL is TRUE)
#' @rdname metHeatmap
#' @export
#' @examples
#' metSet <- metRead(system.file("extdata/promet/metabo",
#' package="metabolis"))
#' metSet <- metHeatmap(metSet)
#' head(fData(metSet))
setMethod("metHeatmap", signature(x = "ExpressionSet"),
          function(x,
                   dissymC = c("euclidean",
                               "maximum",
                               "manhattan",
                               "canberra",
                               "binary",
                               "minkowski",
                               "1-cor",
                               "1-abs(cor)")[7],
                   correlC = c("pearson",
                               "kendall",
                               "spearman"),
                   aggloC = c("ward.D",
                              "ward.D2",
                              "single",
                              "complete",
                              "average",
                              "mcquitty",
                              "median",
                              "centroid")[2],
                   clustVi = c(2, 2),
                   fig.pdfC = NA,
                   plotCexVn = c(1, 1),
                   plotPaletteC = c("blueOrangeRed",
                                    "redBlackGreen")[1],
                   plotScaleL = TRUE,
                   info.txtC = NA) {
            
            if (!is.null(info.txtC) &&
                !is.na(info.txtC))
              sink(info.txtC, append = TRUE)            
            
            x <- .metHeatmap(x,
                              dissymC,
                              correlC,
                              aggloC,
                              clustVi,
                              fig.pdfC,
                              plotCexVn,
                              plotPaletteC,
                              plotScaleL)
            
            validObject(x)
            
            if (!is.null(info.txtC) &&
                !is.na(info.txtC))
              sink()
            
            return(invisible(x))
            
          })
            
.metHeatmap <- function(eSet,
                         disC,    ## dissimilarity
                         corC, ## correlation method
                         aggC, ## agglomeration method
                         cutVi,
                         fig.pdfC,
                         ploCexVn,
                         ploPalC,    ## color scale
                         ploScaL) {
  
  proMN <- t(Biobase::exprs(eSet))
  
  ncaN <- 14 ## Sample and variable name truncature for display
  
  if (disC %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) {
    
    obsHcl <- hclust(dist(proMN, method = disC),
                     method = aggC)
    
    feaHcl <- hclust(dist(t(proMN), method = disC),
                     method = aggC)
    
  } else if (disC == "1-cor") {
    
    obsHcl <- hclust(as.dist(1 - cor(t(proMN),
                                   method = corC,
                                   use = "pairwise.complete.obs")),
                     method = aggC)
    
    feaHcl <- hclust(as.dist(1 - cor(proMN,
                                   method = corC,
                                   use = "pairwise.complete.obs")),
                     method = aggC)
    
  } else if (disC == "1-abs(cor)") {
    
    obsHcl <- hclust(as.dist(1 - abs(cor(t(proMN),
                                       method = corC,
                                       use = "pairwise.complete.obs"))),
                     method = aggC)
    
    feaHcl <- hclust(as.dist(1 - abs(cor(proMN,
                                       method = corC,
                                       use = "pairwise.complete.obs"))),
                     method = aggC)
    
  }
  
  heaMN <- proMN <- proMN[obsHcl[["order"]], feaHcl[["order"]]]
  
  if (ploScaL)
    heaMN <- scale(heaMN)
  
  heaMN <- heaMN[, rev(1:ncol(heaMN)), drop = FALSE]
  
  switch(ploPalC,
         blueOrangeRed = {
           imaPalVn <- colorRampPalette(c("blue", "orange", "red"),
                                        space = "rgb")(5)[1:5]
         },
         redBlackGreen = {
           imaPalVn <- colorRampPalette(c("red", "black", "green"),
                                        space = "rgb")(5)[1:5]
         })
  
  ## figure
  
  if (!is.na(fig.pdfC))
    pdf(fig.pdfC)
  
  layout(matrix(1:4, nrow = 2),
         widths = c(1, 4), heights = c(1, 4))
  
  ## Color scale
  
  scaN <- length(imaPalVn)
  
  par(mar = c(0.6, 0.6, 0.6, 4.1))
  
  ylimVn <- c(0, scaN)
  ybottomVn <- 0:(scaN - 1)
  ytopVn <- 1:scaN
  
  plot(x = 0,
       y = 0,
       bty = "n",
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
  
  rect(xleft = 0.8,
       ybottom = ybottomVn,
       xright = 1,
       ytop = ytopVn,
       col = imaPalVn,
       border = NA)
  
  prtVn <- pretty(range(heaMN, na.rm = TRUE))
  axis(at = scaN / diff(range(prtVn)) * (prtVn - min(prtVn)),
       font = 2,
       font.axis = 2,
       labels = prtVn,
       las = 1,
       lwd = 2,
       lwd.ticks = 2,
       side = 4,
       xpd = TRUE)
  
  arrows(par("usr")[2],
         par("usr")[4],
         par("usr")[2],
         par("usr")[3],
         code = 0,
         lwd = 2,
         xpd = TRUE)
  
  ## Feature dendrogram
  
  par(mar = c(4.1, 0.6, 0, 0.1),
      lwd = 2)
  
  plot(rev(as.dendrogram(feaHcl)), horiz = TRUE,
       leaflab = "none",
       main = "", xaxs = "i", yaxs = "i",
       xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  
  revFeaHcl <- list(merge = cbind(feaHcl[["merge"]][, 2], feaHcl[["merge"]][, 1]),
                    height = feaHcl[["height"]],
                    order = rev(feaHcl[["order"]]),
                    labels = feaHcl[["labels"]])
  
  if (cutVi[2] > 1) {
    cluFeaVn <- cutree(revFeaHcl, k = cutVi[2])[revFeaHcl[["order"]]]
    cutFeaVn <- which(abs(diff(cluFeaVn)) > 0)
    cutFeaTxtVn <- c(cutFeaVn[1] / 2, cutFeaVn + diff(c(cutFeaVn, length(cluFeaVn))) / 2) + 0.5
    cutFeaLinVn <- cutFeaVn + 0.5
    text(par("usr")[1] + 0.2 * diff(par("usr")[1:2]),
         cutFeaTxtVn,
         labels = unique(cluFeaVn),
         cex = 2,
         font = 2,
         las = 2)
  }
  
  ## Observation dendrogram
  
  par(mar = c(0.1, 0, 0.6, 4.1),
      lwd = 2)
  
  plot(as.dendrogram(obsHcl), leaflab = "none",
       main = "", xaxs = "i", yaxs = "i",
       yaxt = "n", xlab = "", ylab = "")
  
  if (cutVi[1] > 1) {
    cluObsVn <- cutree(obsHcl, k = cutVi[1])[obsHcl[["order"]]]
    cutObsVn <- which(abs(diff(cluObsVn)) > 0)
    cutObsTxtVn <- c(cutObsVn[1] / 2, cutObsVn + diff(c(cutObsVn, length(cluObsVn))) / 2) + 0.5
    cutObsLinVn <- cutObsVn + 0.5
    text(cutObsTxtVn,
         0.8 * par("usr")[4],
         labels =  unique(cluObsVn),
         cex = 2,
         font = 2)
  }
  
  ## Heatmap
  
  par(mar = c(4.1, 0, 0, 4.1))
  
  image(x = 1:nrow(heaMN),
        y = 1:ncol(heaMN),
        z = round(heaMN),
        col = imaPalVn,
        font.axis = 2,
        font.lab = 2,
        xaxt = "n",
        yaxt = "n",
        xlab = "",
        ylab = "")
  
  obsOrdVc <- obsHcl[["labels"]][obsHcl[["order"]]]
  obsOrdLenVn <- sapply(obsOrdVc, nchar)
  obsOrdVc <- substr(obsOrdVc, 1, ncaN)
  obsOrdVc <- paste0(obsOrdVc, ifelse(obsOrdLenVn > ncaN, ".", ""), " ")
  
  mtext(obsOrdVc,
        at = 1:nrow(heaMN),
        cex = ploCexVn[1],
        las = 2,
        side = 1)
  
  feaOrdVc <- feaHcl[["labels"]][feaHcl[["order"]]]
  feaOrdLenVn <- sapply(feaOrdVc, nchar)
  feaOrdVc <- substr(feaOrdVc, 1, ncaN)
  feaOrdVc <- paste0(" ", feaOrdVc, ifelse(feaOrdLenVn > ncaN, ".", ""))
  
  mtext(feaOrdVc,
        at = ncol(heaMN):1,
        cex = ploCexVn[2],
        las = 2,
        side = 4)
  
  if (cutVi[2] > 1)
    abline(h = cutFeaLinVn)
  if (cutVi[1] > 1)
    abline(v = cutObsLinVn)
  
  box()
  
  if (!is.na(fig.pdfC))
    dev.off()
  
  ## Returning
  
  if (cutVi[1] > 1) ## number of sample clusters
    Biobase::pData(eSet)[, "heat_clust"] <- cutree(obsHcl, k = cutVi[1])
#  obsDF <- obsDF[obsHcl[["order"]], , drop = FALSE]
  
  if (cutVi[2] > 1) ## number of variable clusters
    Biobase::fData(eSet)[, "heat_clust"] <- cutree(feaHcl, k = cutVi[2])
#  feaDF <- feaDF[feaHcl[["order"]], , drop = FALSE]
  
  return(invisible(eSet))
  
}
