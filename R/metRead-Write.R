#' metRead
#'
#' Reading datasets as subfolders containing the 3 tables 'dataMatrix.tsv',
#' 'sampleMetadata.tsv' and 'variableMetadata.tsv'
#'
#' @param dirC Character: directory containing the 3 .tsv files (single
#' dataset), or containing several subdirectories with 3 .tsv files
#' (multiple datasets)
#' @param subsetVc Vector of characters: specifying a subset of the
#' subdirectories to be included in the MultiDataSet (by default, all
#' subdirectories containing the 3 tables will be considered as datasets)
#' @param filesLs List: if dirC is set to NA, the full names of the
#' individual files can be provided; in case of an ExpressionSet, the
#' names of the list must be 'dataMatrix.tsvC', 'sampleMetadata.tsvC',
#' and 'variableMetadata.tsvC' with the corresponding file full names;
#' in case of a MultiDataSet, the list must consists of one such sublist
#' per dataset
#' @param info.txtC Character: File name for the printed results (call to
#' 'sink()'); if NA (default), messages will be printed on the screen; if NULL,
#' no verbose will be generated
#' @return MultiDataSet (multiple dataset) instance
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#' prometDirC <- system.file("extdata/promet", package="metabolis")
#' ## 1) Single set
#' metEset <- metRead(file.path(prometDirC, "metabo"))
#' # or
#' metEset <- metRead(NA,
#'                     filesLs = list(dataMatrix.tsvC = file.path(prometDirC, "metabo", "dataMatrix.tsv"),
#'                                    sampleMetadata.tsvC = file.path(prometDirC, "metabo", "sampleMetadata.tsv"),
#'                                    variableMetadata.tsvC = file.path(prometDirC, "metabo", "variableMetadata.tsv")))
#' ## 2) Multiple sets
#' prometMset <- metRead(prometDirC)
#' metMset <- metRead(prometDirC, subsetVc = "metabo")
#' # or
#' prometMset <- metRead(NA,
#'                        filesLs = list(metabo = list(dataMatrix.tsvC = file.path(prometDirC, "metabo", "dataMatrix.tsv"),
#'                                                     sampleMetadata.tsvC = file.path(prometDirC, "metabo", "sampleMetadata.tsv"),
#'                                                     variableMetadata.tsvC = file.path(prometDirC, "metabo", "variableMetadata.tsv")),
#'                                       proteo = list(dataMatrix.tsvC = file.path(prometDirC, "proteo", "dataMatrix.tsv"),
#'                                                     sampleMetadata.tsvC = file.path(prometDirC, "proteo", "sampleMetadata.tsv"),
#'                                                     variableMetadata.tsvC = file.path(prometDirC, "proteo", "variableMetadata.tsv"))))
#' @rdname metRead
#' @export
metRead <- function(dirC,
                     filesLs = NULL,
                     subsetVc = NA,
                     info.txtC = NA) {
  
  if (!is.null(info.txtC) &&
      !is.na(info.txtC))
    sink(info.txtC, append = TRUE)

  x <- NULL

  ## Creating the ExpressionSet or building the list for the MultiDataSet

  if (!is.na(dirC)) {

    if (!file.exists(dirC))
      stop("Directory '", dirC, "' was not found.",
           call. = FALSE)

    if (!file.info(dirC)[, "isdir"])
      stop(dirC, "' is not a directory.",
           call. = FALSE)

    dirVc <- dir(dirC, full.names = TRUE)

    dirVl <- file.info(dirVc)[, "isdir"]

    subDirVc <- dirVc[dirVl]

    if (length(subDirVc) == 0) { ## ExpressionSet

      x <- .metRead(dirC,
                    dataMatrix.tsvC = NA,
                    sampleMetadata.tsvC = NA,
                    variableMetadata.tsvC = NA)

    } else {## MultiDataSet

      names(subDirVc) <- basename(subDirVc)

      subDirVl <- sapply(subDirVc,
                         function(subDirC)
                           file.exists(file.path(subDirC,
                                                 "dataMatrix.tsv")))

      if (sum(subDirVl) == 0) {

        stop("None of the subfolders contains a 'dataMatrix.tsv' file:\n",
             paste(subDirVc, collapse = "\n"),
             "\n", call. = FALSE)

      } else if (sum(!subDirVl) > 0) {

        if (!is.null(info.txtC))
          cat("No 'dataMatrix.tsv' file was found in the following subfolders:\n",
              paste(subDirVc[!subDirVl], collapse = "\n"),
              "\nThe corresponding datasets will be skipped.\n",
              sep = "")

        subDirVc <- subDirVc[subDirVl]

      }

      filesLs <- vector(mode = "list", length = length(subDirVc))
      names(filesLs) <- names(subDirVc)

      for (setC in names(filesLs)) {

        filesLs[setC] <- list(file.path(subDirVc[setC],
                                           c("dataMatrix.tsv",
                                             "sampleMetadata.tsv",
                                             "variableMetadata.tsv")))
        names(filesLs[[setC]]) <- c("dataMatrix.tsvC",
                                       "sampleMetadata.tsvC",
                                       "variableMetadata.tsvC")

      }

    }

  } else if (is.na(dirC)) {

    subNamVc <- names(filesLs)

    if (sum(sapply(filesLs, is.list)) == 0) { ## ExpressionSet

      if (length(subNamVc) == 3 &&
         identical(subNamVc, c("dataMatrix.tsvC",
                               "sampleMetadata.tsvC",
                               "variableMetadata.tsvC"))) {

        x <- .metRead(NA,
                      dataMatrix.tsvC = filesLs[["dataMatrix.tsvC"]],
                      sampleMetadata.tsvC = filesLs[["sampleMetadata.tsvC"]],
                      variableMetadata.tsvC = filesLs[["variableMetadata.tsvC"]])

      } else {

        stop("'filesLs does not contain any sublist nor is a list with names 'dataMatrix.tsvC', 'sampleMetadata.tsvC' and 'variableMetadata.tsvC' giving the corresponding file full names.",
             call. = FALSE)

      }

    } else {## MultiDataSet

      for (setC in names(filesLs)) {

        setLs <- filesLs[[setC]]

        if (!identical(names(setLs),
                      c("dataMatrix.tsvC",
                        "sampleMetadata.tsvC",
                        "variableMetadata.tsvC"))) {

          if (!is.null(info.txtC))
            cat("The names of the following sublist are not 'dataMatrix.tsvC', 'sampleMetadata.tsvC' and 'variableMetadata.tsvC':\n",
                setC,
                "\nThe corresponding datasets will be skipped.\n",
                sep = "")

          filesLs[[setC]] <- NULL

        } else if (!file.exists(setLs[["dataMatrix.tsvC"]])) {

          if (!is.null(info.txtC))
            cat("No 'dataMatrix.tsv' file was found in the following sublist:\n",
                setC,
                "\nThe corresponding datasets will be skipped.\n",
                sep = "")

          filesLs[[setC]] <- NULL

        }

      }

      if (length(filesLs) == 0)
        stop("None of the provided sublists meets the requirements for the creation of a dataset.",
             call. = FALSE)

    }

  }

  ## Creating the MultiDataSet

  if (is.null(x)) {

    if (!any(is.na(subsetVc))) {

      misSetVc <- subsetVc[!(subsetVc %in% names(filesLs))]

      if (length(misSetVc)) {
        stop("The following selected subsets were not found in the subfolder(s) or sublist(s):\n",
             paste(misSetVc, collapse = "\n"),
             call. = FALSE)
      } else
        filesLs <- filesLs[subsetVc]

    }

    # MultiDataSet

    x <- MultiDataSet::createMultiDataSet()

    for (setC in names(filesLs)) {

      setLs <- filesLs[[setC]]

      if (!is.null(info.txtC))
        cat("Reading the '", setC, "' dataset\n",
            sep = "")

      ese <- .metRead(NA,
                      setLs[["dataMatrix.tsvC"]],
                      setLs[["sampleMetadata.tsvC"]],
                      setLs[["variableMetadata.tsvC"]])

      ese$id <- rownames(Biobase::pData(ese))

      x <- MultiDataSet::add_eset(x, ese,
                                  dataset.type = setC,
                                  GRanges = NA)

    }

  }

  validObject(x)

  if (!is.null(info.txtC))
    print(x)
  
  if (!is.null(info.txtC) &&
      !is.na(info.txtC))
    sink()

  return(invisible(x))

}


.metRead <- function(dirC,
                     dataMatrix.tsvC = NA,
                     sampleMetadata.tsvC = NA,
                     variableMetadata.tsvC = NA) {

  if (!is.na(dirC) && !is.na(dataMatrix.tsvC))
    stop("Either 'dirC' or 'dataMatrix.tsvC' argument must be set to NA",
         call. = FALSE)

  if (!is.na(dirC)) {

    tabFilVc <- c(dataMatrix = file.path(dirC, "dataMatrix.tsv"),
                  sampleMetadata = file.path(dirC, "sampleMetadata.tsv"),
                  variableMetadata = file.path(dirC, "variableMetadata.tsv"))

  } else {

    tabFilVc <- c(dataMatrix = dataMatrix.tsvC,
                  sampleMetadata = sampleMetadata.tsvC,
                  variableMetadata = variableMetadata.tsvC)

  }

  for (tabC in names(tabFilVc)) {

    tabFilC <- tabFilVc[tabC]

    if (tabC == "dataMatrix") {

      if (is.na(tabFilC) || !file.exists(tabFilC))
        stop("The provided dataMatrix file was not found:\n",
             tabFilC,
             call. = FALSE)

    } else if (!file.exists(tabFilC)) {

      tabFilC <- NA

      cat("The following '",
          tabC,
          "' file was not found:\n",
          tabFilC,
          "\nThe corresponding '",
          ifelse(tabC == "sampleMetadata", 'phenoData', 'featureData'),
          "' slot will be empty.\n",
          sep = "")

    }

    if (!is.na(tabFilC)) {

      ## R standards for row and column names in matrices and data frames
      .checkRformat(tabFilC)
      
      tabDF <- data.frame(data.table::fread(tabFilC,
                                            header = TRUE,
                                            sep = "\t"),
                          check.names = FALSE,
                          row.names = 1)

      # tabDF <- read.table(tabFilC,
      #                     check.names = FALSE,
      #                     header = TRUE,
      #                     row.names = 1,
      #                     sep = "\t",
      #                     stringsAsFactors = FALSE)

      switch(tabC,
             dataMatrix = {
               tdatMN <- as.matrix(tabDF)
             },
             sampleMetadata = {
               samDF <- tabDF
             },
             variableMetadata = {
               varDF <- tabDF
             })

    } else {

      switch(tabC,
             sampleMetadata = {
               samDF <- data.frame(row.names = colnames(tdatMN))
             },
             variableMetadata = {
               varDF <- data.frame(row.names = rownames(tdatMN))
             })

    }
  }

  chkLs <- .checkW4Mformat(t(tdatMN), samDF, varDF)

  if (!chkLs[["chkL"]]) {
    stop("Sample and/or variable names do not match between your tables. Use the 'info.txtC = NA' argument to get more feedback.",
         call. = FALSE)
  } else if (chkLs[["ordL"]]) {
    tdatMN <- t(chkLs[["datMN"]])
  }

  eset <- Biobase::ExpressionSet(assayData = tdatMN,
                                 phenoData = new("AnnotatedDataFrame",
                                                 data = samDF),
                                 featureData = new("AnnotatedDataFrame",
                                                   data = varDF),
                                 experimentData = new("MIAME",
                                                      title = ifelse(!is.na(dirC),
                                                                     basename(dirC),
                                                                     "")))

  validObject(eset)

  return(eset)

}

# metWrite (mset) ####

#' Exporting a MultiDataSet instance into dataset subfolders,
#' each containing the 3 tabulated files 'dataMatrix.tsv',
#' 'sampleMetadata.tsv', 'variableMetadata.tsv'
#'
#' Note that the \code{dataMatrix} is transposed before
#' export (e.g., the samples are written column wise in the 'dataMatrix.tsv'
#' exported file).
#'
#' @param x An S4 object of class \code{MultiDataSet}
#' @param dirC Character: directory where each dataset subfolder should be created
#' @param fileLs List: alternatively to the dirC argument, the full names of the files can be provided as a list
#' @param overwriteL Logical: should existing files be overwritten?
#' @param info.txtC Character: File name for the printed results (call to
#' 'sink()'); if NA (default), messages will be printed on the screen; if NULL,
#' no verbose will be generated
#' @return No object returned.
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#' prometMset <- metRead(system.file("extdata/promet",package="metabolis"))
#'\dontrun{
#' metWrite(prometMset, dirC = file.path(getwd(), "promet"))
#' # alternatively
#' metWrite(prometMset,
#'          dirC = NA,
#'          filesLs = list(metabo = list(dataMatrix.tsvC = file.path(getwd(), "met_dataMatrix.tsv"),
#'                                       sampleMetadata.tsvC = file.path(getwd(), "met_sampleMetadata.tsv"),
#'                                       variableMetadata.tsvC = file.path(getwd(), "met_variableMetadata.tsv")),
#'                         proteo = list(dataMatrix.tsvC = file.path(getwd(), "pro_dataMatrix.tsv"),
#'                                       sampleMetadata.tsvC = file.path(getwd(), "pro_sampleMetadata.tsv"),
#'                                       variableMetadata.tsvC = file.path(getwd(), "pro_variableMetadata.tsv"))))
#'}
#'
#' @rdname metWrite
#' @export
setMethod("metWrite", "MultiDataSet",
          function(x,
                   dirC,
                   filesLs = NULL,
                   overwriteL = FALSE,
                   info.txtC = NA) {
            
            if (!is.null(info.txtC) &&
                !is.na(info.txtC))
              sink(info.txtC, append = TRUE)
            
            infTxtC <- info.txtC
            if (!is.null(infTxtC))
              infTxtC <- NA

            setVc <- names(x)

            if (!is.na(dirC)) {

              setDirVc <- file.path(dirC, setVc)
              names(setDirVc) <- setVc

              if (file.exists(dirC) && file.info(dirC)[, "isdir"]) {

                dirVc <- dir(dirC, full.names = TRUE)

                dirVl <- file.info(dirVc)[, "isdir"]

                subDirVc <- dirVc[dirVl]

                subDupVc <- intersect(setDirVc,
                                      subDirVc)

                if (length(subDupVc) && !overwriteL)
                  stop("The following subfolder(s) were detected in your directory, please remove them or specify another parent directory to avoid overwriting:\n",
                       paste(subDupVc, collapse = "\n"),
                       call. = FALSE)

              } else {

                dir.create(dirC,
                           showWarnings = !is.null(info.txtC))

              }

              for (setC in names(setDirVc)) {

                if (!is.null(info.txtC))
                  cat("Writing the '", setC, "' dataset\n",
                      sep = "")

                dir.create(setDirVc[setC],
                           showWarnings = !is.null(info.txtC))
                


                metabolis::metWrite(x[[setC]],
                                    setDirVc[setC],
                                    overwriteL = overwriteL,
                                    info.txtC = infTxtC)

              }

              if (!is.null(info.txtC))
                cat("The subfolders have been written in the directory:\n",
                    dirC,
                    "\n")

            } else if (is.na(dirC)) {

              if (is.null(filesLs))
                stop("'filesLs' must be provided when 'dirC' is set to NA",
                     call. = FALSE)

              if (is.null(names(filesLs)) || any(is.na(names(filesLs))))
                stop("All names of the sublists must be provided (they should match the names of the MultiDataSet datasets)",
                     call. = FALSE)

              filLisVl <- sapply(filesLs, is.list)

              if (!all(filLisVl))
                stop("The following element(s) of 'filesLs' is/are not sublist(s):\n",
                     paste(names(filLisVl)[!filLisVl], collapse = "\n"),
                     call. = FALSE)

              if (!identical(setVc, names(filLisVl)))
                stop("The name(s) of the 'x' MultiDataSet:\n",
                     paste(setVc, collapse = ", "),
                     "\ndo(es) not match the names of the sublists:\n",
                     paste(names(filLisVl), collapse = ", "),
                     call. = FALSE)

              for (setC in setVc) {

                filLs <- filesLs[[setC]]

                if (length(filLs) != 3 ||
                   !identical(names(filLs), c("dataMatrix.tsvC",
                                              "sampleMetadata.tsvC",
                                              "variableMetadata.tsvC")))
                  stop("The '",
                       setC,
                       "' sublist of 'filesLs does not have the names 'dataMatrix.tsvC', 'sampleMetadata.tsvC' and 'variableMetadata.tsvC'",
                       call. = FALSE)

                if (is.na(filLs[["dataMatrix.tsvC"]]))
                  stop("The 'dataMatrix.tsvC' file name from the '", setC,
                       "' sublist is missing (ie set to NA).",
                       call. = FALSE)

                for (filC in names(filLs)) {

                  filFulNamC <- filLs[[filC]]

                  if (!is.na(filFulNamC)
                      && file.exists(filFulNamC)
                      && !overwriteL)
                    stop("The following file from the '",
                         setC, "' sublist already exists:\n",
                         filFulNamC,
                         "\nPlease remove it or choose another name to avoid overwriting.",
                         call. = FALSE)

                }

              }

              for (setC in setVc) {

                filLs <- filesLs[[setC]]

                if (!is.null(info.txtC))
                  cat("Writing the '", setC, "' dataset\n",
                      sep = "")

                metabolis::metWrite(x[[setC]],
                                    NA,
                                    filLs[["dataMatrix.tsvC"]],
                                    filLs[["sampleMetadata.tsvC"]],
                                    filLs[["variableMetadata.tsvC"]],
                                    overwriteL = overwriteL,
                                    info.txtC = infTxtC)

              }
            }
            
            if (!is.null(info.txtC) &&
                !is.na(info.txtC))
              sink()

          })


# metWrite (eset) ####

#' Exporting an ExpressionSet instance into dataset subfolders,
#' each containing the 3 tabulated files 'dataMatrix.tsv',
#' 'sampleMetadata.tsv', 'variableMetadata.tsv'
#'
#' Note that the \code{dataMatrix} is transposed before
#' export (e.g., the samples are written column wise in the 'dataMatrix.tsv'
#' exported file).
#'
#' @param x An S4 object of class \code{ExpressionSet}
#' @param dirC Character: directory where each dataset should be written
#' @param dataMatrix.tsvC Character: alternatively, the FULL name of the dataMatrix, sampleMetadata, and variableMetadata files can be specified here and in the following two arguments (in that case, the dirC argument must be set to NA); note that if the sampleMetadata and/or variableMetadata file names are left to NA, the corresponding phenoData and featureData slots will be empty in the output ExpressionSet.
#' @param sampleMetadata.tsvC Character: full name of the sampleMetadata file (when dirC is set to NA)
#' @param variableMetadata.tsvC Character: full name of the variableMetadata file (when dirC is set to NA)
#' @param overwriteL Logical: should existing files be overwritten?
#' @param info.txtC Character: File name for the printed results (call to
#' 'sink()'); if NA (default), messages will be printed on the screen; if NULL,
#' no verbose will be generated
#' @return No object returned.
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#' metSet <- metabolis::metRead(system.file("extdata/promet/metabo", package="metabolis"))
#'\dontrun{
#' metWrite(metSet, dirC = file.path(getwd(), "metabo"))
#'}
#' @rdname metWrite
setMethod("metWrite", "ExpressionSet",
          function(x,
                   dirC,
                   dataMatrix.tsvC = NA,
                   sampleMetadata.tsvC = NA,
                   variableMetadata.tsvC = NA,
                   overwriteL = FALSE,
                   info.txtC = NA){

            if (!is.na(dirC) && !is.na(dataMatrix.tsvC))
              stop("Either 'dirC' or 'dataMatrix.tsvC' argument must be set to NA",
                   call. = FALSE)
            
            if (!is.null(info.txtC) &&
                !is.na(info.txtC))
              sink(info.txtC, append = TRUE)

            if (!is.na(dirC)) {

              if (!(file.exists(dirC) && file.info(dirC)[, "isdir"]))
                dir.create(dirC,
                           showWarnings = !is.null(info.txtC))

              tabFilVc <- c(dataMatrix = file.path(dirC, "dataMatrix.tsv"),
                            sampleMetadata = file.path(dirC, "sampleMetadata.tsv"),
                            variableMetadata = file.path(dirC, "variableMetadata.tsv"))

              for (tabC in names(tabFilVc)) {

                tabFilC <- tabFilVc[tabC]

                if (file.exists(tabFilC) && !overwriteL)
                  stop("The following file already exists:\n", tabFilC,
                       "\nPlease remove all 'dataMatrix.tsv', 'sampleMetadata.tsv' and 'variableMetadata.tsv' files, or specify another directory to avoid overwriting.",
                       call. = FALSE)

              }

            } else if (is.na(dirC)) {

              tabFilVc <- c(dataMatrix = dataMatrix.tsvC,
                            sampleMetadata = sampleMetadata.tsvC,
                            variableMetadata = variableMetadata.tsvC)

              for (tabC in names(tabFilVc)) {

                tabFilC <- tabFilVc[tabC]

                if (!is.na(tabFilC) && file.exists(tabFilC) && !overwriteL)
                  stop("The following file already exists:\n", tabFilC,
                       "\nPlease specify another file name.",
                       call. = FALSE)

              }

            }

            ## Writing

            tdatMN <- Biobase::exprs(x)
            samDF <- Biobase::pData(x)
            varDF <- Biobase::fData(x)
            chkLs <- .checkW4Mformat(t(tdatMN), samDF, varDF)

            if (!chkLs[["chkL"]]) {
              stop("Sample and/or variable names do not match between your tables. Use the 'info.txtC = NA' argument to get more feedback.",
                   call. = FALSE)
            } else if (chkLs[["ordL"]]) {
              tdatMN <- t(chkLs[["datMN"]])
            }

            datDF <- cbind.data.frame(dataMatrix = rownames(tdatMN),
                                      as.data.frame(tdatMN))

            write.table(datDF,
                        file = tabFilVc['dataMatrix'],
                        quote = FALSE,
                        row.names = FALSE,
                        sep = "\t")

            if (!is.na(dirC) || !is.na(sampleMetadata.tsvC)) {

              samDF <- cbind.data.frame(sampleMetadata = rownames(samDF),
                                        samDF)
              write.table(samDF,
                          file = tabFilVc['sampleMetadata'],
                          quote = FALSE,
                          row.names = FALSE,
                          sep = "\t")

            }

            if (!is.na(dirC) || !is.na(variableMetadata.tsvC)) {

              varDF <- cbind.data.frame(variableMetadata = rownames(varDF),
                                        varDF)
              write.table(varDF,
                          file = tabFilVc['variableMetadata'],
                          quote = FALSE,
                          row.names = FALSE,
                          sep = "\t")

            }

            if (!is.null(info.txtC))
              cat("The following file(s) have been written:\n",
                  paste(tabFilVc[!is.na(basename(tabFilVc))], collapse = "\n"),
                  "\n")
            
            if (!is.null(info.txtC) &&
                !is.na(info.txtC))
              sink()

          })


.checkRformat <- function(filCa) {

  rowVc <- data.table::fread(filCa,
                             header = TRUE,
                             sep = "\t")[[1]]
  
  # rowVc <- read.table(filCa,
  #                     check.names = FALSE,
  #                     header = TRUE,
  #                     sep = "\t",
  #                     stringsAsFactors = FALSE)[, 1]

  colVc <- unlist(data.table::fread(filCa,
                                    header = FALSE,
                                    nrows = 1,
                                    sep = "\t")[1])[-1]
  
  # colVc <- unlist(read.table(filCa,
  #                            check.names = FALSE,
  #                            nrows = 1,
  #                            sep = "\t",
  #                            stringsAsFactors = FALSE))[-1]

  if (any(duplicated(rowVc)))
    stop("The following ",
         ifelse(names(filCa) == 'sampleMetadata', 'sample', 'variable'),
         " name(s) is/are duplicated in the ",
         names(filCa),
         ": '",
         paste(rowVc[duplicated(rowVc)], collapse = "', '"), "'",
         call. = FALSE)

  if (any(duplicated(colVc)))
    stop("The following ", ifelse(names(filCa) == 'sampleMetadata', 'variable', 'sample'), " name(s) is/are duplicated in the ",
         names(filCa),
         ": '",
         paste(colVc[duplicated(colVc)], collapse = "', '"), "'",
         call. = FALSE)

}


.checkW4Mformat <- function(datMNw, samDFw, varDFw) {

  chkL <- TRUE
  ordL <- FALSE

  if (mode(datMNw) != "numeric") {
    cat("The dataMatrix is not of the 'numeric' type\n")
    chkL <- FALSE
  }

  if (!identical(rownames(datMNw), rownames(samDFw))) {
    ## checking sample names

    datSamDifVc <- setdiff(rownames(datMNw), rownames(samDFw))

    if (length(datSamDifVc)) {
      cat("The following samples were found in the dataMatrix column names but not in the sampleMetadata row names:\n", sep = "")
      print(cbind.data.frame(col = as.numeric(sapply(datSamDifVc,
                                                     function(samC) which(rownames(datMNw) == samC))),
                             name = datSamDifVc))
      chkL <- FALSE
    }

    samDatDifVc <- setdiff(rownames(samDFw), rownames(datMNw))

    if (length(samDatDifVc)) {
      cat("The following samples were found in the sampleMetadata row names but not in the dataMatrix column names:\n",
          sep = "")
      print(cbind.data.frame(row = as.numeric(sapply(samDatDifVc, function(samC) which(rownames(samDFw) == samC))),
                             name = samDatDifVc))
      chkL <- FALSE
    }

    if (nrow(datMNw) != nrow(samDFw)) {
      cat("The dataMatrix has ", nrow(datMNw), " columns (ie samples) whereas the sampleMetadata has ", nrow(samDFw), " rows\n",
          sep = "")
      chkL <- FALSE
    } else if (identical(gsub("^X", "", rownames(datMNw)), rownames(samDFw))) {
      cat("The dataMatrix column names start with an 'X' but not the sampleMetadata row names\n", sep = "")
      chkL <- FALSE
    } else if (identical(gsub("^X", "", rownames(samDFw)), rownames(datMNw))) {
      cat("The sampleMetadata row names start with an 'X' but not the dataMatrix column names\n", sep = "")
      chkL <- FALSE
    } else if (identical(sort(rownames(datMNw)), sort(rownames(samDFw)))) {
      cat("Message: Re-ordering dataMatrix sample names to match sampleMetadata\n")
      datMNw <- datMNw[rownames(samDFw), , drop = FALSE]
      stopifnot(identical(sort(rownames(datMNw)), sort(rownames(samDFw))))
      ordL <- TRUE
    } else {
      cat("The dataMatrix column names and the sampleMetadata row names are not identical:\n", sep = "")
      print(cbind.data.frame(indice = 1:nrow(datMNw),
                             dataMatrix_columnnames = rownames(datMNw),
                             sampleMetadata_rownames = rownames(samDFw))[rownames(datMNw) != rownames(samDFw), , drop = FALSE])
      chkL <- FALSE
    }

  }

  if (!identical(colnames(datMNw), rownames(varDFw))) {
    ## checking variable names

    datVarDifVc <- setdiff(colnames(datMNw), rownames(varDFw))

    if (length(datVarDifVc)) {
      cat("The following variables were found in the dataMatrix row names but not in the variableMetadata row names:\n", sep = "")
      print(cbind.data.frame(row = as.numeric(sapply(datVarDifVc, function(varC) which(colnames(datMNw) == varC))),
                             name = datVarDifVc))
      chkL <- FALSE
    }

    varDatDifVc <- setdiff(rownames(varDFw), colnames(datMNw))

    if (length(varDatDifVc)) {
      cat("The following variables were found in the variableMetadata row names but not in the dataMatrix row names:\n", sep = "")
      print(cbind.data.frame(row = as.numeric(sapply(varDatDifVc, function(varC) which(rownames(varDFw) == varC))),
                             name = varDatDifVc))
      chkL <- FALSE
    }

    if (ncol(datMNw) != nrow(varDFw)) {
      cat("The dataMatrix has ",
          nrow(datMNw),
          " rows (ie variables) whereas the variableMetadata has ",
          nrow(varDFw),
          " rows\n",
          sep = "")
      chkL <- FALSE
    } else if (identical(sort(colnames(datMNw)), sort(rownames(varDFw)))) {
      cat("Message: Re-ordering dataMatrix variable names to match variableMetadata\n")
      datMNw <- datMNw[, rownames(varDFw), drop = FALSE]
      stopifnot(identical(sort(colnames(datMNw)), sort(rownames(varDFw))))
      ordL <- TRUE
    } else {
      cat("\n\nThe dataMatrix row names and the variableMetadata row names are not identical:\n",
          sep = "")
      print(cbind.data.frame(row = 1:ncol(datMNw),
                             dataMatrix_rownames = colnames(datMNw),
                             variableMetadata_rownames = rownames(varDFw))[colnames(datMNw) != rownames(varDFw), , drop = FALSE])
      chkL <- FALSE
    }
  }

  return(list(chkL = chkL,
              ordL = ordL,
              datMN = datMNw))

}















