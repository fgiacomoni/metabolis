testthat::context("Testing 'read-write'")

testthat::test_that(".metRead", {
  
  metDirC <- system.file("extdata/promet/metabo",
                         package = "metabolis")
  
  metSet1 <- metabolis:::.metRead(metDirC)
  
  testthat::expect_true(class(metSet1) == "ExpressionSet")
  
  testthat::expect_equal(exprs(metSet1)[1, 1],
                         5.576,
                         tolerance = 1e-3)
  # alternatively
  metSet2 <- metabolis:::.metRead(NA,
                                  file.path(metDirC, "dataMatrix.tsv"),
                                  file.path(metDirC, "sampleMetadata.tsv"),
                                  file.path(metDirC, "variableMetadata.tsv"))
  
  testthat::expect_true(class(metSet2) == "ExpressionSet")
  
  testthat::expect_equal(Biobase::exprs(metSet2)[1, 1],
                         5.576,
                         tolerance = 1e-3)
  
  
  testthat::expect_error(metabolis:::.metRead(NA,
                                              file.path(metDirC, "dataMatrix.tsv_XXX"),
                                              file.path(metDirC, "sampleMetadata.tsv"),
                                              file.path(metDirC, "variableMetadata.tsv")))
  
})

testthat::test_that("metRead_ExpressionSet", {
  
  metDirC <- system.file("extdata/promet/metabo",
                         package = "metabolis")
  
  metSet2 <- metabolis::metRead(metDirC)
  
  testthat::expect_true(class(metSet2) == "ExpressionSet")
  
  testthat::expect_equal(Biobase::exprs(metSet2)[1, 1],
                         5.576,
                         tolerance = 1e-3)
  # alternatively
  metSet3 <- metabolis::metRead(NA,
                                filesLs = list(dataMatrix.tsvC = file.path(metDirC, "dataMatrix.tsv"),
                                               sampleMetadata.tsvC = file.path(metDirC, "sampleMetadata.tsv"),
                                               variableMetadata.tsvC = file.path(metDirC, "variableMetadata.tsv")))
  
  testthat::expect_true(class(metSet3) == "ExpressionSet")
  
  testthat::expect_equal(exprs(metSet3)[1, 1],
                         5.576,
                         tolerance = 1e-3)
  
  
  testthat::expect_error(metabolis::metRead(NA,
                                            filesLs = list(dataMatrix.tsvC = file.path(metDirC, "dataMatrix.tsv_XXX"),
                                                           sampleMetadata.tsvC = file.path(metDirC, "sampleMetadata.tsv"),
                                                           variableMetadata.tsvC = file.path(metDirC, "variableMetadata.tsv"))))
  
})


testthat::test_that("metRead_MultiDataSet", {
  
  prometDirC <- system.file("extdata/promet",
                            package = "metabolis")
  
  prometMset1 <- metabolis::metRead(prometDirC)
  
  testthat::expect_true(class(prometMset1) == "MultiDataSet")
  
  testthat::expect_identical(names(prometMset1), c("metabo", "proteo"))
  
  testthat::expect_equal(exprs(prometMset1[["metabo"]])[1, 1],
                         5.576,
                         tolerance = 1e-3)
  
  
  prometMset2 <- metabolis::metRead(NA,
                                    filesLs = list(metabo = list(dataMatrix.tsvC = file.path(prometDirC, "metabo", "dataMatrix.tsv"),
                                                                 sampleMetadata.tsvC = file.path(prometDirC, "metabo", "sampleMetadata.tsv"),
                                                                 variableMetadata.tsvC = file.path(prometDirC, "metabo", "variableMetadata.tsv")),
                                                   proteo = list(dataMatrix.tsvC = file.path(prometDirC, "proteo", "dataMatrix.tsv"),
                                                                 sampleMetadata.tsvC = file.path(prometDirC, "proteo", "sampleMetadata.tsv"),
                                                                 variableMetadata.tsvC = file.path(prometDirC, "proteo", "variableMetadata.tsv"))))
  
  testthat::expect_true(class(prometMset2) == "MultiDataSet")
  
  testthat::expect_identical(names(prometMset2), c("metabo", "proteo"))
  
  testthat::expect_equal(exprs(prometMset2[["metabo"]])[1, 1],
                         5.576,
                         tolerance = 1e-3)
  
  testthat::expect_identical(colnames(pData(prometMset2[["metabo"]])), c("gene", "id"))
  
  
  
  metMset1 <- metabolis::metRead(prometDirC, subsetVc = "metabo")
  
  testthat::expect_true(class(metMset1) == "MultiDataSet")
  
  testthat::expect_identical(names(metMset1), "metabo")
  
  testthat::expect_equal(exprs(metMset1[["metabo"]])[1, 1],
                         5.576,
                         tolerance = 1e-3)
  
  
  metMset2 <- metabolis::metRead(NA,
                                 filesLs = list(metabo = list(dataMatrix.tsvC = file.path(prometDirC, "metabo", "dataMatrix.tsv"),
                                                              sampleMetadata.tsvC = file.path(prometDirC, "metabo", "sampleMetadata.tsv"),
                                                              variableMetadata.tsvC = file.path(prometDirC, "metabo", "variableMetadata.tsv")),
                                                proteo = list(dataMatrix.tsvC = file.path(prometDirC, "proteo", "dataMatrix.tsv"),
                                                              sampleMetadata.tsvC = file.path(prometDirC, "proteo", "sampleMetadata.tsv"),
                                                              variableMetadata.tsvC = file.path(prometDirC, "proteo", "variableMetadata.tsv"))),
                                 subsetVc = "metabo")
  
  testthat::expect_true(class(metMset2) == "MultiDataSet")
  
  testthat::expect_identical(names(metMset2), "metabo")
  
  testthat::expect_equal(exprs(metMset2[["metabo"]])[1, 1],
                         5.576,
                         tolerance = 1e-3)
  
  
  testthat::expect_error(metabolis::metRead(NA,
                                            list(metabo = list(dataMatrix.tsvC = file.path(prometDirC, "metabo", "dataMatrix.tsv_XXX"),
                                                               sampleMetadata.tsvC = file.path(prometDirC, "metabo", "sampleMetadata.tsv"),
                                                               variableMetadata.tsvC = file.path(prometDirC, "metabo", "variableMetadata.tsv")),
                                                 proteo = list(dataMatrix.tsvC = file.path(prometDirC, "proteo", "dataMatrix.tsv_XXX"),
                                                               sampleMetadata.tsvC = file.path(prometDirC, "proteo", "sampleMetadata.tsv"),
                                                               variableMetadata.tsvC = file.path(prometDirC, "proteo", "variableMetadata.tsv")))))
  
  prometMset3 <- metabolis::metRead(NA,
                                    list(metabo = list(dataMatrix.tsvC = file.path(prometDirC, "metabo", "dataMatrix.tsv"),
                                                       sampleMetadata.tsvC = file.path(prometDirC, "metabo", "sampleMetadata.tsv_XXX"),
                                                       variableMetadata.tsvC = file.path(prometDirC, "metabo", "variableMetadata.tsv")),
                                         proteo = list(dataMatrix.tsvC = file.path(prometDirC, "proteo", "dataMatrix.tsv"),
                                                       sampleMetadata.tsvC = file.path(prometDirC, "proteo", "sampleMetadata.tsv"),
                                                       variableMetadata.tsvC = file.path(prometDirC, "proteo", "variableMetadata.tsv"))))
  
  testthat::expect_true(colnames(pData(prometMset3[["metabo"]])) == "id")
  
  
})

testthat::test_that("metWrite_MultiDataSet", {
  
  prometDirC <- system.file("extdata/promet",
                            package = "metabolis")
  
  prometMset4 <- metRead(prometDirC)
  
  testthat::expect_error(metabolis::metWrite(prometMset4,
                                             dirC = NA,
                                             filesLs = list(metabo = list(dataMatrix.tsvC = NA,
                                                                          sampleMetadata.tsvC = file.path(getwd(), "metabo_sampleMetadata.tsv"),
                                                                          variableMetadata.tsvC = file.path(getwd(), "metabo_variableMetadata.tsv")),
                                                            proteo = list(dataMatrix.tsvC = file.path(getwd(), "proteo_dataMatrix.tsv"),
                                                                          sampleMetadata.tsvC = file.path(getwd(), "proteo_sampleMetadata.tsv"),
                                                                          variableMetadata.tsvC = file.path(getwd(), "proteo_variableMetadata.tsv")))))
  
})
