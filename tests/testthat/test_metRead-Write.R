testthat::context("Testing 'read-write'")

testthat::test_that(".metRead", {
  
  sacDirC <- system.file("extdata/sacurine",
                         package = "metabolis")
  
  sacSet1 <- metabolis:::.metRead(sacDirC)
  
  testthat::expect_true(class(sacSet1) == "ExpressionSet")
  
  testthat::expect_equal(exprs(sacSet1)[1, 1],
                         477491,
                         tolerance = 1e-3)
  # alternatively
  sacSet2 <- metabolis:::.metRead(NA,
                                  file.path(sacDirC, "dataMatrix.tsv"),
                                  file.path(sacDirC, "sampleMetadata.tsv"),
                                  file.path(sacDirC, "variableMetadata.tsv"))
  
  testthat::expect_true(class(sacSet2) == "ExpressionSet")
  
  testthat::expect_equal(Biobase::exprs(sacSet2)[1, 1],
                         477491,
                         tolerance = 1e-3)
  
  
  testthat::expect_error(metabolis:::.metRead(NA,
                                              file.path(sacDirC, "dataMatrix.tsv_XXX"),
                                              file.path(sacDirC, "sampleMetadata.tsv"),
                                              file.path(sacDirC, "variableMetadata.tsv")))
  
})

testthat::test_that("metRead_ExpressionSet", {
  
  sacDirC <- system.file("extdata/sacurine",
                         package = "metabolis")
  
  sacSet2 <- metabolis::metRead(sacDirC)
  
  testthat::expect_true(class(sacSet2) == "ExpressionSet")
  
  testthat::expect_equal(Biobase::exprs(sacSet2)[1, 1],
                         477491,
                         tolerance = 1e-3)
  # alternatively
  sacSet3 <- metabolis::metRead(NA,
                                filesLs = list(dataMatrix.tsvC = file.path(sacDirC, "dataMatrix.tsv"),
                                               sampleMetadata.tsvC = file.path(sacDirC, "sampleMetadata.tsv"),
                                               variableMetadata.tsvC = file.path(sacDirC, "variableMetadata.tsv")))
  
  testthat::expect_true(class(sacSet3) == "ExpressionSet")
  
  testthat::expect_equal(exprs(sacSet3)[1, 1],
                         477491,
                         tolerance = 1e-3)
  
  
  testthat::expect_error(metabolis::metRead(NA,
                                            filesLs = list(dataMatrix.tsvC = file.path(sacDirC, "dataMatrix.tsv_XXX"),
                                                           sampleMetadata.tsvC = file.path(sacDirC, "sampleMetadata.tsv"),
                                                           variableMetadata.tsvC = file.path(sacDirC, "variableMetadata.tsv"))))
  
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
