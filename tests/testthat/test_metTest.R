testthat::context("Testing 'metTest'")

testthat::test_that("ttest-eset", {
  
  sacSet <- metabolis::metRead(system.file("extdata/sacurine",
                                           package = "metabolis"))
  
  sacSet <- metabolis::metCorrect(sacSet)
  sacSet <- sacSet[, Biobase::pData(sacSet)[, "sampleType"] != "pool"]
  sacSet <- metabolis::metTransform(sacSet)
  sacSet <- sacSet[, Biobase::sampleNames(sacSet) != "HU_neg_096_b2"]
  
  sacSet <- metabolis::metTest(sacSet,
                               factorC = "gender",
                               testC = "ttest",
                               fig.pdfC = NULL,
                               info.txtC = NULL)
  
  testthat::expect_equal(Biobase::fData(sacSet)["1,7-Dimethyluric acid",
                                                "gender_ttest_Female.Male_BH"],
                         0.5868704,
                         tolerance = 1e-6)
  
  sacSet <- metabolis::metTest(sacSet,
                               factorC = "gender",
                               testC = "nonparam",
                               fig.pdfC = NULL,
                               info.txtC = NULL)
  
  testthat::expect_equal(Biobase::fData(sacSet)["1-Methylxanthine",
                                                "gender_wilcoxon_Female.Male_BH"],
                         0.03904366,
                         tolerance = 1e-6)
  
  sacSet <- metabolis::metTest(sacSet,
                               factorC = "gender",
                               testC = "limma",
                               fig.pdfC = NULL,
                               info.txtC = NULL)
  
  testthat::expect_equal(Biobase::fData(sacSet)["1-Methylxanthine",
                                                "gender_limma_Female.Male_BH"],
                         0.0701083,
                         tolerance = 1e-6)
  
})

testthat::test_that("ttest-mset", {
  
  prometMset <- metabolis::metRead(system.file("extdata/promet",
                                              package = "metabolis"))
  
  prometMset <- metabolis::metTest(prometMset,
                                     factorC = "gene",
                                     testC = "ttest",
                                     fig.pdfC = NULL,
                                     info.txtC = NULL)
  
  testthat::expect_identical(sapply(fData(prometMset),
                                    function(fdaDF) {
                                      sum(fdaDF[, "gene_ttest_KO.WT_sig"])
                                    }),
                             c(metabo = 2, proteo = 7))
  
})
