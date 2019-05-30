testthat::context("Testing 'metTest'")

testthat::test_that("ttest-eset", {
  
  proSet <- metabolis::metRead(system.file("extdata/promet/proteo",
                                          package = "metabolis"))
  
  proSet <- metabolis::metTest(proSet,
                                 factorC = "gene",
                                 testC = "ttest",
                                 fig.pdfC = NULL,
                                 info.txtC = NULL)
  
  testthat::expect_equal(Biobase::fData(proSet)["v3", "gene_ttest_KO.WT_BH"],
                         0.006818127,
                         tolerance = 1e-6)
  
  proSet <- metabolis::metTest(proSet,
                                 factorC = "gene",
                                 testC = "nonparam",
                                 fig.pdfC = NULL,
                                 info.txtC = NULL)
  
  testthat::expect_equal(Biobase::fData(proSet)["v3", "gene_wilcoxon_KO.WT_BH"],
                         0.007156371,
                         tolerance = 1e-6)
  
  proSet <- metabolis::metTest(proSet,
                                 factorC = "gene",
                                 testC = "limma",
                                 fig.pdfC = NULL,
                                 info.txtC = NULL)
  
  testthat::expect_equal(Biobase::fData(proSet)["v3", "gene_limma_KO.WT_BH"],
                         0.004541396,
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
