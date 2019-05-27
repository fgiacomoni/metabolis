testthat::context("Testing 'metTest'")

testthat::test_that("ttest-eset", {
  
  proSet <- metabolis::metRead(system.file("extdata/promet/proteo",
                                          package = "metabolis"))
  
  proSet <- metabolis::metTest(proSet,
                                 factorC = "gene",
                                 testC = "ttest",
                                 fig.pdfC = NULL,
                                 info.txtC = NULL)
  
  testthat::expect_equal(Biobase::fData(proSet)["v3", "gene_ttest_L.W_BH"],
                         0.006009,
                         tolerance = 1e-6)
  
  proSet <- metabolis::metTest(proSet,
                                 factorC = "gene",
                                 testC = "nonparam",
                                 fig.pdfC = NULL,
                                 info.txtC = NULL)
  
  testthat::expect_equal(Biobase::fData(proSet)["v3", "gene_wilcoxon_L.W_BH"],
                         0.008120,
                         tolerance = 1e-6)
  
  proSet <- metabolis::metTest(proSet,
                                 factorC = "gene",
                                 testC = "limma",
                                 fig.pdfC = NULL,
                                 info.txtC = NULL)
  
  testthat::expect_equal(Biobase::fData(proSet)["v3", "gene_limma_L.W_BH"],
                         0.005378,
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
                                      sum(fdaDF[, "gene_ttest_L.W_sig"])
                                    }),
                             c(metabo = 1138, proteo = 348))
  
})
