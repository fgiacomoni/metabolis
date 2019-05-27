testthat::context("Testing 'metView'")

testthat::test_that("metView-eset", {
  
  proSet <- metabolis::readSet(system.file("extdata/promet/proteo",
                                           package = "metabolis"))
  
  metabolis::metView(proSet,
                     fig.pdfC = "figures/test-metView-eset.pdf",
                     info.txtC = NULL)
  metabolis::metView(proSet,
                     factorC = "gene",
                     fig.pdfC = "figures/test-metView-eset_gene.pdf",
                     info.txtC = NULL)
  
})

testthat::test_that("metView-mset", {
  
  prometMset <- metabolis::readSet(system.file("extdata/promet",
                                               package = "metabolis"))
  
  metabolis::metView(prometMset,
                     fig.pdfC = "figures/test-metView-mset.pdf",
                     info.txtC = NULL)
  metabolis::metView(prometMset,
                     factorC = "gene",
                     plotOverviewL = FALSE,
                     fig.pdfC = "figures/test-metView-mset_gene.pdf",
                     info.txtC = NULL)
  
})
