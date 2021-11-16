test_that("findChain works", {

    expect_output( findBeta( pathFileInput = "C:/Users/Shashank/Desktop/peptides_second rep.csv",
                             pathDirOutput = "C:/Users/Shashank/Downloads/", TRUE ),
                   "Analysis completed successfully!" )

})
