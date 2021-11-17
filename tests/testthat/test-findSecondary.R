test_that("findSecondary works", {

    pathFileInput = system.file( "extdata", "exampleInputFile.csv", package = "ypssc" )

    temp = tempdir()
    if ( !dir.exists(temp) ) {
        dir.create(temp)
    }
    pathDirOutput = temp

    expect_output( findSecondary( pathFileInput = "C:/Users/Shashank/Desktop/peptides_second rep.csv",
                                  pathDirOutput = "C:/Users/Shashank/Downloads/", TRUE ),
                   "Analysis completed successfully!" )

    unlink( pathDirOutput, recursive = TRUE )

})
