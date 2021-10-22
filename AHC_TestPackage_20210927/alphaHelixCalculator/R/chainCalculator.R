
####################################################################################################################################
####################################################################################################################################
# >>
#' @title chainCalculator xxx xxx xxx
#' @description xxxx xxxx This function does bla bla bla xxxx xxxx.
#'   xxxxx xxxxxx xxxxx xxxx
#'   xxxxx xxxxxx xxxxx xxxx
#'   xxxxx xxxxxx xxxxx xxxx
#'   xxxxx xxxxxx xxxxx xxxx
#' @param pathFileInput input file path from which bla bla bla xxxx xxxx xxxx xxxx xxxx xxxx xxxx.
#' @param pathDirOutput directory path to which the output files will be generated xxxx xxxx xxxx.
#' @return xxxx xxxx xxxx xxxx xxxx xxxx xxxx xxxx
#' @details DETAILS
#' @examples
#' \dontrun{
#' chainCalculator( pathFileInput = "<someInputFileName>",
#'                  pathDirOutput = "<someOutputFolderName>")
#' }
#' @export
# <<
####################################################################################################################################
####################################################################################################################################


####################################################################################################################################
##################################### chainCalculator() ############################################################################
# >>
chainCalculator <- function( pathFileInput = "C:/Users/Shashank/Desktop/peptides_second rep.csv",
                             pathDirOutput = "C:/Users/Shashank/Downloads" ) {

    # Begin >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    print("Started")
    startTime = Sys.time()

    # Getting current working directory >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    originalWorkingDir = getwd()

    # Checking if 'pathDirOutput' is provided >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    if ( is.null( pathDirOutput ) ) {
        pathDirOutput = getwd()
    }

    # Reading the input sample file >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    dataFileInput     = readFileInput( pathFileInput )
    df                = dataFileInput$df
    sampleNames       = dataFileInput$sampleNames
    sampleNamesUpdate = dataFileInput$sampleNamesUpdate

    # Create output folder >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    dateTimeCurrent = creatOutputDir( pathDirOutput )

    # Removing the rows that are not needed >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    df = removeRows( df, dateTimeCurrent )

    # Chain calculation for dataBase >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    chainCalculation     ( df, sampleNames, sampleNamesUpdate )

    # Setting working directory back to original >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    setwd( originalWorkingDir )

    # End >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    endTime   = Sys.time()
    timeTaken = endTime - startTime
    print( paste0( "Time taken for the AHC run: ", format(timeTaken) ) )

    return( invisible(NULL) )
}
# <<
##################################### chainCalculator() ############################################################################
####################################################################################################################################
