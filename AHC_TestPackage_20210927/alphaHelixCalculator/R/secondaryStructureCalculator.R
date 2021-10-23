
####################################################################################################################################
####################################################################################################################################
# >>
#' @title secondaryStructureCalculator xxx xxx xxx
#' @description xxxx xxxx This function does bla bla bla xxxx xxxx.
#'   xxxxx xxxxxx xxxxx xxxx
#'   xxxxx xxxxxx xxxxx xxxx
#'   xxxxx xxxxxx xxxxx xxxx
#'   xxxxx xxxxxx xxxxx xxxx
#' @param pathFileInput input file path from which bla bla bla xxxx xxxx xxxx xxxx xxxx xxxx xxxx.
#' @param pathDirOutput directory path to which the output files will be generated xxxx xxxx xxxx.
#' @import dplyr
#' @import readxl
#' @import stringr
#' @import eulerr
#' @import ggplot2
#' @import Peptides
#' @import utils
#' @import svDialogs
#' @import tcltk
#' @return NULL
#' @details DETAILS
#' @examples
#' \dontrun{
#' secondaryStructureCalculator( pathFileInput = "<someInputFileName>",
#'                               pathDirOutput = "<someOutputFolderName>")
#' }
#' @export
# <<
####################################################################################################################################
####################################################################################################################################


####################################################################################################################################
##################################### secondaryStructureCalculator() ###############################################################
# >>
secondaryStructureCalculator <- function( pathFileInput = "C:/Users/Shashank/Desktop/peptides_second rep_small.csv",
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

    # Writing `dataBase_numOfAA` >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    write.csv( dataBase_numOfAA,
               "dataBase_numOfAA.csv",
               row.names = FALSE )

    # Alpha helix calculation for dataBase >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    alphaHelixCalculation( df, sampleNames, sampleNamesUpdate, dateTimeCurrent )

    # Beta-sheet calculation for dataBase >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    betaSheetCalculation ( df, sampleNames, sampleNamesUpdate, dateTimeCurrent )

    # Chain calculation for dataBase >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    chainCalculation     ( df, sampleNames, sampleNamesUpdate, dateTimeCurrent )

    # End >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    endTime   = Sys.time()
    timeTaken = endTime - startTime
    print( paste0( "Time taken for the AHC run: ", format(timeTaken) ) )

    # Setting working directory back to original >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    setwd( originalWorkingDir )

    return( invisible(NULL) )

}
# <<
##################################### secondaryStructureCalculator() ###############################################################
####################################################################################################################################
