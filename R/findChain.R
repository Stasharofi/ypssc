
####################################################################################################################################
####################################################################################################################################
# >>
#' @title Chain Calculator
#' @description Form bottom-up proteomics data of proteins (peptides), this function determines the sections of proteins (in
#'   percentage) with primary, structure.
#' @param pathFileInput Path of the input csv file generated from MaxQuant. \cr
#'   \cr
#'   MaxQuant is a quantitative proteomics software designed to analyze large mass-spectrometric data. The input of MaxQuant is a
#'   raw file (.raw) from high-resolution mass spectrometers. After analysis of the raw file in MaxQuant, the program generates a
#'   folder named “combined”. \cr
#'   \cr
#'   In this folder there is another folder named “txt” which contains many files with text format (.txt). One of the files called
#'   “peptides” which is the input of the YPSSC to calculate secondary structures. YPSSC has been designed such a way that can
#'   analyzed and extract information regarding the sample regardless of the name that user chosen for the sample.
#' @param pathDirOutput Path of the directory to which the output files will be generated.
#' @return The output of the program is a csv file (.csv) that contains 5 columns, and the number of rows depends on the number of
#'   proteins in the sample. \cr
#'   \cr
#'   First column contains the ID of the identified alpha-helix proteins in the sample, second column contains the number of identified
#'   amino acids from the corresponding protein, third column contains number of identified amino acids with secondary structure,
#'   fourth column contains the number of amino acids that the protein originally has in the SSDYP, and fifth column contains the
#'   number of amino acids in chain structure that the protein originally has in the SSDYP. \cr
#'   \cr
#'   These columns should provide all information that the user needs to know about the protein and its structural information as
#'   well as structural information about the parts of the protein that has been identified in the sample. \cr
#'   \cr
#'   In addition, it also generates 4 more '.csv' files. \cr
#'   1. The no. of proteins found in the sample. \cr
#'   2. The no. of peptides found in the sample. \cr
#'   3. The no. of amino acids for each protein in database. \cr
#'   4. It is the input file from MaxQuant that's been cleaned up for the sole purpose of calculating secondary structures.
#' @examples
#' \dontrun{
#' findChain( pathFileInput = "<somePathInputFile>",
#'            pathDirOutput = "<somePathOutputFolder>")
#' }
#' @seealso [`findSecondary`], [`findAlpha`], [`findBeta`]
#' @export
# <<
####################################################################################################################################
####################################################################################################################################


####################################################################################################################################
##################################### findChain() ##################################################################################
# >>
findChain <- function( pathFileInput = NULL,
                       pathDirOutput = NULL ) {

    # Begin >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    print("Started")
    startTime = Sys.time()

    # Getting current working directory >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    originalWorkingDir = getwd()

    # Checking if 'pathFileInput' is provided >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    if ( is.null( pathFileInput ) ) {
        writeLines("Please provide path of the input file as argument `pathFileInput`")
        return( invisible(NULL) )
    }

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

    dateTimeCurrent = creatOutputDir( pathDirOutput, "chain" )

    # Removing the rows that are not needed >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    df = removeRows( df, dateTimeCurrent )

    # Writing `dataBase_numOfAA` >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    write.csv( dataBase_numOfAA,
               paste0( dateTimeCurrent, "_dataBase_numOfAA.csv" ),
               row.names = FALSE )

    # Chain calculation for dataBase >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    chainCalculation     ( df, sampleNames, sampleNamesUpdate, dateTimeCurrent )

    # End >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    endTime   = Sys.time()
    timeTaken = endTime - startTime
    print( paste0( "Time taken for the ypssc run: ", format(timeTaken) ) )

    # Setting working directory back to original >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    setwd( originalWorkingDir )

    return( invisible(NULL) )

}
# <<
##################################### findChain() ##################################################################################
####################################################################################################################################
