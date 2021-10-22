
####################################################################################################################################
##################################### auxil functions ##############################################################################
####################################################################################################################################


####################################################################################################################################
####################################################################################################################################
# >>
#' @title readFileInput
#' @param pathFileInput input file path from which bla bla bla xxxx xxxx xxxx xxxx xxxx xxxx xxxx.
# <<
####################################################################################################################################
####################################################################################################################################


####################################################################################################################################
##################################### readFileInput() ##############################################################################
# >>
readFileInput <- function( pathFileInput ) {

    # Reading csv file >>

    df = read.csv( pathFileInput )

    # Removing the columns that are not needed and finding the columns containing sample information >>

    df = df[ , -which( names(df) %in% c(    "Sequence","N.term.cleavage.window",
                                            "C.term.cleavage.window","Amino.acid.before",
                                            "First.amino.acid","Second.amino.acid",
                                            "Second.last.amino.acid","Last.amino.acid",
                                            "Amino.acid.after","A.Count","R.Count","N.Count",
                                            "D.Count","C.Count","Q.Count","E.Count",
                                            "G.Count","H.Count","I.Count","L.Count",
                                            "K.Count","M.Count","F.Count","P.Count",
                                            "S.Count","T.Count","W.Count","Y.Count",
                                            "V.Count","U.Count","O.Count","Length",
                                            "Missed.cleavages","Mass",
                                            "Leading.razor.protein","Gene.names",
                                            "Protein.names","Unique..Groups.",
                                            "Unique..Proteins.","Charges","PEP",
                                            "Score","Experiment.ST168.THF.A",
                                            "Experiment.ST168.THF.O",
                                            "Experiment.ST169.DMSO.A","Experiment.ST169.DMSO.O",
                                            "Experiment.ST170.But.A","Experiment.ST170.But.O",
                                            "Intensity.ST168.THF.A",
                                            "Intensity.ST168.THF.O",
                                            "Intensity.ST169.DMSO.A",
                                            "Intensity.ST169.DMSO.O",
                                            "Intensity.ST170.But.A",
                                            "id","Protein.group.IDs","Mod..peptide.IDs",
                                            "Evidence.IDs","MS.MS.IDs","Best.MS.MS",
                                            "Oxidation..M..site.IDs","Taxonomy.IDs",
                                            "MS.MS.Count" ) ) ]

    names = names(df)

    sampleNames         = names[ grepl("Intensity.", names) ]
    sampleNamesUpdate   = gsub( '\\.|Intensity.', ' ', sampleNames )
    names_list          = vector()
    i  = 1
    pb = winProgressBar( title = "progress bar",
                         min   = 0,
                         max   = length(sampleNames),
                         width = 300 )

    for ( i in 1 : length(sampleNames) ) {
        temp       = paste(sampleNamesUpdate[i],' \n \n ')
        temp
        names_list = paste( names_list, temp )
        Sys.sleep(0.9)
        setWinProgressBar( pb, i, title = paste( sampleNamesUpdate[i], '    ', round(i/length(sampleNames)*100, 0), "% done") )
    }

    close(pb)

    # Conformation about sample names from user >>

    sampleNameConfirmation = dlgInput(paste("Identified sample names in the uploaded file:\n \n \n", names_list,
                                            "\nIf it is correct, please enter 'Yes'"))$res
    class(sampleNameConfirmation)

    if ( sampleNameConfirmation == "yes" ) {

        tkmessageBox( title   = "Message",
                      message = "Your analysis in in progress",
                      icon    = "info",
                      type    = "ok" )
        as.character(names(df))

    } else if( sampleNameConfirmation=="no") {

        tkmessageBox( title   = "Message",
                      message = "Please put a 'samples.CSV' file containing just sample names in one column",
                      icon    = "info",
                      type    = "ok" )

        Sample_names = read.csv('samples.csv')
    }

    # Returning multiple variables as a R-list >>

    dataFileInput                   = list()
    dataFileInput$df                = df
    dataFileInput$sampleNames       = sampleNames
    dataFileInput$sampleNamesUpdate = sampleNamesUpdate

    return( dataFileInput )
}
# <<
##################################### readFileInput() ##############################################################################
####################################################################################################################################


####################################################################################################################################
##################################### creatOutputDir() #############################################################################
# >>
creatOutputDir <- function( pathDirOutput ) {

    dateTimeCurrent = format( Sys.time(), "%Y%m%d_%H%M%S" )        # << get current date and time
    nameDirOutput   = paste0( "results_AHC_", dateTimeCurrent )    # << name of the output folder
    pathDirOutput   = paste0( pathDirOutput, "/", nameDirOutput )  # << path of the output folder
    dir.create( pathDirOutput )                                    # creating new folder for output files
    setwd( pathDirOutput )              # << setting working dir to "pathDirOutput" to write output files

    return( dateTimeCurrent )
}
##################################### creatOutputDir() #############################################################################
####################################################################################################################################


####################################################################################################################################
##################################### removeRows() #################################################################################
# >>
removeRows <- function( df, dateTimeCurrent ) {

    removeDoubious = dlgInput( paste0("Do you want to remove the rows containing doubious proteins?\n",
                                      "Rows that have 2 or more protiens assigned to one identified peptide are called doubious\n",
                                      "Answer with yes or no") )$res
    df = filter( df, !grepl( ';', df$Proteins) )
    write.csv( df, paste0( dateTimeCurrent, " ", 'df.csv' ), row.names = FALSE)

    removeReverse  = dlgInput( paste0("Do you want to remove rows that contains peptides that matched to decoy that has reverse ",
                                      "sequnce of real protein?\n",
                                      "Theses proteins are usually removed.\n",
                                      "Answer with yes or no") )$res
    df = filter( df, !grepl( '\\+', df$Reverse) )

    removeReverse  = dlgInput( paste0("Do you want to remove rows that contains peptides that are showing signs of contamination?\n",
                                      "Theses proteins are usually removed.\n",
                                      "Answer with yes or no") )$res
    df = filter( df, !grepl( '\\+', df$Potential.contaminant) )

    removeReverse  = dlgInput( paste("Do you want to remove rows that contains peptides that are not showing any intensity?\n",
                                     "Theses proteins are usually removed.\n",
                                     "Answer with yes or no") )$res
    df = filter( df, df$Intensity > 0 )

    return( df )
}
##################################### removeRows() #################################################################################
####################################################################################################################################


####################################################################################################################################
##################################### auxil functions ##############################################################################
####################################################################################################################################
