
####################################################################################################################################
####################################################################################################################################
# >>
#' @title alphaHelixCalculator xxx xxx xxx
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
#' @return xxxx xxxx xxxx xxxx xxxx xxxx xxxx xxxx
#' @details DETAILS
#' @examples
#' \dontrun{
#' alphaHelixCalculator( pathFileInput    = "<someInputFileName>",
#'                       pathDirOutput = "<someOutputFolderName>")
#' }
#' @export
# <<
####################################################################################################################################
####################################################################################################################################


####################################################################################################################################
##################################### alphaHelixCalculator() #######################################################################
# >>
alphaHelixCalculator <- function( pathFileInput = "C:/Users/Shashank/Desktop/peptides_second rep.csv",
                                  pathDirOutput = "C:/Users/Shashank/Downloads" ) {

    print("Started")

    originalWorkingDir = getwd()   # << getting original current working directory

    # Checking if 'pathDirOutput' is provided ###################################################################################

    if ( is.null( pathDirOutput ) ) {
        pathDirOutput = getwd()
    }

    # Creating separate dataBases for alpha, beta and chain ########################################################################

    dataBase_alpha    = dataBase_small
    dataBase_alpha$id = sub( "(sp\\|)", "",
                             dataBase_alpha$id )
    dataBase_alpha$id = sub( "(\\|.*)", "",
                            dataBase_alpha$id )
    dataBase_alpha = filter(dataBase_alpha,
                            q3=='H')

    dataBase_beta = dataBase_small
    dataBase_beta$id = sub("(sp\\|)","",
                           dataBase_beta$id)
    dataBase_beta$id = sub("(\\|.*)","",
                           dataBase_beta$id)
    dataBase_beta = filter(dataBase_beta,
                           q3=='E')

    dataBase_chain = dataBase_small
    dataBase_chain$id = sub("(sp\\|)","",
                            dataBase_chain$id)
    dataBase_chain$id = sub("(\\|.*)","",
                            dataBase_chain$id)
    dataBase_chain = filter(dataBase_chain,
                            q3=='C')

    # Reading the sample form the user##############################################################################################

    df = read.csv( pathFileInput )

    # Removing the columns that are not needed and finding the columns containing sample information ###############################

    df = df[ , -which(names(df) %in% c("Sequence","N.term.cleavage.window",
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
                                       "MS.MS.Count"))]

    names               = names(df)
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

    # Removing the rows that are not needed ########################################################################################

    dateTimeCurrent  = format( Sys.time(), "%Y%m%d_%H%M%S" )
    nameFolderOutput = paste0( "results_AHC_", dateTimeCurrent )         # << name of the output folder
    pathDirOutput = paste0( pathDirOutput, "/", nameFolderOutput ) # << path of the output folder
    dir.create( pathDirOutput )                                       # creating new folder for output files
    setwd( pathDirOutput )                      # << setting working dir to "pathDirOutput" to write output files

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

    # Relaying the alpha database to data base reduced for local calculation #######################################################

    typeOfAnalysis =  dlgInput( paste( "What type of analysis are you willing to perform?\n \n \n",
                                       "1) Alpha_helix\n\n",
                                       "2) Beta-sheet\n\n",
                                       "3) Chain\n\n",
                                       "4) Alpha_helix and Beta-sheet\n\n",
                                       "5) Alpha_helix and Chain\n\n",
                                       "6) Beta_sheet and chain\n\n",
                                       "7) 1,2 and 3\n\n" ) )$res

    # Alpha helix calculation for dataBase #########################################################################################

    if( typeOfAnalysis == "1" |
        typeOfAnalysis == "4" |
        typeOfAnalysis == "5" |
        typeOfAnalysis == "7" ) {

        dataBase_alpha   = select(dataBase_alpha,c(1,2))
        dataBase_reduced = dataBase_alpha

        num_Pro_aaa    = unique(dataBase_reduced$id)
        protein        = vector()
        num_aaa_pro_DB = vector()

        pb_1 = winProgressBar( title = "progress bar",
                               min   = 0,
                               max   = length(num_Pro_aaa),
                               width = 300 )
        i = 1
        for( i in 1 : length(num_Pro_aaa) ) {
            item                = num_Pro_aaa[i]
            proteins            = filter(dataBase_reduced, id == item)
            num_aaa_pro_DB_temp = length(proteins$id)
            num_aaa_pro_DB      = c(num_aaa_pro_DB_temp,num_aaa_pro_DB)
            protein             = c(unique(proteins$id),protein)
            proteins            = vector()
            num_aaa_pro_DB_temp = vector()

            setWinProgressBar( pb_1, i, title = paste( 'Alpha-helix calculation for database     ',
                                                       round(i/length(num_Pro_aaa)*100, 0),
                                                       "% done") )
        }
        close(pb_1)

        # calculating the number of amino acids for alpha ####

        dataBase_small_2    = dataBase_small
        dataBase_small_2$id = sub("(sp\\|)","",
                                  dataBase_small_2$id)
        dataBase_small_2$id = sub("(\\|.*)","",
                                  dataBase_small_2$id)

        numOfProteinsInDatabase = unique(dataBase_small$id)

        if(typeOfAnalysis=="1"|
           typeOfAnalysis=="4"|
           typeOfAnalysis=="5"|
           typeOfAnalysis=="7"){

            i  = 1
            aa = vector()
            bb = vector()

            pb_3  =  winProgressBar(title = "progress bar",
                                    min = 0,
                                    max = length(numOfProteinsInDatabase),
                                    width = 300)

            for(i in 1:length(numOfProteinsInDatabase)){
                item = numOfProteinsInDatabase[i]
                a = nrow(filter(dataBase_small_2, id==item))
                aa = c(aa,a)
                bb = c(bb,item)
                setWinProgressBar(pb_3, i, title=paste('Processing the DataBase ',
                                                       '    ',
                                                       round(i/length(numOfProteinsInDatabase)*100, 0),
                                                       "% done"))
            }
            close(pb_3)
        }

        dataBase_numOfAA = data.frame(id=bb,numberofAA=aa)

        write.csv(dataBase_numOfAA, "dataBase_numOfAA.csv",
                  row.names = FALSE)

        aaa              = data.frame( id      = protein,
                                       num_aaa = num_aaa_pro_DB )
        cal_for_database = left_join( dataBase_numOfAA,
                                      aaa,
                                      by = 'id' )

        Sys.sleep(0.5)

        # Samples ####

        i = 1
        for( i in 1 : length(sampleNames) ) {

            temp = which( names(df) == sampleNames[i] )

            # Peptides in the sample >>

            sample_peptides = filter( df, df[,temp] > 0 )
            write.csv( sample_peptides,
                       paste0( dateTimeCurrent,
                               " ", 'List of peptides in',
                               sampleNamesUpdate[i], '.csv' ),
                       row.names = FALSE )

            sample = paste( as.character(sampleNamesUpdate[i]), '_ peptides' )
            assign( sample, sample_peptides )

            # Proteins in the sample >>

            sample_proteins = unique( sample_peptides$Proteins )
            write.csv( sample_proteins,
                       paste0( dateTimeCurrent,
                               " ", 'List of proteins in',
                               sampleNamesUpdate[i], '.csv' ),
                       row.names = FALSE )

            sample = paste( as.character(sampleNamesUpdate[i]), '_ proteins' )
            assign( sample, sample_proteins )

            # Calculating alpha helix coverage for samples >>

            startTime = Sys.time()

            proteins_in_s = vector()
            aa_in_s       = vector()
            aaa_in_s      = vector()

            pb_2 = winProgressBar( title = "progress bar",
                                   min   = 0,
                                   max   = length(sample_proteins),
                                   width = 300 )

            j = 1
            for ( j in 1 : length(sample_proteins) ) {

                item      = sample_proteins[j]
                Pro_chunk = filter(sample_peptides, sample_peptides$Proteins == item)

                k         = 1
                list_aa_s = vector()

                for( k in 1 : length(Pro_chunk$Proteins) ) {

                    start = Pro_chunk$Start.position[k]
                    end   = Pro_chunk$End.position[k]
                    list_aa_s_temp = seq(start:end)
                    list_aa_s_temp = list_aa_s_temp+start-1
                    list_aa_s      = c(list_aa_s_temp,list_aa_s)
                    list_aa_s_temp = vector()

                }

                proteins_temp = item
                proteins_in_s = c(proteins_temp,proteins_in_s)
                proteins_temp = vector()

                aa_in_s_temp = length(unique(list_aa_s))
                aa_in_s      = c(aa_in_s_temp,aa_in_s)
                aa_in_s_temp = vector()

                protein_chunk_dataBase = filter(dataBase_reduced, id==item)

                aaa_in_s_temp = unique(list_aa_s)%in%protein_chunk_dataBase$n
                aaa_in_s_temp = sum(aaa_in_s_temp)
                aaa_in_s      = c(aaa_in_s_temp,aaa_in_s)
                aaa_in_s_temp = vector()

                results = data.frame( id                              = proteins_in_s,
                                      num_amino_acids_in_sample       = aa_in_s,
                                      num_alpha_amino_acids_in_sample = aaa_in_s )

                results = left_join( results, cal_for_database, by = 'id' )

                write.csv( results,
                           paste0( dateTimeCurrent,
                                   " ", 'alpha_helix analysis of',
                                   sampleNamesUpdate[i], '.csv' ),
                           row.names = FALSE )

                setWinProgressBar( pb_2, j,
                                   title = paste( 'Alpha-helix calculation for',
                                                  sampleNames[i],
                                                  '    ',
                                                  round(j/length(sample_proteins)*100, 0),
                                                  "% done"))

            }

            close(pb_2)

        }

        if(typeOfAnalysis=="1"){

            setwd( originalWorkingDir )         # << setting working directory back to original

            endTime   = Sys.time()
            timeTaken = endTime - startTime
            print( paste0( "Time taken for the AHC run: ", format(timeTaken) ) )

            return( invisible(NULL) )
        }

    }


    # Beta-sheet calculation for dataBase ##########################################################################################

    if(typeOfAnalysis=="2"|
       typeOfAnalysis=="4"|
       typeOfAnalysis=="6"|
       typeOfAnalysis=="7"){

        dataBase_beta = select(dataBase_beta,c(1,2))
        dataBase_reduced = dataBase_beta

        num_Pro_baa = unique(dataBase_reduced$id)

        protein = vector()
        num_baa_pro_DB = vector()


        pb_1 = winProgressBar(title = "progress bar",
                              min = 0,
                              max = length(num_Pro_baa),
                              width = 300)
        i = 1
        for(i in 1:length(num_Pro_baa)){
            item = num_Pro_baa[i]
            proteins = filter(dataBase_reduced, id==item)
            num_baa_pro_DB_temp = length(proteins$id)
            num_baa_pro_DB = c(num_baa_pro_DB_temp,num_baa_pro_DB)
            protein = c(unique(proteins$id),protein)
            proteins  = vector()
            num_baa_pro_DB_temp = vector()

            setWinProgressBar(pb_1, i, title=paste('Beta-sheet calculation for database     ',
                                                   round(i/length(num_Pro_baa)*100, 0),
                                                   "% done"))
        }
        close(pb_1)

        # Calculating the number of amino acids for beta ####

        dataBase_small_2 = dataBase_small
        dataBase_small_2$id = sub("(sp\\|)","",
                                  dataBase_small_2$id)
        dataBase_small_2$id = sub("(\\|.*)","",
                                  dataBase_small_2$id)
        numOfProteinsInDatabase = unique(dataBase_small_2$id)

        if(typeOfAnalysis=="2"|
           typeOfAnalysis=="6"){

            i = 1
            aa = vector()
            bb = vector()

            pb_3  =  winProgressBar(title = "progress bar",
                                    min = 0,
                                    max = length(numOfProteinsInDatabase),
                                    width = 300)

            for(i in 1:length(numOfProteinsInDatabase)){
                item = numOfProteinsInDatabase[i]
                a = nrow(filter(dataBase_small_2, id==item))
                aa = c(aa,a)
                bb = c(bb,item)
                setWinProgressBar(pb_3, i,
                                  title=paste('Processing the DataBase ',
                                              '    ',
                                              round(i/length(numOfProteinsInDatabase)*100, 0),
                                              "% done"))
            }
            close(pb_3)
        }
        dataBase_numOfAA = data.frame(id=bb,numberofAA=aa)

        baa = data.frame(id=protein,
                         num_baa=num_baa_pro_DB)

        cal_for_database = left_join(dataBase_numOfAA,baa,
                                     by='id')


        Sys.sleep(0.5)

        # Samples ####

        i<-1
        for(i in 1:length(sampleNames)){
            temp<-which(names(df)==sampleNames[i])

            # Peptides in the sample >>

            sample_peptides = filter(df,df[,temp]>0)
            write.csv(sample_peptides,
                      paste('List of peptides in',
                            sampleNamesUpdate[i],
                            '.csv'),
                      row.names = FALSE)

            sample = paste(as.character(sampleNamesUpdate[i]),'_ peptides')
            assign(sample,sample_peptides)

            # Proteins in the sample >>

            sample_proteins = unique(sample_peptides$Proteins)

            write.csv(sample_proteins,
                      paste('List of proteins in',
                            sampleNamesUpdate[i],
                            '.csv'),
                      row.names = FALSE)

            sample = paste(as.character(sampleNamesUpdate[i]),'_ proteins')
            assign(sample,sample_proteins)

            # Calculating beta-sheet coverage for samples >>
            startTime = Sys.time()

            proteins_in_s = vector()
            aa_in_s = vector()
            baa_in_s = vector()

            pb_2  =  winProgressBar(title = "progress bar",
                                    min = 0,
                                    max = length(sample_proteins),
                                    width = 300)

            j = 1
            for( j in 1:length(sample_proteins)){
                item = sample_proteins[j]
                Pro_chunk = filter(sample_peptides,sample_peptides$Proteins==item)

                k<-1
                list_aa_s<-vector()
                for(k in 1:length(Pro_chunk$Proteins)){
                    start<-Pro_chunk$Start.position[k]
                    end<-Pro_chunk$End.position[k]
                    list_aa_s_temp<-seq(start:end)
                    list_aa_s_temp<-list_aa_s_temp+start-1
                    list_aa_s<-c(list_aa_s_temp,list_aa_s)
                    list_aa_s_temp<-vector()


                }

                proteins_temp<-item
                proteins_in_s<-c(proteins_temp,proteins_in_s)
                proteins_temp<-vector()

                aa_in_s_temp<-length(unique(list_aa_s))
                aa_in_s<-c(aa_in_s_temp,aa_in_s)
                aa_in_s_temp<-vector()

                protein_chunk_dataBase<-filter(dataBase_reduced, id==item)

                baa_in_s_temp<-unique(list_aa_s)%in%protein_chunk_dataBase$n
                baa_in_s_temp<-sum(baa_in_s_temp)
                baa_in_s<-c(baa_in_s_temp,baa_in_s)
                baa_in_s_temp<-vector()

                results<-data.frame(id=proteins_in_s,
                                    num_amino_acids_in_sample=aa_in_s,
                                    num_beta_amino_acids_in_sample=baa_in_s)
                results<-left_join(results,cal_for_database,by='id')

                write.csv(results,
                          paste('beta-sheet analysis of',
                                sampleNamesUpdate[i],
                                '.csv'),
                          row.names = FALSE)

                setWinProgressBar(pb_2, j,
                                  title=paste('Beta-sheet calculation for ',
                                              sampleNames[i],
                                              '    ',
                                              round(j/length(sample_proteins)*100, 0),
                                              "% done"))

            }

            close(pb_2)
        }
        if(typeOfAnalysis=="2"|
           typeOfAnalysis=="4"){

            setwd( originalWorkingDir )         # << setting working directory back to original

            endTime   = Sys.time()
            timeTaken = endTime - startTime
            print( paste0( "Time taken for the AHC run: ", format(timeTaken) ) )

            return( invisible(NULL) )
        }

    }

    # Chain calculation for dataBase ###############################################################################################

    if(typeOfAnalysis=="3"|
       typeOfAnalysis=="5"|
       typeOfAnalysis=="6"|
       typeOfAnalysis=="7"){

        dataBase_chain = select(dataBase_chain,c(1,2))
        dataBase_reduced = dataBase_chain

        num_Pro_caa = unique(dataBase_reduced$id)

        protein = vector()
        num_caa_pro_DB = vector()

        pb_1  =  winProgressBar(title = "progress bar",
                                min = 0,
                                max = length(num_Pro_caa),
                                width = 300)
        i<-1
        for(i in 1:length(num_Pro_caa)){
            item = num_Pro_caa[i]
            proteins = filter(dataBase_chain, id==item)
            num_caa_pro_DB_temp = length(proteins$id)
            num_caa_pro_DB = c(num_caa_pro_DB_temp,num_caa_pro_DB)
            protein = c(unique(proteins$id),protein)
            proteins = vector()
            num_caa_pro_DB_temp = vector()

            setWinProgressBar(pb_1, i,
                              title=paste('chain calculation for database     ',
                                          round(i/length(num_Pro_baa)*100, 0),
                                          "% done"))
        }
        close(pb_1)

        # Calculating the number of amino acids for chain ####

        dataBase_small_2 = dataBase_small
        dataBase_small_2$id = sub("(sp\\|)","",
                                  dataBase_small_2$id)
        dataBase_small_2$id = sub("(\\|.*)","",
                                  dataBase_small_2$id)
        numOfProteinsInDatabase = unique(dataBase_small_2$id)

        if(typeOfAnalysis=="3"){

            i = 1
            aa = vector()
            bb = vector()

            pb_3  =  winProgressBar(title = "progress bar",
                                    min = 0,
                                    max = length(numOfProteinsInDatabase),
                                    width = 300)

            for(i in 1:length(numOfProteinsInDatabase)){
                item = numOfProteinsInDatabase[i]
                a = nrow(filter(dataBase_small_2, id==item))
                aa = c(aa,a)
                bb = c(bb,item)
                setWinProgressBar(pb_3, i,
                                  title=paste('Processing the DataBase ',
                                              '    ',
                                              round(i/length(numOfProteinsInDatabase)*100, 0),
                                              "% done"))
            }
            close(pb_3)
        }
        dataBase_numOfAA = data.frame(id=bb,numberofAA=aa)

        caa = data.frame(id=protein,
                         num_caa=num_caa_pro_DB)
        cal_for_database = left_join(dataBase_numOfAA
                                     ,caa,
                                     by='id')


        Sys.sleep(0.5)

        # Samples ####

        i = 1
        for(i in 1:length(sampleNames)){

            temp = which(names(df)==sampleNames[i])

            # Peptides in the sample >>

            sample_peptides = filter(df,df[,temp]>0)
            write.csv(sample_peptides,
                      paste('List of peptides in',
                            sampleNamesUpdate[i],
                            '.csv'),
                      row.names = FALSE)

            sample = paste(as.character(sampleNamesUpdate[i]),'_ peptides')
            assign(sample,sample_peptides)

            # Proteins in the sample >>

            sample_proteins = unique(sample_peptides$Proteins)
            write.csv(sample_proteins,
                      paste('List of proteins in',
                            sampleNamesUpdate[i],
                            '.csv'),
                      row.names = FALSE)

            sample<-paste(as.character(sampleNamesUpdate[i]),'_ proteins')
            assign(sample,sample_proteins)

            # Calculating beta-sheet coverage for samples >>

            startTime = Sys.time()

            proteins_in_s = vector()
            aa_in_s = vector()
            caa_in_s = vector()

            pb_2  =  winProgressBar(title = "progress bar",
                                    min = 0,
                                    max = length(sample_proteins),
                                    width = 300)

            j = 1
            for( j in 1:length(sample_proteins)){
                item = sample_proteins[j]
                Pro_chunk = filter(sample_peptides,sample_peptides$Proteins==item)

                k = 1
                list_aa_s = vector()
                for(k in 1:length(Pro_chunk$Proteins)){
                    start = Pro_chunk$Start.position[k]
                    end = Pro_chunk$End.position[k]
                    list_aa_s_temp = seq(start:end)
                    list_aa_s_temp = list_aa_s_temp+start-1
                    list_aa_s = c(list_aa_s_temp,list_aa_s)
                    list_aa_s_temp = vector()


                }

                proteins_temp = item
                proteins_in_s = c(proteins_temp,proteins_in_s)
                proteins_temp = vector()

                aa_in_s_temp = length(unique(list_aa_s))
                aa_in_s = c(aa_in_s_temp,aa_in_s)
                aa_in_s_temp = vector()

                protein_chunk_dataBase = filter(dataBase_reduced, id==item)

                caa_in_s_temp = unique(list_aa_s)%in%protein_chunk_dataBase$n
                caa_in_s_temp = sum(caa_in_s_temp)
                caa_in_s = c(caa_in_s_temp,caa_in_s)
                caa_in_s_temp = vector()

                results = data.frame(id=proteins_in_s,
                                     num_amino_acids_in_sample=aa_in_s,
                                     num_chain_amino_acids_in_sample=caa_in_s)
                results = left_join(results,cal_for_database,by='id')

                write.csv(results,
                          paste('chain analysis of',
                                sampleNamesUpdate[i],
                                '.csv'),
                          row.names = FALSE)

                setWinProgressBar(pb_2, j,
                                  title=paste('chain calculation for ',
                                              sampleNames[i],'    ',
                                              round(j/length(sample_proteins)*100, 0),
                                              "% done"))

            }

            close(pb_2)
        }

        if ( typeOfAnalysis == "3" |
             typeOfAnalysis == "5" |
             typeOfAnalysis == "6" |
             typeOfAnalysis == "7" ) {

            setwd( originalWorkingDir )         # << setting working directory back to original

            endTime   = Sys.time()
            timeTaken = endTime - startTime
            print( paste0( "Time taken for the AHC run: ", format(timeTaken) ) )

            return( invisible(NULL) )
        }
    }
}
# <<
##################################### alphaHelixCalculator() #######################################################################
####################################################################################################################################



####################################################################################################################################
##################################### betaSheetCalculator() ########################################################################
# >>
betaSheetCalculator = function( pathFileInput    = "C:/Users/Shashank/Desktop/peptides_second rep.csv",
                                pathDirOutput = "C:/Users/Shashank/Downloads" ) {

    print("Started")

    originalWorkingDir = getwd()   # << getting original current working directory

    # Checking if 'pathDirOutput' is provided ###################################################################################

    if ( is.null( pathDirOutput ) ) {
        pathDirOutput = getwd()
    }

    # Start Here >>>>>>>

}
# <<
##################################### betaSheetCalculator() ########################################################################
####################################################################################################################################



####################################################################################################################################
##################################### chainCalculator() ############################################################################
# >>
chainCalculator = function( pathFileInput    = "C:/Users/Shashank/Desktop/peptides_second rep.csv",
                            pathDirOutput = "C:/Users/Shashank/Downloads" ) {

    print("Started")

    originalWorkingDir = getwd()   # << getting original current working directory

    # Checking if 'pathDirOutput' is provided ###################################################################################

    if ( is.null( pathDirOutput ) ) {
        pathDirOutput = getwd()
    }

    # Start Here >>>>>>>


}
# <<
##################################### chainCalculator() ############################################################################
####################################################################################################################################
