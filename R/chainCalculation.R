
####################################################################################################################################
####################################################################################################################################
# >>
#' @title chainCalculation
#' @param df Dataframe for input file data.
#' @param sampleNames Names of the samples found in the input file.
#' @param sampleNamesUpdate Updated names of the samples found in the input file.
#' @param dateTimeCurrent Date and time at the time the simulation began.
#' @noRd
# <<
####################################################################################################################################
####################################################################################################################################


####################################################################################################################################
##################################### chainCalculation() ###########################################################################
# >>
chainCalculation <- function( df, sampleNames, sampleNamesUpdate, dateTimeCurrent ) {

    dataBase_chain   = select(dataBase_chain,c(1,2))
    dataBase_reduced = dataBase_chain

    num_Pro_caa    = unique(dataBase_reduced$id)
    protein        = vector()
    num_caa_pro_DB = vector()

    pb_1  =  winProgressBar( title = "progress bar",
                             min   = 0,
                             max   = length(num_Pro_caa),
                             width = 300)
    i = 1
    for( i in 1 : length(num_Pro_caa) ) {

        item                = num_Pro_caa[i]
        proteins            = filter( dataBase_chain, id == item )
        num_caa_pro_DB_temp = length(proteins$id)
        num_caa_pro_DB      = c(num_caa_pro_DB_temp,num_caa_pro_DB)
        protein             = c(unique(proteins$id),protein)
        proteins            = vector()
        num_caa_pro_DB_temp = vector()

        setWinProgressBar( pb_1, i,
                           title = paste( 'chain calculation for database     ',
                                          round( i/length(num_Pro_caa)*100, 0 ),
                                          "% done") )
    }

    close(pb_1)

    # Calculating the number of amino acids for chain ####

    caa              = data.frame( id      = protein,
                                   num_caa = num_caa_pro_DB )
    cal_for_database = left_join(  dataBase_numOfAA,
                                   caa,
                                   by = 'id' )

    Sys.sleep(0.5)

    # Samples ####

    i = 1
    for( i in 1 : length(sampleNames) ) {

        temp = which( names(df) == sampleNames[i] )

        # Peptides in the sample >>

        sample_peptides = filter( df, df[,temp] > 0 )
        write.csv( sample_peptides,
                   gsub( " ", "_", paste0( dateTimeCurrent,
                                           " ", 'List of peptides in',
                                           sampleNamesUpdate[i],
                                           '.csv' ), fixed = FALSE ),
                   row.names = FALSE )

        sample = paste( as.character(sampleNamesUpdate[i]), '_ peptides' )
        assign( sample, sample_peptides )

        # Proteins in the sample >>

        sample_proteins = unique(sample_peptides$Proteins)
        write.csv( sample_proteins,
                   gsub( " ", "_", paste0( dateTimeCurrent,
                                           " ", 'List of proteins in',
                                           sampleNamesUpdate[i],
                                           '.csv' ), fixed = FALSE ),
                   row.names = FALSE )

        sample = paste( as.character(sampleNamesUpdate[i]), '_ proteins' )
        assign( sample, sample_proteins )

        # Calculating beta-sheet coverage for samples >>

        proteins_in_s = vector()
        aa_in_s       = vector()
        caa_in_s      = vector()

        pb_2  =  winProgressBar( title = "progress bar",
                                 min   = 0,
                                 max   = length(sample_proteins),
                                 width = 300 )

        j = 1
        for( j in 1 : length(sample_proteins) ) {

            item      = sample_proteins[j]
            Pro_chunk = filter( sample_peptides, sample_peptides$Proteins == item )

            k = 1
            list_aa_s = vector()

            for( k in 1 : length(Pro_chunk$Proteins) ) {

                start          = Pro_chunk$Start.position[k]
                end            = Pro_chunk$End.position[k]
                list_aa_s_temp = seq(start:end)
                list_aa_s_temp = list_aa_s_temp+start-1
                list_aa_s      = c( list_aa_s_temp, list_aa_s )
                list_aa_s_temp = vector()

            }

            proteins_temp = item
            proteins_in_s = c( proteins_temp, proteins_in_s )
            proteins_temp = vector()

            aa_in_s_temp  = length( unique(list_aa_s) )
            aa_in_s       = c( aa_in_s_temp, aa_in_s )
            aa_in_s_temp  = vector()

            protein_chunk_dataBase = filter( dataBase_reduced, id == item )

            caa_in_s_temp = unique(list_aa_s)%in%protein_chunk_dataBase$n
            caa_in_s_temp = sum(caa_in_s_temp)
            caa_in_s      = c( caa_in_s_temp, caa_in_s )
            caa_in_s_temp = vector()

            results = data.frame( id                             = proteins_in_s,
                                  num_amino_acids_in_sample      = aa_in_s,
                                  num_chain_amino_acids_in_sample= caa_in_s )

            results = left_join( results, cal_for_database, by = 'id' )

            setWinProgressBar( pb_2, j,
                               title = paste( 'chain calculation for ',
                                              sampleNames[i],
                                              '    ',
                                              round( j/length(sample_proteins)*100, 0 ),
                                              "% done") )

        }

        write.csv( results,
                   gsub( " ", "_", paste0( dateTimeCurrent,
                                           " ", "chain analysis of",
                                           sampleNamesUpdate[i],
                                           ".csv" ), fixed = FALSE ),
                   row.names = FALSE )

        close(pb_2)

    }

    return( invisible(NULL) )
}
# <<
##################################### chainCalculation() ###########################################################################
####################################################################################################################################
