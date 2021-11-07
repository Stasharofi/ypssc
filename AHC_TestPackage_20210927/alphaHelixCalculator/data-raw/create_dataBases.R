
####################################################################################################################################
##################################### Reading the database #########################################################################
# >>

library(dplyr)

# `dataBase_small` database >>

dataBase_small    = read.csv("C:/Users/Shashank/Dropbox/Projects/Transmembrane-alpha-helix-calculator/dataBases/dataBase_small.csv")
dataBase_small$id = sub( "(\\|.*)", "",
                         sub( "(sp\\|)", "",
                              dataBase_small$id ) )

# `Alpha`, `beta`, `chain` databases >>

dataBase_alpha    = filter( dataBase_small, q3 == 'H' )
dataBase_beta     = filter( dataBase_small, q3 == 'E' )
dataBase_chain    = filter( dataBase_small, q3 == 'C' )

# `dataBase_numOfAA` database >>

numOfProteinsInDatabase = unique(dataBase_small$id)
i  = 1
aa = vector()
bb = vector()

pb_3 = winProgressBar(  title = "progress bar",
                        min   = 0,
                        max   = length(numOfProteinsInDatabase),
                        width = 300 )

for( i in 1 : length(numOfProteinsInDatabase) ) {

    item = numOfProteinsInDatabase[i]
    a    = nrow( filter(dataBase_small, id == item) )
    aa   = c(aa, a)
    bb   = c(bb, item)
    setWinProgressBar( pb_3, i, title = paste( 'Processing the DataBase ',
                                               '    ',
                                               round( i/length(numOfProteinsInDatabase)*100, 0 ),
                                               "% done") )
}

close(pb_3)

dataBase_numOfAA = data.frame( id = bb, numberofAA = aa )

# Saving `dataBase_small`, `dataBase_alpha`, `dataBase_beta`,`dataBase_chain`, `dataBase_numOfAA` database >>

  # Note: If you have more than one internal data, you must save them together.

usethis::use_data( dataBase_small,
                   dataBase_alpha,
                   dataBase_beta,
                   dataBase_chain,
                   dataBase_numOfAA,
                   overwrite = TRUE,
                   internal  = TRUE )

# <<
##################################### Reading the database #########################################################################
####################################################################################################################################
