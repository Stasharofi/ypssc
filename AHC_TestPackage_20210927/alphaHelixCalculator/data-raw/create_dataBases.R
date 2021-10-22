
####################################################################################################################################
##################################### Reading the database #########################################################################
# >>

# Main database >>

dataBase_small    = read.csv("C:/Users/Shashank/Dropbox/Projects/Transmembrane-alpha-helix-calculator/dataBases/dataBase_small.csv")
dataBase_small$id = sub( "(\\|.*)", "",
                         sub( "(sp\\|)", "",
                              dataBase_small$id ) )
usethis::use_data( dataBase_small, overwrite = TRUE )      # << writing main database

# Alpha, beta, chain databases >>

dataBase_alpha    = filter( dataBase_small, q3 == 'H' )
dataBase_beta     = filter( dataBase_small, q3 == 'E' )
dataBase_chain    = filter( dataBase_small, q3 == 'C' )

usethis::use_data( dataBase_alpha, overwrite = TRUE )      # << writing alpha database
usethis::use_data( dataBase_beta,  overwrite = TRUE )      # << writing beta  database
usethis::use_data( dataBase_chain, overwrite = TRUE )      # << writing chain database

# <<
##################################### Reading the database #########################################################################
####################################################################################################################################
