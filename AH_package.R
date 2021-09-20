install.packages('svMisc')
library(tidyverse)
library(readxl)
library(stringr)
library(eulerr)
library(ggplot2)
library(Peptides)
library(utils)
library(svDialogs)
library(tcltk)
library(svMisc)

# setting the working directory to downloads so we can download the database from Dropbox

setwd('C:/Users/sxt4104/OneDrive - University of Texas at Arlington/Documents/AH_Package')

# reading the reduced database file form GitHub
dataBase_reduced<-read.csv("https://raw.githubusercontent.com/Stasharofi/Transmembrane-alpha-helix-calculator/main/dataBase_reduced.csv")

# reading the database containing number of AA for each protein file form GitHub
dataBase_numOfAA<-read.csv("https://raw.githubusercontent.com/Stasharofi/Transmembrane-alpha-helix-calculator/main/dataBase_numOfAA.csv")

# reading the sample form the user

df<-read.csv("peptides_second rep.csv")
df[c(61:69)]<-NULL
names<-names(df)
sampleNames<-names[grepl("Intensity.", names)]
sampleNamesUpdate<-gsub('\\.|Intensity',' ',sampleNames)

names_list<-vector()
i<-1
pb <- winProgressBar(title = "progress bar", min = 0,
                     max = length(sampleNames), width = 300)
for(i in 1:length(sampleNames)){
  temp<-paste(sampleNamesUpdate[i],' \n \n ')
  temp
  names_list<-paste(names_list,temp)
  Sys.sleep(0.9)
  setWinProgressBar(pb, i, title=paste(sampleNamesUpdate[i],'    ', round(i/length(sampleNames)*100, 0),
                                        "% done"))
}
close(pb)



sampleNameConfirmation<- dlgInput(paste("Identified sample names in the uploaded file:\n \n \n", names_list,
                                         "\nIf it is correct, please enter 'Yes'"))$res
class(sampleNameConfirmation)
if(sampleNameConfirmation=="yes"){
  print('ok')
  tkmessageBox(title = "Message",
               message = "Your analysis in in progress", icon = "info", type = "ok")
 as.character(names(df))

}else if(sampleNameConfirmation=="no"){
  print('notok')
  tkmessageBox(title = "Message",
               message = "Please put a 'samples.CSV' file containing just sample names in one column", icon = "info", type = "ok")

  Sample_names<-read.csv('samples.csv')
}

# removing the rows that are not needed

removeDoubious <- dlgInput(paste("Do you want to remove the rows containing doubious proteins?\n
Rows that have 2 or more protiens assigned to one identified peptide are called doubious\n
                                 Answer with yes or no"))$res
df<-filter(df, !grepl(';',Proteins ))
write.csv(df,'df.csv',row.names = FALSE)

removeReverse <- dlgInput(paste("Do you want to remove rows that contains peptides that matched to decoy that has reverse sequnce of real protein?\n
Theses proteins are usually removed.\n
                                 Answer with yes or no"))$res
df<-filter(df, !grepl('\\+',Reverse ))

removeReverse <- dlgInput(paste("Do you want to remove rows that contains peptides that are showing signs of contamination?\n
Theses proteins are usually removed.\n
                                 Answer with yes or no"))$res
df<-filter(df, !grepl('\\+',Potential.contaminant ))

removeReverse <- dlgInput(paste("Do you want to remove rows that contains peptides that are not showing any intensity?\n
Theses proteins are usually removed.\n
                                 Answer with yes or no"))$res
df<-filter(df, Intensity>0)

# alpha helix calculation for dataBase

num_Pro_aaa<-unique(dataBase_reduced$id)
protein<-vector()
num_aaa_pro_DB<-vector()


pb_1 <- winProgressBar(title = "progress bar", min = 0,
                     max = length(num_Pro_aaa), width = 300)
i<-1
for(i in 1:length(num_Pro_aaa)){
  item<-num_Pro_aaa[i]
  proteins<-filter(dataBase_reduced, id==item)
  num_aaa_pro_DB_temp<-length(proteins$id)
  num_aaa_pro_DB<-c(num_aaa_pro_DB_temp,num_aaa_pro_DB)
  protein<-c(unique(proteins$id),protein)
  proteins<-vector()
  num_aaa_pro_DB_temp<-vector()

  setWinProgressBar(pb_1, i, title=paste('Alpha-helix calculation for database     ', round(i/length(num_Pro_aaa)*100, 0),
                                        "% done"))
}
close(pb_1)
aaa<-data.frame(id=protein,num_aaa=num_aaa_pro_DB)
cal_for_database<-left_join(dataBase_numOfAA,aaa, by='id')

Sys.sleep(0.5)

# samples



i<-1
for(i in 1:length(sampleNames)){
  temp<-which(names(df)==sampleNames[i])

  # Peptides in the sample

  sample_peptides<-filter(df,df[,temp]>0)
  write.csv(sample_peptides,paste('List of peptides in',sampleNamesUpdate[i],'.csv'),row.names = FALSE)

  sample<-paste(as.character(sampleNamesUpdate[i]),'_ peptides')
  assign(sample,sample_peptides)

  # Proteins in the sample

  sample_proteins<-unique(sample_peptides$Proteins)
  write.csv(sample_proteins,paste('List of proteins in',sampleNamesUpdate[i],'.csv'),row.names = FALSE)

  sample<-paste(as.character(sampleNamesUpdate[i]),'_ proteins')
  assign(sample,sample_proteins)

  # calculating alpha helix coverage for samples

  startTime<-Sys.time()

  proteins_in_s<-vector()
  aa_in_s<-vector()
  aaa_in_s<-vector()

  pb_2 <- winProgressBar(title = "progress bar", min = 0,
                         max = length(sample_proteins), width = 300)

  j<-1
  for( j in 1:length(sample_proteins)){
    item<-sample_proteins[j]
    Pro_chunk<-filter(sample_peptides,Proteins==item)

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

    aaa_in_s_temp<-unique(list_aa_s)%in%protein_chunk_dataBase$n
    aaa_in_s_temp<-sum(aaa_in_s_temp)
    aaa_in_s<-c(aaa_in_s_temp,aaa_in_s)
    aaa_in_s_temp<-vector()

    results<-data.frame(id=proteins_in_s,
                        num_amino_acids_in_sample=aa_in_s,
                        num_alpha_amino_acids_in_sample=aaa_in_s)
    results<-left_join(results,cal_for_database,by='id')

    write.csv(results,paste('alpha_helix analysis of',sampleNamesUpdate[i],'.csv'),row.names = FALSE)

    setWinProgressBar(pb_2, j, title=paste('Alpha-helix calculation for ',sampleNames[i],'    ', round(j/length(sample_proteins)*100, 0),
                                           "% done"))

  }
  close(pb_2)

}

endTime<-Sys.time()
endTime-startTime


















