################################
library(tidyverse) 
library(readxl)
library(readr)
library(dplyr)
#### here we have a list of selected proteins in each phase that we are looking for 
####I put aq, org alternatively and the progrma will read each list and search it in
####corresponding data base and peptide list

lop_tot<-read_excel("lop_tot.xlsx")
#this is the data base containing the part of protein that are TM, it hase been acuired from TMHMM
protter<-read_excel("protter.xlsx")
# this is a list of peptides identified in MaxQaunt
sample_totl<-read.csv('peptides.csv')
names(sample_totl)
#I am just removing the columns that I don't need for this caclculation
sample_totl[,c(2:34,36,39:51,56:76)]<-NULL
#I am renaming the column that contain proteins from Proteins to name;because this progrma recognizes
#the column name of "name"for proteins
sample_totl<-rename(sample_totl,name=Proteins)
#renaming Start.position to satrt
sample_totl<-rename(sample_totl,start=Start.position)
#renaming End.position to end
sample_totl<-rename(sample_totl,end=End.position)
names(sample_totl)
# I am putting the list of names of samples from 'prptides' file to a vector to be used
#later for naming the results accordingly
names<-names(sample_totl)
#I am filtering the rows of peptides that are TMhelix for the data base
protter<-filter(protter, position=="TMhelix")
names
#in this loop I am dreading the list of 2 of 3 proteins for each phase
# and telling the program to use it for all related replicates for that sample and phase
tt<-4
j=1
for(j in 1:4){
  tt<-tt+1
  sample<-filter(sample_totl, sample_totl[,tt]>0)
print(j)  
 
lop<-lop_tot[,j]
lop<-rename(lop,name=names[tt])
    
#in this loop we are filtering the peptides of a protein form the database and comparing it 
#with the peptides foud in sample to calculate coverage
lop_rows<-nrow(lop)
i=1
y_3<-vector()
y_2<-vector()
for(i in 1:lop_rows){
  
  item<-lop$name[i]
  
  protter_1<-filter(protter, name==item)
  protter_1_rows<-nrow(protter_1)
  sample_1<-filter(sample, name==item)
  ind<-order(sample_1$start)
  sample_1<-data.frame(
    name=sample_1$name[ind],start=sample_1$start[ind],
    end=sample_1$end[ind])
  sample_1_rows<-nrow(sample_1)
  
  if(protter_1_rows>0 &sample_1_rows>0){
#I am creating a database of start and end positions of TM proteins from the main database
  y<-vector()
  for(i in 1:protter_1_rows){
    a<-protter_1$start[i]
    b<-protter_1$end[i]
    x<-seq(a:b)
    x<-x+a-1
    z<-c(x,y)
    y<-z
    x<-vector()
  }
  z<-sort(z)
#I am crating a list of peptide positions in the sample to be compared with the previous database
  y_1<-vector()
  for(i in 1:sample_1_rows){
    a_1<-sample_1$start[i]
    b_1<-sample_1$end[i]
    x_1<-seq(a_1:b_1)
    x_1<-x_1+a_1-1
    z_1<-c(x_1,y_1)
    y_1<-z_1
    x_1<-vector()
  }

  z_1<-sort(z_1)
  f<-z%in%z_1
  g<-sum(f)
  length<-length(z)
  coverage<-g/length

  x_2<-item
  z_2<-c(x_2,y_2)
  y_2<-z_2
  x_2<-vector()
 
  x_3<-coverage
  z_3<-c(x_3,y_3)
  y_3<-z_3
  x_3<-vector()
  
  }
}
#after ordering the data frame created, I am writing the results
results<-data.frame(protein_name=z_2,coverage=z_3)
results<-results[order(-results$coverage),]
write.csv(results,file=paste(names[tt],'_results.csv'))
}
##################################
