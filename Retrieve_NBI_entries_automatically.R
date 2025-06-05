
library(rentrez)
library(openxlsx)


# read the meat-information: --------
data<-read.xlsx("/Home/johanna/Dokumente/1000_bull_genomes/Data/Run9-TAUINDraw-animalDistributionList-20210630.xlsx")


# extract NCBI information where ever possible: -------------------
# prepare he dataframe to store the NCBI information
data_add_ncbi=data

data_add_ncbi$ncbi_infraspecies<-vector(length=6191)
data_add_ncbi$ncbi_title<-vector(length=6191)
data_add_ncbi$ncbi_species<-vector(length=6191)

start=Sys.time()

for(i in 1:nrow(data_add_ncbi)){
  
  #make sure in there is an ID available
  if(!(is.na(data$BioSampleID[i]))){
    #extract the correct biosample information
    bp2=entrez_search(db="biosample", data_add_ncbi$BioSampleID[i])
    
    #make sure that the ID is for the BioSample database and not another
    if(length(bp2$ids) > 0){
      #extract the sample from NCBI:
      e_summary <- entrez_summary(db = "biosample", id = bp2$ids)
      
      #Deal with cases in which there are two entries for one sample,
      # second one seems to be the correct upon tests:
      if(length(e_summary) > 1 & length(e_summary) < 3){
        e_sub=e_summary[[2]]
        
        e_title = e_sub$title
        e_breed = e_sub$infraspecies
        e_species= e_sub$organism
      }else{
        #here go entries with a list length of 1
        e_breed = e_summary$infraspecies
        e_title = e_summary$title
        e_species= e_summary$organism
      }
      
      #here goes everything for lists that are not empty:
      #fill the elements with the information extracted before:
      if(length(e_title) >= 1){
        data_add_ncbi$ncbi_title[i] = e_title
      }
      if(length(e_breed) >= 1){
        data_add_ncbi$ncbi_infraspecies[i] = e_breed
      }
      if(length(e_species) >= 1){
        data_add_ncbi$ncbi_species[i] = e_species
      }
      
      #here goes how to deal with the empty parts of the NCBI entries:
      
      if(data_add_ncbi$ncbi_title[i] == FALSE){
        data_add_ncbi$ncbi_title[i] = NA
      }
      if(data_add_ncbi$ncbi_infraspecies[i] == FALSE){
        data_add_ncbi$ncbi_infraspecies[i] = NA
      }
      if(data_add_ncbi$ncbi_species[i] == FALSE){
        data_add_ncbi$ncbi_species[i] = NA
      }
      
      
    }else if(length(bp2$ids) == 0){
      data_add_ncbi$ncbi_infraspecies[i] = NA
      data_add_ncbi$ncbi_title[i] = NA
      data_add_ncbi$ncbi_species[i] = NA 
    }#closes the length > 0 for the ncbi entries
    
  }else if(is.na(data$BioSampleID[i])){
    data_add_ncbi$ncbi_infraspecies[i] = NA
    data_add_ncbi$ncbi_title[i] = NA
    data_add_ncbi$ncbi_species[i] = NA 
    
  } #closes the not-na logical 
  
} #closes the for-loop

end=Sys.time()


print(end-start)

write.table(x=data_add_ncbi,
            file="/Home/johanna/Dokumente/1000_bull_genomes/Data/Run9-TAUIND_add_NCBI_information.txt",
            sep="\t",
            row.names=F, 
            quote=F)
