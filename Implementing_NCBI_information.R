
library(ggplot2)
library(ggpubr)
library(rentrez)
library(openxlsx)

# Read the ncbi enriched meta-data: ------------
data=read.delim("/Home/johanna/Dokumente/1000_bull_genomes/Data/Run9-TAUIND_add_NCBI_information.txt",
                     header=T,
                     sep="\t")
data_plus=data
dim(data_plus) #correct number of rows

#examine the ncbi entries: -------
table(data_plus$ncbi_species)

#re-format variables
data_plus$Breed<-as.character(data_plus$Breed)
data_plus$ncbi_infraspecies<-as.character(data_plus$ncbi_infraspecies)
data_plus$ncbi_species<-as.character(data_plus$ncbi_species)
data_plus$ncbi_title<-as.character(data_plus$ncbi_title)




# Step 1 retrieval = extract breeds from infraspcies -----------

for(i in which(data_plus$Breed =="Unknown") ){
  
  if(length(grep("breed", data_plus$ncbi_infraspecies[i])) > 0){

      a = unlist(strsplit(data_plus$ncbi_infraspecies[i], split="breed: "))
      data_plus$Breed[i] = a[2]
      print(a[2])
    }
  }



# Step 2 retrieval= extract breeds from titles -------------------

for(i in which(data_plus$Breed =="Unknown")){
  
  if(length(grep("cattle", data_plus$ncbi_title[i])) > 0){
  
    a = unlist(strsplit(data_plus$ncbi_title[i], split=" "))
    
    if(length(a) == 2){
      data_plus$Breed[i] = a[1]
      #print(a[1])
      
    }else if(length(a) == 3 ){
    
    data_plus$Breed[i] = a[2]
    #print(a[2])
    }
  }   
}



# Step 3 retrieval: individuals with entries that cannot be automated: -------

#View(data_plus[which(data_plus$Breed =="Unknown"), c("ncbi_title", "ncbi_infraspecies", "ncbi_species")])


table(data[which(data$Breed =="Unknown"), "Species"])
table(data_plus[which(data_plus$Breed =="Unknown"), "Species"])

dim(as.table(table(as.character(data[which(data$Species =="Indicus"), "Breed"]))))
dim(as.table(table(as.character(data_plus[which(data_plus$Species =="Indicus"), "Breed"]))))


for(i in which(data_plus$Breed =="Unknown")){
  
  data_plus[grep("Gloucester", data_plus$ncbi_title), "Breed" ]= "Gloucester"
  data_plus[grep("Pustertaler Sprinzen", data_plus$ncbi_title), "Breed" ]= "Pustertaler Sprinzen"
  data_plus[grep("Irish moiled", data_plus$ncbi_title), "Breed" ]= "Irish moiled"
  data_plus[grep("British White", data_plus$ncbi_title), "Breed" ]= "British White"
  data_plus[grep("NDama", data_plus$ncbi_title), "Breed" ]= "NDama"
}



# Manuel corrections: make breed writings uniform --------------------

## BohaiBlack two individuals
data_plus[grep("Bohai Black", data_plus$ncbi_infraspecies),"Breed" ]= "BohaiBlack"



write.table(data_plus,
            "/Home/johanna/Dokumente/1000_bull_genomes/Data/Run9-TAUIND_recovery_by_NCBI_information.txt",
            col.names=T,
            row.names = ,
            quote = F,
           sep="\t")


# Complicated cases of hidden breed information in the NCBI entries: -----------

data<-read.delim("/Home/johanna/Dokumente/1000_bull_genomes/Data/Run9-TAUIND_recovery_by_NCBI_information.txt")
data$ncbi_title=as.character(data$ncbi_title)
data$Breed=as.character(data$Breed)

# left over samples:
sort(table(data[which(data$Breed == "Unknown"), "BioProjectID"]), decreasing=T)


##Bioproject PRJEB38981: 13 individuals -------------
data[grep("Asturiana",data$ncbi_title), "Breed" ]= "Asturiana"

##BioProject PRJNA283480: 75 individuals ----------
breed_codes<-read.xlsx("/Home/johanna/Dokumente/1000_bull_genomes/Data/MI_et_al_2018_Breed_codes/SuppTables .xlsx",
                       sheet=1,
                       colNames = T,
                       startRow = 2)


breed_codes=breed_codes[,c(1,2)]
breed_codes$Breed_name = rep("No_name", len=length(breed_codes$Breed))


for(i in 1:151){
  
  if(!is.na(breed_codes$Breed[i]) & is.na(breed_codes$Breed[i+1])){
    
    breed_codes$Breed[i+1] = breed_codes$Breed[i]
  }
  
  
  #recover the breed names with information from the publication of Mei et al.:
  if(breed_codes$Breed[i] == "QCC"){
    breed_codes$Breed_name[i] = "Qinchuan"
  }else if(breed_codes$Breed[i] == "NYC"){
    breed_codes$Breed_name[i] = "Nanyang"
  }else if(breed_codes$Breed[i] == "LXC"){
    breed_codes$Breed_name[i] = "Luxi"
  }else if(breed_codes$Breed[i] == "YBC"){
    breed_codes$Breed_name[i] = "Yanbian"
  }else if(breed_codes$Breed[i] == "YNC"){
    breed_codes$Breed_name[i] = "Yunnan"
  }else if(breed_codes$Breed[i] == "LQC"){
    breed_codes$Breed_name[i] = "Leiqiong"
  }else if(breed_codes$Breed[i] == "JBC"){
    breed_codes$Breed_name[i] = "Japanese black"
  }else if(breed_codes$Breed[i] == "RAN"){
    breed_codes$Breed_name[i] = "Red Angus"
  }
  
  
  # remove the notes on the usage of the samples:
  breed_codes$Sample = gsub("*#", "", breed_codes$Sample )
  breed_codes$Sample = gsub(" *", "", breed_codes$Sample )
  breed_codes$Sample = gsub("*", "", breed_codes$Sample )
}



mixed_samples=which(data$BioProjectID == "PRJNA283480")

#The following line is used in the controle for the function of the breed name entry
#collect=c()

for(i in 1:length(which(data$BioProjectID == "PRJNA283480"))){
  
  bp2=entrez_search(db="biosample", data[mixed_samples[i], "BioSampleID"])
  
  e_summary <- entrez_summary(db = "biosample", id = bp2$ids)
  
  sample=unlist(strsplit(unlist(strsplit(e_summary$identifiers, split=";"))[2], split="Sample name: BGI-"))[2]
  

  
  if(length(breed_codes[grep(sample, breed_codes$Sample), "Breed_name"])== 0){
    full_breed = breed_codes[grep(gsub("QC", "QX", sample), breed_codes$Sample), "Breed_name"]
  }else{
    full_breed = breed_codes[grep(sample, breed_codes$Sample), "Breed_name"]
  }
  
  #the breed names must not contain spaces:
  full_breed=gsub(" ", "",full_breed)
  
  #Fill in the Breeds in the meta-data:
  data[mixed_samples[i], "Breed"] = full_breed
  
}


write.table(data,
           "/Home/johanna/Dokumente/1000_bull_genomes/Data/Run9-TAUIND_recovery_by_NCBI_information_v2.txt",
            col.names=T,
            row.names = ,
            quote = F,
           sep="\t")



# Recovery IranAdmixed and Uganda Admixed: Breed grouping -------------
data<-read.delim("/Home/johanna/Dokumente/1000_bull_genomes/Data/Run9-TAUIND_recovery_by_NCBI_information_v2.txt")

data$ncbi_infraspecies=as.character(data$ncbi_infraspecies)
data$Breed=as.character(data$Breed)

data[grep("Rashoki", data$ncbi_infraspecies), "Breed" ]= "Rashoki"

## Determine number of breeds in the Uganda Admixed:
table(data[grep("UgandaAdmixed", data$Breed), "ncbi_infraspecies"])

Uganda_admixed=grep("UgandaAdmixed", data$Breed)
for(i in 1:length(Uganda_admixed)){
  if(length(grep("breed", data$ncbi_infraspecies[Uganda_admixed[i]])) > 0){
    a = unlist(strsplit(data$ncbi_infraspecies[Uganda_admixed[i]], split="breed: "))
    #print(gsub(" ", "",a[2]))
    data$Breed[Uganda_admixed[i]] = gsub(" ", "",a[2])
  }
}




## Japanese Native : breed grouping  ------
data[grep("Mishima", data$ncbi_title), "Breed" ]= "Mishima"
data[grep("Mishima", data$ncbi_title), "Species" ]= "Taurus"

data[grep("Kuchinoshima", data$ncbi_title), "Breed" ]= "Kuchinoshima"
data[grep("Kuchinoshima", data$ncbi_title), "Species" ]= "Taurus"




write.table(data,
  "/Home/johanna/Dokumente/1000_bull_genomes/Data/Run9-TAUIND_recovery_by_NCBI_information_v3.txt",
   col.names=T,
   row.names = ,
   quote = F,
   sep="\t")

