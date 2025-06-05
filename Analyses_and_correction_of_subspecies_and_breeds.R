
library(openxlsx)   # read the original meta-data
library(BGData)     # memory map the bed file
library(ggplot2)    # visulalization of results
library(viridis)    # color-blind friendly color palette for visuals
library(ggforce)    # for the facet zoom in PCA plot
library(reshape2)
library(ggpubr)
library(forcats)    #used for the re-ordering in the heatmaps
library(dplyr)
library(ggalluvial) # used for the visualization of the subspecies changes

#Load fam and corrected-meta files: ------------------------------------------------------------

fam_file<-read.table("/media/Scratch3/Johanna_1000Bull/ImputierteDateien/Chr1_BeaglePhased.fam")
colnames(fam_file)<-c("FID", "IID", "MAT", "PAT", "SEX", "Pheno")

#corrected Meta-File:
data<-read.delim("/Home/johanna/Dokumente/1000_bull_genomes/Data/Run9-TAUIND_recovery_by_NCBI_information_v2.txt",
                 header=T,
                 sep="\t",
                 stringsAsFactors = F)

#original meta file:
meta_original<-read.xlsx("/Home/johanna/Dokumente/1000_bull_genomes/Data/Run9-TAUINDraw-animalDistributionList-20210630.xlsx")

# For the visualisations we have to keep the order of the fam file and add the recovered Breeds from the NCBI:
fam_file_plus= merge( y=data[,c(1,8,9)], x=fam_file[,c(1:2)], by.x="IID", by.y = "SampleID")

#we want to re-assure that the fam-file order does not change:
fam_file_plus = fam_file_plus[fam_file$IID, ]
colnames(fam_file_plus)[1] = "SampleID"

all(fam_file_plus$SampleID == fam_file$IID)
all(fam_file_plus$FID == fam_file$FID)
all(fam_file_plus$Breed == fam_file$FID) #ok, this cannot be equal due to recoveries from the NCBI

# In several of the functions rownames have to be equal to the indices: 
rownames(fam_file_plus) = 1:nrow(fam_file_plus)


#Create PCA file: ------------------------------------------------------------------------------

genome<-data.frame(matrix(nrow=29, ncol=3))
colnames(genome)<-c("Chromosome", "Variants", "Animals")

single_chroms<-list() #create the list objects to store the vitual matrices in

#start<-Sys.time()
setwd("/media/Scratch3/Johanna_1000Bull/ImputierteDateien/")
for(i in 1:29){
  
  single_chroms[[i]] <- BEDMatrix(paste0("Chr",i,"_BeaglePhased.bed"), simple_names = TRUE) 
  
  genome$Chromosome[i]<-i
  genome$Variants[i]<-ncol(single_chroms[[i]])
  genome$Animals[i]<-nrow(single_chroms[[i]])
  
  cat("Chromosome_",i,"_done", "\n")
}

comb_chroms<- as.ColumnLinkedMatrix(single_chroms)
rm(single_chroms)


#define the sample size of variants:
sample_size=50000

#sample variants across all autosomes
random_sample=sample(1:ncol(comb_chroms), size=sample_size)
length(random_sample)

pca_allIndiv_50k=prcomp(comb_chroms[,random_sample])

#save the results for later, to base the pictures on the same pca
#saveRDS(pca_allIndiv_50k, "/Home/johanna/Dokumente/1000_bull_genomes/Data/2025-March-21_PCA_all_individuals_50k.RData")

#clear the workspace from files no longer important
#rm(comb_chroms, genome , pca_allIndiv_50k)

#Load the PCA in later usages:
pca_allIndiv=readRDS("/Home/johanna/Dokumente/1000_bull_genomes/Data/2025-March-21_PCA_all_individuals_50k.RData")

# Load the genomic distances files and process the names: -------------------------------------
sim_by_ibs=read.delim("/media/Scratch3/Johanna_1000Bull/Genomic_distances/GRM_by_1-IBS_all_individuals.mdist")
sim_by_ibs_names=read.table("/media/Scratch3/Johanna_1000Bull/Genomic_distances/GRM_by_1-IBS_all_individuals.mdist.id")

colnames(sim_by_ibs)=sim_by_ibs_names$V1

rownames(sim_by_ibs)=sim_by_ibs_names$V1[-1] #first is eliminated since there is no diagonal


# Large PCA for Bohai Black individuals : ---------------------------------------------------------------
fam_sub=fam_file_plus
breed1="BohaiBlack"
breed2="Wenshan"

individuals <- fam_file_plus$SampleID[which(fam_file_plus$Breed == breed1 | fam_file_plus$Breed == breed2)]


target_rows<-which(fam_file_plus$SampleID %in% individuals)

fam_sub[which(!(rownames(fam_file_plus) %in% target_rows)), "Breed"] <-"others"

df<-data.frame(Species=fam_sub$Breed, x=pca_allIndiv$x[,1], y=pca_allIndiv$x[,2])
df$Species<-as.character(df$Species)


df<-data.frame(Breed=fam_sub$Breed,  x=pca_allIndiv$x[,1], y=pca_allIndiv$x[,2])
breed3_sub<-subset(df,  Breed=="Wenshan" )
breed4_sub<-subset(df,  Breed=="BohaiBlack")


pca_variances<-pca_allIndiv$sdev^2 / sum(pca_allIndiv$sdev^2)*100


pca_plot_b1<-ggplot(df, aes(x=x, y=y)) +
  geom_point(aes(color= Breed) )+
  
  scale_colour_manual("Breed", breaks=c(  "Wenshan", "BohaiBlack", "others"),
                      values = c("red","blue",  "grey")) +
  
  theme(legend.position = "none") +
  geom_point( data=breed3_sub, aes(x=x, y=y), colour="red")+ #Wenshan
  geom_point( data=breed4_sub, aes(x=x, y=y), colour="blue")+ #Bohai
  labs(x=paste0("<- taurine                          ",
                "PC1  ", round(pca_variances[1], digits=2), "%",
                "                         indicine ->"),
       y=paste0("PC2  ",round(pca_variances[2], digits=2), "%"))



ggsave(file="/Home/johanna/Dokumente/1000_bull_genomes/Graphics/For_publication/PCA_species_for_paper_Bohai_Wenshan.png",
       plot=pca_plot_b1,
       dpi=300, 
       width=3000,
       height=1500, unit="px")




# Large PCA plot for original sub-species of all individuals: --------------------------------------------

# For the visualisations we have to keep the order of the fam file and add the recovered Breeds from the NCBI:
fam_original = merge( y = meta_original[,c(1,8,9)], x=fam_file[,c(1:2)], by.x="IID", by.y = "SampleID")

fam_original = fam_original[fam_file$IID, ]
colnames(fam_original)[1] = "SampleID"

all(fam_original$SampleID == fam_file$IID)
all(as.character(fam_original$FID) == as.character(fam_file$FID))
rm(fam_file)

# In several of the functions rownames have to be equal to the indices: 
rownames(fam_original) = 1:nrow(fam_original)

fam_sub = fam_original

individuals <- fam_original$SampleID[which(fam_original$Species == "Taurus" |fam_original$Species == "Indicus" )]
length(individuals)

target_rows<-which(fam_original$SampleID %in% individuals)

fam_sub[which(!(rownames(fam_original) %in% target_rows)), "Species"] <-"others"

df<-data.frame(Species=fam_sub$Species, x=pca_allIndiv$x[,1], y=pca_allIndiv$x[,2])
df$Species<-as.character(df$Species)


pca_variances<-pca_allIndiv$sdev^2 / sum(pca_allIndiv$sdev^2)*100

pca_species<-ggplot(data=df,aes(x=x, y=y, colour=Species)) + 
  geom_point() +
  scale_color_manual(values=c("others"="grey", "Indicus" = "blue", "Taurus" = "red"))+
  theme(legend.position = "none") +
  labs(x=paste0("<- taurine                          ",
                "PC1  ", round(pca_variances[1], digits=2), "%",
                "                         indicine ->"),
       y=paste0("PC2  ",round(pca_variances[2], digits=2), "%"))




figure_plots<-ggarrange(pca_species,                                   
                        plot,    #from the SCOPE script
                        ncol = 1, 
                        labels = c("A", "B"),
                        nrow = 2) 

ggsave(file="/Home/johanna/Dokumente/1000_bull_genomes/Graphics/For_publication/PCA_species_for_paper_Mar_23_2025.png",
       plot=figure_plots,
       dpi=300, 
       width=2500,
       height=2000, unit="px")



## Mass-application for all breeds or species questionable: ---------------------------------

pca_function = function(target_breed = "no", target_species ="no", legend_option="Species"){
  
  df<-data.frame(Breed=fam_file_plus$Breed, Species=fam_file_plus$Species, x=pca_allIndiv$x[,1], y=pca_allIndiv$x[,2])
  
  df$Breed<-as.character(df$Breed)
  
  
  pca_variances<-pca_allIndiv$sdev^2 / sum(pca_allIndiv$sdev^2)*100
  
  
  if(target_breed != "no" & target_species == "no"){
    df2 <- df %>%
      arrange(Breed == target_breed)  # FALSE (target_breed) comes after TRUE (others)
    
    df2$Target = if_else(df2$Breed == target_breed, T,F)
    
    plot=ggplot(df2, aes(x = x, y = y, color = Target)) +
      geom_point() +
      scale_color_manual(
        values = c("FALSE" = "gray", "TRUE" = "red" ),   # put the target second to place it ontop of the other points
        labels = c("FALSE" = "Others", "TRUE" = target_breed),
        name = "Breed")+
      theme(legend.position = "top",
            legend.title = element_text(size=16),
            legend.text = element_text(size=14))+
      #average quantile1 and qunatile3 positions from the scope analysis
      geom_vline(xintercept = c(10.03, 53.74), linetype="dotted", size = 0.4)+
      annotate("text",
               x = 10.03 , y = 20,
               label = "Q1",
               color = "black")+
      annotate("text",
               x = 53.74, y = 20,
               label = "Q3", color = "black")+
      labs(x=paste0("PC1  ", round(pca_variances[1], digits=2), "%"),
           y=paste0("PC2  ",round(pca_variances[2], digits=2), "%"))
    
  }else if(target_breed=="no" & target_species!="no"){
    df2 <- df %>%
      arrange(Species == target_species)  # FALSE (target_breed) comes after TRUE (others)
    
    df2$Target = if_else(df2$Species == target_species, T,F)
    
    plot=ggplot(df2, aes(x = x, y = y, color = Target)) +
      geom_point() +
      scale_color_manual(
        values = c("FALSE" = "gray", "TRUE" = "red" ),   # put the target second to place it ontop of the other points
        labels = c("FALSE" = "Others", "TRUE" = target_species),
        name = "Species")+
      theme(legend.position = "top",
            legend.title = element_text(size=16),
            legend.text = element_text(size=14)) +
      #average quantile1 and qunatile3 positions from the scope analysis
      geom_vline(xintercept = c(10.03, 53.74), linetype="dotted", size = 0.4)+
      annotate("text",
               x = 10.03 , y = 20,
               label = "Q1",
               color = "black")+
      annotate("text",
               x = 53.74, y = 20,
               label = "Q3", color = "black")+
      labs(x=paste0("PC1  ", round(pca_variances[1], digits=2), "%"),
           y=paste0("PC2  ",round(pca_variances[2], digits=2), "%"))
  }else if(target_breed!="no" & target_species!="no"){
    if(legend_option=="Species"){
      df2 <- df %>%
        arrange(Species == target_species)  # FALSE (target_breed) comes after TRUE (others)
      
      df2$Target = if_else(df2$Species == target_species & df2$Breed == target_breed, T,F)
      
      plot=ggplot(df2, aes(x = x, y = y, color = Target)) +
        geom_point() +
        scale_color_manual(
          values = c("FALSE" = "gray", "TRUE" = "red" ),   # put the target second to place it ontop of the other points
          labels = c("FALSE" = "Others", "TRUE" = target_species),
          name = "Species")+
        theme(legend.position = "top",
              legend.title = element_text(size=16),
              legend.text = element_text(size=14)) +
        #average quantile1 and qunatile3 positions from the scope analysis
        geom_vline(xintercept = c(10.03, 53.74), linetype="dotted", size = 0.4)+
        annotate("text",
                 x = 10.03 , y = 20,
                 label = "Q1",
                 color = "black")+
        annotate("text",
                 x = 53.74, y = 20,
                 label = "Q3", color = "black")+
        labs(x=paste0("PC1  ", round(pca_variances[1], digits=2), "%"),
             y=paste0("PC2  ",round(pca_variances[2], digits=2), "%"))
    }else if(legend_option=="Breed"){
      df2 <- df %>%
        arrange(Species == target_species)  # FALSE (target_breed) comes after TRUE (others)
      
      df2$Target = if_else(df2$Species == target_species & df2$Breed == target_breed, T,F)
      
      plot=ggplot(df2, aes(x = x, y = y, color = Target)) +
        geom_point() +
        scale_color_manual(
          values = c("FALSE" = "gray", "TRUE" = "red" ),   # put the target second to place it ontop of the other points
          labels = c("FALSE" = "Others", "TRUE" = target_breed),
          name = "Breed")+
        theme(legend.position = "top",
              legend.title = element_text(size=16),
              legend.text = element_text(size=14)) +
        #average quantile1 and qunatile3 positions from the scope analysis
        geom_vline(xintercept = c(10.03, 53.74), linetype="dotted", size = 0.4)+
        annotate("text",
                 x = 10.03 , y = 20,
                 label = "Q1",
                 color = "black")+
        annotate("text",
                 x = 53.74, y = 20,
                 label = "Q3", color = "black")+
        labs(x=paste0("PC1  ", round(pca_variances[1], digits=2), "%"),
             y=paste0("PC2  ",round(pca_variances[2], digits=2), "%"))
    }
  }
  return(plot)
}



# Breeds with uncertain or wrong species assignments PCA: -------------------

breeds_in_question = c("Ankole", "Benishangul", "BohaiBlack","Crossbreed", "Butana", "Dehong", 
                       "Fujian", "Gourounsi", "JapaneseNative","JiaxianRed",
                       "Jian", "Jinjiang", "Kenana", "Leiqiong", "Lingnan", "Maure",
                       "Nanyang","NariMaster", "Nuer","Qinchuan","SichuanIndigenous", "Weining", "Wenshan")

length(breeds_in_question)

# Generate PCA plots 
pca_collection=list() 

for(i in 1:length(breeds_in_question)){
  if(breeds_in_question[i]=="Crossbreed"){
    pca_collection[[i]] = pca_function(target_breed = breeds_in_question[i],
                                       target_species = "other",
                                       legend_option = "Breed")
  }else{
    pca_collection[[i]] = pca_function(target_breed = breeds_in_question[i])
  }
}

# Save the plots
Breed_plots=ggarrange(plotlist=pca_collection[1:12],
                      ncol=3,
                      nrow=4,
                      labels = LETTERS[1:12])

Breed_plots2 = ggarrange(plotlist=pca_collection[13:23],
                         ncol=3,
                         nrow=4,
                         labels = LETTERS[1:11])

ggsave(file="/Home/johanna/Dokumente/1000_bull_genomes/Graphics/For_publication/PCA_Breeds_nop_species_for_paper_1.png",
       plot=Breed_plots,
       dpi=300, 
       width=5000,
       height=4000,
       unit="px")

ggsave(file="/Home/johanna/Dokumente/1000_bull_genomes/Graphics/For_publication/PCA_Breeds_nop_species_for_paper_2.png",
       plot=Breed_plots2,
       dpi=300, 
       width=5000,
       height=4000,
       unit="px")

# Groups with missing subspecies : ---------------------

species_in_question = c("other", "unknown", "PredictedTaurus", "ancient")
breed_in_question = c("Unknown", "Unknown", "Unknown", "Ancient")

length(species_in_question)

# Generate PCA plots 
pca_collection=list() 

for(i in 1:length(species_in_question)){
  pca_collection[[i]] = pca_function(target_species =  species_in_question[i], target_breed = breed_in_question[i])
}

# Save the plots
species_plots=ggarrange(plotlist=pca_collection[1:length(species_in_question)],
                        ncol=2,
                        nrow=2,
                        labels = LETTERS[1:length(species_in_question)])

ggsave(file="/Home/johanna/Dokumente/1000_bull_genomes/Graphics/For_publication/PCA_Breeds_no_species_approximation.png",
       plot=species_plots,
       dpi=300, 
       width=2500,
       height=2000,
       unit="px")



# Function for admixture plotting from Felix: -------------
preparePlotData_SCOPE = function(famDF, admixtureDF){
  
  breed_counts=table(famDF$Breed)
  bc=c()
  for(i in 1:length(breed_counts)){
    bc=c(bc,rep(breed_counts[i], len=breed_counts[i]*2))
  }
  
  combinedAdmixtureDF = rbind(as.character(famDF$Breed), admixtureDF)
  plotdf_admixture = data.frame(stringsAsFactors = F)
  for(i in 1:ncol(combinedAdmixtureDF)){
    breed = combinedAdmixtureDF[1,i]
    for(j in 2:nrow(combinedAdmixtureDF)){
      plotdf_admixture = rbind(plotdf_admixture, data.frame(i, breed, (j-1), as.numeric(combinedAdmixtureDF[j,i])))
    }
  }
  colnames(plotdf_admixture) = c("ID","Breed", "Subspecies", "Percent")
  plotdf_admixture$Subspecies = as.factor(plotdf_admixture$Subspecies)
  plotdf_admixture$ID = as.factor(plotdf_admixture$ID)
  plotdf_admixture$Percent = round(plotdf_admixture$Percent, digits = 2)
  
  plotdf_admixture= plotdf_admixture[order(as.character(plotdf_admixture$Breed)),]
  plotdf_admixture$Breed = paste0(as.character(plotdf_admixture$Breed), ", ", "n=", bc)
  
  
  return(plotdf_admixture)
}


createAdmixturePlot = function(plotDF){
  p =ggplot(data = plotDF, aes(x = ID, y = Percent, fill = Subspecies)) +
    geom_bar(stat = "identity") +
    scale_fill_brewer(palette = "Paired", labels=c("B.indicus","B.taurus")) +
    facet_wrap(~ Breed, scales = "free_x")+
    ylim(0,1)+
    theme(
      axis.text.x = element_blank(),  # Remove x-axis text labels
      axis.title.x = element_blank(), # Remove x-axis title
      panel.grid = element_blank(),
      axis.ticks.x = element_blank()
    )
  return(p)
}




# Breeds with uncertain or wrong species assignments SCOPE: -------------------

breeds_in_question_scope =  c("Ankole", "Benishangul", "BohaiBlack", "Butana", "Dehong", 
                              "Fujian", "Gourounsi", "JapaneseNative","JiaxianRed", "Jian", "Jinjiang", "Kenana", 
                              "Leiqiong", "Lingnan", "Maure", "Nanyang","NariMaster", "Nuer","Qinchuan",
                              "SichuanIndigenous", "Weining", "Wenshan")

#There is one crossbreed individual that has no subspecies assignment, add it by the BioSample ID
odd_individual_crossbreed = which(data$BioSampleID == "SAMN01915343")


df<-data.frame(ID=fam_file_plus$SampleID, x=pca_allIndiv$x[,1], y=pca_allIndiv$x[,2])

#we want to compare the breeds without assignments to the polar ends of the PC1:
taurine = which(df$x < -3.009293 ) #section1 end
indicine =which(df$x >= 68.990707  ) #final section start

s_unclear = c(which(fam_file_plus$Breed %in% breeds_in_question_scope), odd_individual_crossbreed)
length(s_unclear)

#We want to ensure that we do not sample duplicates
s_taurine = sample(taurine[!(taurine %in% s_unclear)], size=100, replace=F)
s_indicine = sample(indicine[!(indicine %in% s_unclear)], size=100, replace=F)


sample_ids = c(s_taurine, s_indicine, s_unclear)
length(sample_ids) #cross-check the number of IDs

individuals = as.character(fam_file_plus$SampleID[sample_ids])
length(unique(individuals)) #cross-check the number of individuals

#PLINK requires a two-column file for --keep, create this by duplicating the IDs
sample_fam <- data.frame(individuals, individuals) 

#write.table(sample_fam, "/media/Scratch3/Johanna_1000Bull/ScopeAnalysis/Breeds_unclear_species/Filter_file_PCA_unclear_species.txt",
#          quote=F,col.names=F, row.names=F, sep="\t") 




# Pre-process the scope Outputs to fit th PCA order:
setwd("/media/Scratch3/Johanna_1000Bull/")
ffam_o2 = read.table("/media/Scratch3/Johanna_1000Bull/ScopeAnalysis/Breeds_unclear_species/Chr_Merged_BeaglePhased_PCA_unclear_species.fam")
admixtureDF= read.table("/media/Scratch3/Johanna_1000Bull/ScopeAnalysis/Breeds_unclear_species/Chr_Merged_BeaglePhased_PCA_unclear_species_USV_k2Qhat.txt")

df<-data.frame(IID=fam_file_plus$SampleID, x=pca_allIndiv$x[,1], y=pca_allIndiv$x[,2], meta= fam_file_plus$Species)

df_to_order = df[which(df$IID %in% ffam_o2$V1), ] #subset the PCA df


## add the x coordinate to the Fam-file:
ffam_with_x = merge( ffam_o2, df_to_order[, c("IID", "x")], by.y="IID", by.x="V1")

all(as.character(ffam_with_x$v2)==as.character(ffam_o2$V1)) #order is same

ffam_with_x$V1 = fam_file_plus[which(fam_file_plus$SampleID %in% as.character(ffam_with_x$V2)), "Breed"]

all(as.character(ffam_with_x$IID) == as.character(ffam_o2$V1)) #order is same


#Subset to breeds:
famDF = ffam_with_x[which(ffam_with_x$V1 == "Ankole" |
                            ffam_with_x$V1 == "Benishangul" |
                            ffam_with_x$V1 == "BohaiBlack" | 
                            ffam_with_x$V1 == "Butana" |
                            ffam_with_x$V1 == "Crossbreed" |
                            ffam_with_x$V1 == "Dehong" |
                            ffam_with_x$V1 == "Fujian" |
                            ffam_with_x$V1 == "Gourounsi" |
                            ffam_with_x$V1 == "JapaneseNative" |
                            ffam_with_x$V1 == "JiaxianRed" |
                            ffam_with_x$V1 == "Jian" |
                            ffam_with_x$V1 == "Jinjiang" |
                            ffam_with_x$V1 == "Kenana" |
                            ffam_with_x$V1 == "Kuchinoshima" | 
                            ffam_with_x$V1 == "Leiqiong" |
                            ffam_with_x$V1 == "Lingnan" |
                            ffam_with_x$V1 == "Maure" |
                            ffam_with_x$V1 == "Nanyang" |
                            ffam_with_x$V1 == "NariMaster" |
                            ffam_with_x$V1 == "Nuer" |
                            ffam_with_x$V1 == "Qinchuan" |
                            ffam_with_x$V1 == "SichuanIndigenous" |
                            ffam_with_x$V1 == "Weining" |
                            ffam_with_x$V1 == "Wenshan"),]



colnames(famDF)[1]<-"Breed"

#subset the SCOPE results
admixtureDF=admixtureDF[  ,which(ffam_with_x$V1 == "Ankole" |
                                   ffam_with_x$V1 == "Benishangul" |
                                   ffam_with_x$V1 == "BohaiBlack" | 
                                   ffam_with_x$V1 == "Butana" |
                                   ffam_with_x$V1 == "Crossbreed" |
                                   ffam_with_x$V1 == "Dehong" |
                                   ffam_with_x$V1 == "Fujian" |
                                   ffam_with_x$V1 == "Gourounsi" |
                                   ffam_with_x$V1 == "JapaneseNative" |
                                   ffam_with_x$V1 == "JiaxianRed" |
                                   ffam_with_x$V1 == "Jian" |
                                   ffam_with_x$V1 == "Jinjiang" |
                                   ffam_with_x$V1 == "Kenana" |
                                   ffam_with_x$V1 == "Kuchinoshima" | 
                                   ffam_with_x$V1 == "Leiqiong" |
                                   ffam_with_x$V1 == "Lingnan" |
                                   ffam_with_x$V1 == "Maure" |
                                   ffam_with_x$V1 == "Nanyang" |
                                   ffam_with_x$V1 == "NariMaster" |
                                   ffam_with_x$V1 == "Nuer" |
                                   ffam_with_x$V1 == "Qinchuan" |
                                   ffam_with_x$V1 == "SichuanIndigenous" |
                                   ffam_with_x$V1 == "Weining" |
                                   ffam_with_x$V1 == "Wenshan")]


plotdf_admixture_species=preparePlotData_SCOPE(famDF, admixtureDF)

adm_plot_species <- createAdmixturePlot(plotdf_admixture_species)


ggsave(plot=adm_plot_species,
       filename= "/Home/johanna/Dokumente/1000_bull_genomes/Graphics/For_publication/For_Publication_ADM_PCA_unclear_species.png",
       dpi=600,
       units="px", 
       height=3000,
       width= 6000
)


## Percent taurine and indicine per breeds on average for table in publication:
breeds_in_question_scope = c(breeds_in_question_scope, "Crossbreed")
ancestry=data.frame(matrix(ncol=3, nrow=length(breeds_in_question_scope)))

colnames(ancestry) = c("Breed", "Taurus", "Indicus")

for(i in 1:length(breeds_in_question_scope)){
  avg=admixtureDF[,which(famDF$Breed == breeds_in_question_scope[i])]
  avg=t(avg)
  
  ancestry[i,1] = breeds_in_question_scope[i]
  ancestry[i,2] = round(mean(avg[,2]), digits=2) #taurine
  ancestry[i,3] = round(mean(avg[,1]), digits=2) #indicine
  
}


for(i in 1:nrow(ancestry)){
  cat(ancestry[i,1], "&", ancestry[i,2], "&", ancestry[i,3], "\n")
}

ancestry=ancestry[order(ancestry$Indicus, decreasing=F),]
rownames(ancestry)=1:nrow(ancestry)
ancestry[,c(1,3)]



# Breed misassignments: ------------------------------------------------------------------------

## Function PCA graphics for two or three breeds: --------------------
pca_function_Groups = function(breed1="no", breed2="no", breed3="no", species="no", legend=TRUE) {
  
  fam_sub = fam_file_plus
  #Deteime which individuals belong to the breeds investigated
  if (breed3 == "no") {
    individuals <- fam_file_plus$SampleID[fam_file_plus$Breed %in% c(breed1, breed2)]
    breed_list <- c("others", breed1, breed2)
    color_list <- c("grey", "red", "blue")
  } else {
    individuals <- fam_file_plus$SampleID[fam_file_plus$Breed %in% c(breed1, breed2, breed3)]
    breed_list <- c("others", breed1, breed2, breed3)
    color_list <- c("grey", "red", "blue", "orange")
  }
  
  # Identify rows corresponding to selected individuals
  target_rows <- which(fam_file_plus$SampleID %in% individuals)
  
  # Later we want all other individulas to appear in grey, so form one group of "others" here
  fam_sub[!(rownames(fam_file_plus) %in% target_rows), "Breed"] <- "others"
  
  
  df <- data.frame(Breed = fam_sub$Breed, x = pca_allIndiv$x[,1], y = pca_allIndiv$x[,2])
  
  # Later we want to ensure that the colored points end up on top of the grey
  breed1_sub <- subset(df, Breed == breed1)
  breed2_sub <- subset(df, Breed == breed2)
  if (breed3 != "no") {
    breed3_sub <- subset(df, Breed == breed3)
  }
  
  # For the axis labeling
  pca_variances <- (pca_allIndiv$sdev^2 / sum(pca_allIndiv$sdev^2)) * 100
  
  # PCA plot basis to layer individulas on top
  pca_plot <- ggplot(df, aes(x = x, y = y, color = Breed)) +
    geom_point() +
    scale_colour_manual("Breed", breaks = breed_list, values = color_list)
  
  if (breed3 != "no") {
    pca_plot <- pca_plot + 
      geom_point(data = breed3_sub, aes(x = x, y = y), colour = "orange")+
      geom_point(data = breed1_sub, aes(x = x, y = y), colour = "red") +
      geom_point(data = breed2_sub, aes(x = x, y = y), colour = "blue")
  }else{
    pca_plot <- pca_plot +
      geom_point(data = breed1_sub, aes(x = x, y = y), colour = "red") +
      geom_point(data = breed2_sub, aes(x = x, y = y), colour = "blue")
    
  }
  
  # Add general labs 
  pca_plot <- pca_plot +
    labs(x = paste0("PC1  ", round(pca_variances[1], 2), "%"),
         y = paste0("PC2  ", round(pca_variances[2], 2), "%"))
  
  #Determine weather to use a legend or not
  if(legend==FALSE){
    pca_plot=pca_plot + theme(legend.position = "none") 
  }else if (legend==TRUE){
    pca_plot=pca_plot + theme(legend.position = "top") 
  }
  
  #Write out the plot object:
  return(pca_plot)
}


## Function for the the SCOPE Plots for the breeds, ammended by samples size: ----------------------------------------
preparePlotData_SCOPE = function(famDF, admixtureDF, n){
  combinedAdmixtureDF = rbind(as.character(famDF$Breed), admixtureDF)
  plotdf_admixture = data.frame(stringsAsFactors = F)
  for(i in 1:ncol(combinedAdmixtureDF)){
    breed = combinedAdmixtureDF[1,i]
    for(j in 2:nrow(combinedAdmixtureDF)){
      plotdf_admixture = rbind(plotdf_admixture, data.frame(i, breed, (j-1), as.numeric(combinedAdmixtureDF[j,i])))
    }
  }
  colnames(plotdf_admixture) = c("ID","Breed", "Pop", "Percent")
  plotdf_admixture$Pop = as.factor(plotdf_admixture$Pop)
  plotdf_admixture$ID = as.factor(plotdf_admixture$ID)
  plotdf_admixture$Percent = round(plotdf_admixture$Percent, digits = 2)
  
  
  for(n in 1:n){
    diff = 1.00 - sum(plotdf_admixture$Percent[plotdf_admixture$ID == n])
    plotdf_admixture$Percent[plotdf_admixture$ID == n][which.max(plotdf_admixture$Percent[plotdf_admixture$ID == n])] = max(plotdf_admixture$Percent[plotdf_admixture$ID == n]) + diff
  }
  
  
  return(plotdf_admixture)
}

createAdmixturePlot = function(plotDF){
  p =ggplot(data = plotDF, aes(x = ID, y = Percent, fill = Pop)) +
    geom_bar(stat = "identity") +
    scale_fill_brewer(palette = "Paired") +
    facet_wrap(~ Breed, scales = "free_x")+
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),  # Remove x-axis text labels
      axis.title.x = element_blank(), # Remove x-axis title
      panel.grid = element_blank(),
      strip.text = element_text(size = 10)
    )+
    scale_y_continuous(breaks=c(0.0, 0.25, 0.5, 0.75, 1.0), labels = c( 0.0, 25, 50, 75, 100))
  
  return(p)
}



## Function writing the IBS plot alone: -----------

GroupIBSPlot <- function(breedA, breedB, breedC="no", fam_file_plus, sim_by_ibs ) {
  
  # Identify individuals belonging to the two breeds
  individuals_breedA = fam_file_plus$SampleID[which(fam_file_plus$Breed == breedA)]
  individuals_breedB = fam_file_plus$SampleID[which(fam_file_plus$Breed == breedB)]
  if(breedC != "no"){
    individuals_breedC = fam_file_plus$SampleID[which(fam_file_plus$Breed == breedC)]
  }
  
  
  breedA_samples = which(colnames(sim_by_ibs) %in% individuals_breedA)
  breedB_samples = which(colnames(sim_by_ibs) %in% individuals_breedB)
  if(breedC != "no"){
    breedC_samples = which(colnames(sim_by_ibs) %in% individuals_breedC)
    individuals = c(breedA_samples, breedB_samples, breedC_samples)
  }else{
    individuals = c(breedA_samples, breedB_samples)
  }
  
  # Create IBS similarity matrix
  sim_matrix = as.matrix(sim_by_ibs[individuals-1, individuals])
  
  
  
  # Rename column names with generic identifiers
  newnames_breedA = paste0(rep(paste0(breedA, "-"), times=length(breedA_samples)), seq(1:length(breedA_samples)))
  newnames_breedB = paste0(rep(paste0(breedB, "-"), times=length(breedB_samples)), seq(1:length(breedB_samples)))
  colnames(sim_matrix)[which(colnames(sim_matrix) %in% individuals_breedA)] = newnames_breedA
  colnames(sim_matrix)[which(colnames(sim_matrix) %in% individuals_breedB)] = newnames_breedB
  if(breedC != "no"){
    newnames_breedC = paste0(rep(paste0(breedC, "-"), times=length(breedC_samples)), seq(1:length(breedC_samples)))
    colnames(sim_matrix)[which(colnames(sim_matrix) %in% individuals_breedC)] = newnames_breedC
  }
  
  
  # Convert to long format for plotting
  data1 = melt(sim_matrix)
  data1$Var1 = as.factor(data1$Var1)
  data1$Var2 = as.factor(data1$Var2)
  
  #Re-orede the dataframe to obtain the nice lower-triangle layout
  data1=data1 %>%
    mutate(name = fct_reorder(Var1, Var2, .fun="length"))
  
  if(breedC != "no"){
    # Heatmap plot
    heatmap_plot = ggplot(data1, aes(x = Var1, y = Var2, fill = value)) +
      geom_tile() + 
      scale_fill_viridis(discrete = FALSE, limits = c(0, 0.11), name = "Distance", direction= -1) +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),  
        axis.title.x = element_blank(), 
        axis.text.y = element_text(vjust = 1, hjust = 1, size = 10),
        axis.title.y= element_blank(),
        panel.grid = element_blank()
      ) +
      scale_y_discrete(breaks = c(newnames_breedA[1],newnames_breedB[1], newnames_breedC[1]))
  }else{
    # Heatmap plot
    heatmap_plot = ggplot(data1, aes(x = Var1, y = Var2, fill = value)) +
      geom_tile() + 
      scale_fill_viridis(discrete = FALSE, limits = c(0, 0.11), name = "Distance", direction= -1) +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),  
        axis.title.x = element_blank(), 
        axis.text.y = element_text(vjust = 1, hjust = 1, size = 10),
        axis.title.y= element_blank(),
        panel.grid = element_blank()
      ) +
      scale_y_discrete(breaks = c(newnames_breedA[1],newnames_breedB[1]))
  }
  
  
  return(heatmap_plot)
}




## Function calling all plots and arrange them: ------------------------------------------------

createBreedComparisonPlot <- function(breedA, breedB, breedC="no", fam_file_plus, sim_by_ibs, famDF, admixtureDF_k3, admixtureDF_k2, n_individuals, legend_decision) {
  
  # PCA Analysis
  PCA_plot = pca_function_Groups(breed1 = breedA, breed2 = breedB, breed3 = breedC, legend=legend_decision)+
    theme(plot.margin = margin(10, 15, 5, 20))
  
  #IBS Plots:
  heatmap_plot = GroupIBSPlot(breedA = breedA, breedB = breedB, breedC = breedC,
                              fam_file_plus = fam_file_plus,
                              sim_by_ibs = sim_by_ibs)+
    theme(plot.margin = margin(10, 15, 5, 20))  # all in pt
  
  # Generate admixture plots
  plotdf_admixture_k3 = preparePlotData_SCOPE(famDF, admixtureDF_k3, n = n_individuals)
  adm_plot_k3 <- createAdmixturePlot(plotdf_admixture_k3)+
    theme(plot.margin = margin(10, 15, 5, 20))  # all in pt
  
  plotdf_admixture_k2 = preparePlotData_SCOPE(famDF, admixtureDF_k2, n = n_individuals)
  adm_plot_k2 <- createAdmixturePlot(plotdf_admixture_k2)+
    theme(plot.margin = margin(10, 15, 5, 20))  # all in pt
  
  # Arrange all plots into a final figure
  final_plot = ggarrange(PCA_plot, heatmap_plot, adm_plot_k2, adm_plot_k3,
                         nrow = 2, ncol = 2,
                         labels = c("A", "B", "C", "D"),
                         heights = c(1, 1.5))
  
  
  
  return(final_plot)
}




## Function calling all plots and arrange them: ------------------------------------------------

createBreedComparisonPlot_small <- function(breedA, breedB, breedC="no", fam_file_plus, sim_by_ibs, lettering=TRUE, legend_decision, label1, label2) {
  
  # PCA Analysis
  PCA_plot = pca_function_Groups(breed1 = breedA, breed2 = breedB, breed3 = breedC, legend=legend_decision)
  
  # IBS plot
  heatmap_plot = GroupIBSPlot(breedA, breedB, breedC, fam_file_plus, sim_by_ibs)
  
  # Arrange all plots into a final figure, with or without letters
  if (lettering == TRUE){
    final_plot = ggarrange(PCA_plot, heatmap_plot,
                           nrow = 1, ncol = 2,
                           labels = c(label1, label2))
  }else if (lettering == FALSE){
    final_plot = ggarrange(PCA_plot, heatmap_plot,
                           nrow = 1, ncol = 2)
  }
  
  return(final_plot)
}


# Genrate the plot arrangements for breeds arising quastions: ------------------------------------

### FinnishAyrshire ----------
setwd("/Home/johanna/Dokumente/1000_bull_genomes/SCopr_ancestral/Ayrshire/")
famDF = read.table("Chr_Merged_BeaglePhased_Ayrshire.fam")
colnames(famDF) = c("Breed", "ID","V3","V4","V5", "V6")
admixtureDF2= read.table("Chr_Merged_BeaglePhased_Ayrshire_k2_USVQhat.txt")
admixtureDF3= read.table("Chr_Merged_BeaglePhased_Ayrshire_k3_USVQhat.txt")


final_graphic_FinAyr <- createBreedComparisonPlot(
  breedB = "FinnishAyrshire",
  breedA = "AyrshireFinnish",
  fam_file_plus = fam_file_plus,
  sim_by_ibs = sim_by_ibs,
  famDF = famDF,
  admixtureDF_k3 = admixtureDF3,
  admixtureDF_k2 = admixtureDF2,
  n_individuals = 57,
  legend_decision = FALSE
)


ggsave(plot=final_graphic_FinAyr,
       filename= "/Home/johanna/Dokumente/1000_bull_genomes/Graphics/For_publication/For_Publication_FinAyr.png",
       dpi=600,
       units="px", 
       height=3500,
       width= 6000
)




### Red Angus: --------
query_ids<-which(fam_file_plus$Breed == "RedAngus" | fam_file_plus$Breed=="AngusRed" )
length(query_ids)

filter_file<- data.frame(fam_file_plus$SampleID[query_ids], fam_file_plus$SampleID[query_ids])
dim(filter_file)

write.table(filter_file, "/media/Scratch3/Johanna_1000Bull/ScopeAnalysis/Mixed_Breeds/Chr_Merged_BeaglePhased_Filter_RedAngus.txt",
            quote=F,col.names=F, row.names=F, sep="\t")



setwd("/media/Scratch3/Johanna_1000Bull/ScopeAnalysis/Mixed_Breeds/")
famDF = read.table("Chr_Merged_BeaglePhased_RedAngus.fam")
colnames(famDF) = c("Breed", "ID","V3","V4","V5", "V6")
famDF$Breed = as.factor(fam_file_plus[which(as.character(fam_file_plus$SampleID) %in% as.character(famDF$ID)), "Breed"])

admixtureDF2= read.table("Chr_Merged_BeaglePhased_RedAngus_USV_k2Qhat.txt")
admixtureDF3= read.table("Chr_Merged_BeaglePhased_RedAngus_USV_k3Qhat.txt")


final_graphic_RedAng <- createBreedComparisonPlot(
  breedB = "RedAngus",
  breedA = "AngusRed",
  breedC="no",
  fam_file_plus = fam_file_plus,
  sim_by_ibs = sim_by_ibs,
  famDF = famDF,
  admixtureDF_k3 = admixtureDF3,
  admixtureDF_k2 = admixtureDF2,
  n_individuals = 53,
  legend_decision = FALSE
)


ggsave(plot=final_graphic_RedAng,
filename= "/Home/johanna/Dokumente/1000_bull_genomes/Graphics/For_publication/For_Publication_RedAngus.png",
dpi=600,
units="px", 
height=3500,
width= 6000
)




### Wagyu ----------
setwd("/Home/johanna/Dokumente/1000_bull_genomes/SCopr_ancestral/Wagyu/")
famDF = read.table("Chr_Merged_BeaglePhased_wagyu.fam")
colnames(famDF) = c("Breed", "ID","V3","V4","V5", "V6")
admixtureDF2= read.table("Chr_Merged_BeaglePhased_wagyu_k2_USVQhat.txt")
admixtureDF3= read.table("Chr_Merged_BeaglePhased_wagyu_k3_USVQhat.txt")

final_graphic_Wag <- createBreedComparisonPlot(
  breedA = "Wagyu",
  breedB = "WagyuModern",
  fam_file_plus = fam_file_plus,
  sim_by_ibs = sim_by_ibs,
  famDF = famDF,
  admixtureDF_k3 = admixtureDF3,
  admixtureDF_k2 = admixtureDF2,
  n_individuals = 30,
  legend_decision = FALSE
)


ggsave(plot=final_graphic_Wag,
       filename= "/Home/johanna/Dokumente/1000_bull_genomes/Graphics/For_publication/For_Publication_Wagyu.png",
       dpi=600,
       units="px", 
       height=3500,
       width= 6000
)


### Tyrolean breeds: -------
query_ids<-which(fam_file_plus$FID == "GreyCattle" | fam_file_plus$FID =="TyroleanGrauvieh" | fam_file_plus$FID =="TyroleanGrey")


filter_file<- data.frame(fam_file_plus$SampleID[query_ids], fam_file_plus$SampleID[query_ids])

#write.table(filter_file, "/media/Scratch3/Johanna_1000Bull/ScopeAnalysis/Mixed_Breeds/Chr_Merged_BeaglePhased_Filter_TyrolGrey.txt",
#            quote=F,col.names=F, row.names=F, sep="\t")


setwd("/media/Scratch3/Johanna_1000Bull/ScopeAnalysis/Mixed_Breeds/")
famDF = read.table("Chr_Merged_BeaglePhased_TyrolGrey.fam")
colnames(famDF) = c("Breed", "ID","V3","V4","V5", "V6")
famDF$Breed = as.factor(fam_file_plus[which(as.character(fam_file_plus$SampleID) %in% as.character(famDF$ID)), "Breed"])

admixtureDF2= read.table("Chr_Merged_BeaglePhased_TyrolGrey_k2Qhat.txt")
admixtureDF3= read.table("Chr_Merged_BeaglePhased_TyrolGrey_k4Qhat.txt")

final_graphic_TyGrey <- createBreedComparisonPlot(
  breedC = "TyroleanGrey",
  breedB =  "TyroleanGrauvieh",
  breedA = "GreyCattle",
  fam_file_plus = fam_file_plus,
  sim_by_ibs = sim_by_ibs,
  famDF = famDF,
  admixtureDF_k3 = admixtureDF3,
  admixtureDF_k2 = admixtureDF2,
  n_individuals = 32,
  legend_decision = FALSE
)

ggsave(plot=final_graphic_TyGrey,
       filename= "/Home/johanna/Dokumente/1000_bull_genomes/Graphics/For_publication/For_Publication_TyrolGrey.png",
       dpi=600,
       units="px", 
       height=3500,
       width= 6000
)


### Angler breeds: --------------------------------------------------------
##prep file for plink and Scope:
query_ids<-which(fam_file_plus$FID == "ModernAngler" | fam_file_plus$FID =="TraditionalAngler" | fam_file_plus$FID =="GermanRedAngler")

filter_file<- data.frame(fam_file_plus$SampleID[query_ids], fam_file_plus$SampleID[query_ids])

#write.table(filter_file, "/media/Scratch3/Johanna_1000Bull/ScopeAnalysis/Mixed_Breeds/Chr_Merged_BeaglePhased_Filter_Angler.txt",
#           quote=F,col.names=F, row.names=F, sep="\t")

##### Plot of the results:
setwd("/media/Scratch3/Johanna_1000Bull/ScopeAnalysis/Mixed_Breeds/")
famDF = read.table("Chr_Merged_BeaglePhased_Angler.fam")
colnames(famDF) = c("Breed", "ID","V3","V4","V5", "V6")
famDF$Breed = as.factor(fam_file_plus[which(as.character(fam_file_plus$SampleID) %in% as.character(famDF$ID)), "Breed"])

admixtureDF2= read.table("Chr_Merged_BeaglePhased_Angler_k2Qhat.txt")
admixtureDF3= read.table("Chr_Merged_BeaglePhased_Angler_k4Qhat.txt")

final_graphic_Angler <- createBreedComparisonPlot(
  breedA = "ModernAngler",
  breedB = "GermanRedAngler",
  breedC = "TraditionalAngler",
  fam_file_plus = fam_file_plus,
  sim_by_ibs = sim_by_ibs,
  famDF = famDF,
  admixtureDF_k3 = admixtureDF3,
  admixtureDF_k2 = admixtureDF2,
  n_individuals = 31,
  legend_decision = FALSE
)


ggsave(plot=final_graphic_Angler,
       filename= "/Home/johanna/Dokumente/1000_bull_genomes/Graphics/For_publication/For_Publication_Angler.png",
       dpi=600,
       units="px", 
       height=3500,
       width= 6000
)


### Danish Red : ---------------------------------------------
query_ids<-which(fam_file_plus$FID == "ModernDanishRed" |
                   fam_file_plus$FID =="TraditionalDanishRed" |
                   fam_file_plus$FID =="DanishRedDairy")

filter_file<- data.frame(fam_file_plus$SampleID[query_ids], fam_file_plus$SampleID[query_ids])

#write.table(filter_file, "/media/Scratch3/Johanna_1000Bull/ScopeAnalysis/Mixed_Breeds/Chr_Merged_BeaglePhased_Filter_DanishRed.txt",
#           quote=F,col.names=F, row.names=F, sep="\t")

############## Plot of the results:
setwd("/media/Scratch3/Johanna_1000Bull/ScopeAnalysis/Mixed_Breeds/")
famDF = read.table("Chr_Merged_BeaglePhased_DanishRed.fam")
colnames(famDF) = c("Breed", "ID","V3","V4","V5", "V6")
famDF$Breed = as.factor(fam_file_plus[which(as.character(fam_file_plus$SampleID) %in% as.character(famDF$ID)), "Breed"])

admixtureDF2= read.table("Chr_Merged_BeaglePhased_DanishRed_k2Qhat.txt")
admixtureDF3= read.table("Chr_Merged_BeaglePhased_DanishRed_k4Qhat.txt")

final_graphic_DanRed <- createBreedComparisonPlot(
  breedC = "ModernDanishRed",
  breedA = "DanishRedDairy",
  breedB = "TraditionalDanishRed",
  fam_file_plus = fam_file_plus,
  sim_by_ibs = sim_by_ibs,
  famDF = famDF,
  admixtureDF_k3 = admixtureDF3,
  admixtureDF_k2 = admixtureDF2,
  n_individuals = 73,
  legend_decision = FALSE
)


ggsave(plot=final_graphic_DanRed,
       filename= "/Home/johanna/Dokumente/1000_bull_genomes/Graphics/For_publication/For_Publication_DanRed.png",
       dpi=600,
       units="px", 
       height=4500,
       width= 8000
)


## Holstein -- HolsteinFriesian: ---------------
query_ids_hf = which(fam_file_plus$FID =="HolsteinFriesian")
query_ids_hol = sample(which(fam_file_plus$FID =="Holstein"), size=length(query_ids_hf)*2)

query_ids = c(query_ids_hf, query_ids_hol)

filter_file<- data.frame(fam_file_plus$SampleID[query_ids], fam_file_plus$SampleID[query_ids])


#write.table(filter_file, "/media/Scratch3/Johanna_1000Bull/ScopeAnalysis/Mixed_Breeds/Chr_Merged_BeaglePhased_Filter_Holstein.txt",
#           quote=F,col.names=F, row.names=F, sep="\t")

############## Plot of the results:
setwd("/media/Scratch3/Johanna_1000Bull/ScopeAnalysis/Mixed_Breeds/")
famDF = read.table("Chr_Merged_BeaglePhased_Holstein.fam")
colnames(famDF) = c("Breed", "ID","V3","V4","V5", "V6")
famDF$Breed = as.factor(fam_file_plus[which(as.character(fam_file_plus$SampleID) %in% as.character(famDF$ID)), "Breed"])

admixtureDF2= read.table("Chr_Merged_BeaglePhased_Holstein_USV_k2Qhat.txt")
admixtureDF3= read.table("Chr_Merged_BeaglePhased_Holstein_USV_k3Qhat.txt")

final_graphic_Holstein <- createBreedComparisonPlot(
  breedA = "Holstein",
  breedB = "HolsteinFriesian",
  fam_file_plus = fam_file_plus,
  sim_by_ibs = sim_by_ibs,
  famDF = famDF,
  admixtureDF_k3 = admixtureDF3,
  admixtureDF_k2 = admixtureDF2,
  n_individuals = 126,
  legend_decision = FALSE
)


ggsave(plot=final_graphic_Holstein,
       filename= "/Home/johanna/Dokumente/1000_bull_genomes/Graphics/For_publication/For_Publication_Holstein.png",
       dpi=600,
       units="px", 
       height=3500,
       width= 6000
)


# Small Taurine breeds : --------
# Since some of the breeds wre discovered in the first step of the analysis framework,
# the improved metadata has to be loaded for further analyses

fam_file<-read.table("/media/Scratch3/Johanna_1000Bull/ImputierteDateien/Chr1_BeaglePhased.fam")
colnames(fam_file)<-c("FID", "IID", "MAT", "PAT", "SEX", "Pheno")

#meta file:
data = read.delim("/Home/johanna/Dokumente/1000_bull_genomes/Data/Run9-TAUIND_recovery_by_NCBI_information_v3.txt")

# For the visualizations we have to keep the order of the fam file and add the recovered Breeds from the NCBI:
fam_file_plus= merge( y=data[,c(1,8,9)], x=fam_file[,c(1:2)], by.x="IID", by.y = "SampleID")

#we want to re-assure that the fam-file order does not change:
fam_file_plus = fam_file_plus[fam_file$IID, ]
colnames(fam_file_plus)[1] = "SampleID"

all(fam_file_plus$SampleID == fam_file$IID)
all(fam_file_plus$FID == fam_file$FID)

# In several of the functions rownames have to be equal to the indices: 
rownames(fam_file_plus) = 1:nrow(fam_file_plus)

# Change the breed column to characters, for the PCA plotting to work appropriately
fam_file_plus$Breed = as.character(fam_file_plus$Breed)


## Admixture small Taurine breeds: -------------
#Smaller breeds have to be grouped for the analysis with SCOPE die to the sample size limiation

query_ids<-which(fam_file_plus$Breed == "FjällCattle" | fam_file_plus$Breed =="Fjäll" |
                   fam_file_plus$Breed == "Kalmyk" | fam_file_plus$Breed =="Kalmykian" |
                   fam_file_plus$Breed =="BelgianBlue" | fam_file_plus$Breed =="BelgiumBlue" |
                   fam_file_plus$Breed =="LithuanianRed" | fam_file_plus$Breed =="TraditionalLithuanianRed"|
                   fam_file_plus$Breed =="Muturu" | fam_file_plus$Breed =="Mututu")

length(query_ids)


filter_file<- data.frame(fam_file_plus$SampleID[query_ids], fam_file_plus$SampleID[query_ids])
filter_file


#write.table(filter_file, "/media/Scratch3/Johanna_1000Bull/ScopeAnalysis/Mixed_Breeds/Chr_Merged_BeaglePhased_Filter_SmallTaurineBreeds.txt",
quote=F,col.names=F, row.names=F, sep="\t")


### Produce the graphic outputs: 
setwd("/media/Scratch3/Johanna_1000Bull/ScopeAnalysis/Mixed_Breeds/")
famDF = read.table("Chr_Merged_BeaglePhased_SmallTaurine.fam")
colnames(famDF) = c("Breed", "ID","V3","V4","V5", "V6")
famDF$Breed = as.factor(fam_file_plus[which(as.character(fam_file_plus$SampleID) %in% as.character(famDF$ID)), "Breed"])

admixtureDF2= read.table("Chr_Merged_BeaglePhased_SmallTaurineBreeds_USV_k5Qhat.txt")
admixtureDF3= read.table("Chr_Merged_BeaglePhased_SmallTaurineBreeds_USV_k11Qhat.txt")

k5_plotDF=preparePlotData_SCOPE(famDF, admixtureDF2, n=70)
k11_plotDF=preparePlotData_SCOPE(famDF, admixtureDF3, n=70)

k5=createAdmixturePlot(k5_plotDF)
k11=createAdmixturePlot(k11_plotDF)


final_graphic_smallBreeds = ggarrange(k5, k11)

print(final_graphic_smallBreeds)


ggsave(plot=final_graphic_smallBreeds,
       filename= "/Home/johanna/Dokumente/1000_bull_genomes/Graphics/For_publication/For_Publication_SmallTaurineBreeds_k5_11.png",
       dpi=300,
       units="px", 
       height=3000,
       width= 6000
)


## PCA and distance (1-IBS) visualizations for potential breed merges in small breeds -----------------
## Fjäll -- Fjäll Cattle -----------
final_graphic_Fjaell <- createBreedComparisonPlot_small(
  breedA = "Fjäll",
  breedB = "FjällCattle",
  fam_file_plus = fam_file_plus,
  sim_by_ibs = sim_by_ibs,
  legend_decision = FALSE,
  lettering = TRUE,
  label1 = "A",
  label2 = "B"
)
print(final_graphic_Fjaell)



## Mututu and Muturu: --------------
final_graphic_Muturu <- createBreedComparisonPlot_small(
  breedA = "Muturu",
  breedB = "Mututu",
  fam_file_plus = fam_file_plus,
  sim_by_ibs = sim_by_ibs,
  legend_decision = FALSE,
  lettering = TRUE,
  label1 = "C",
  label2 = "D"
  
)
print(final_graphic_Muturu)


## Kalmyk -- Kalmykian -----------------
final_graphic_Kalmyk <- createBreedComparisonPlot_small(
  breedA = "Kalmyk",
  breedB = "Kalmykian",
  fam_file_plus = fam_file_plus,
  sim_by_ibs = sim_by_ibs,
  legend_decision = FALSE,
  lettering = TRUE,
  label1 = "E",
  label2 = "F"
)
print(final_graphic_Kalmyk)

## BelgianBlue -- BelgiumBlue ------
final_graphic_Belgian <- createBreedComparisonPlot_small(
  breedB = "BelgianBlue",
  breedA = "BelgiumBlue",
  fam_file_plus = fam_file_plus,
  sim_by_ibs = sim_by_ibs,
  legend_decision = FALSE,
  lettering = TRUE,
  label1 = "G",
  label2 = "H"
)
print(final_graphic_Belgian)

## Lithuanian Red -- Traditional Lithuanian Red:
final_graphic_LithuanianRed <- createBreedComparisonPlot_small(
  breedA = "LithuanianRed",
  breedB = "TraditionalLithuanianRed",
  fam_file_plus = fam_file_plus,
  sim_by_ibs = sim_by_ibs,
  legend_decision = FALSE,
  lettering = TRUE,
  label1 = "I",
  label2 = "J"
)
print(final_graphic_LithuanianRed)  



final_graphic_Muturu_Fjaell_Kalmyk_Belgian_lithuanian = ggarrange(final_graphic_Fjaell, 
                                                                  final_graphic_Muturu,
                                                                  final_graphic_Kalmyk,
                                                                  final_graphic_Belgian,
                                                                  final_graphic_LithuanianRed,
                                                                  nrow=5,
                                                                  ncol=1)

print(final_graphic_Muturu_Fjaell_Kalmyk_Belgian_lithuanian)

ggsave(plot= final_graphic_Muturu_Fjaell_Kalmyk_Belgian_lithuanian,
       filename= "/Home/johanna/Dokumente/1000_bull_genomes/Graphics/For_publication/For_Publication_Muturu_Fjaell_Kalmyk_Belgian_Lituanian.png",
       dpi=300,
       units="px", 
       height=3000,
       width= 6000
)



## small indicine breeds: ----------------------
query_ids<-which(fam_file_plus$FID == "Tharparkar" | fam_file_plus$FID =="TharparkerModern" |
                   fam_file_plus$FID =="Sahiwal" | fam_file_plus$FID =="Shaiwal" |
                   fam_file_plus$FID =="Nelore" | fam_file_plus$FID =="Ogaden")

dim(fam_file_plus[query_ids, ])

filter_file<- data.frame(fam_file_plus$SampleID[query_ids], fam_file_plus$SampleID[query_ids])
filter_file


#write.table(filter_file, "/media/Scratch3/Johanna_1000Bull/ScopeAnalysis/Mixed_Breeds/Chr_Merged_BeaglePhased_Filter_SmallIndicineBreeds.txt",
quote=F,col.names=F, row.names=F, sep="\t")

### Produce the graphic outputs: 
setwd("/media/Scratch3/Johanna_1000Bull/ScopeAnalysis/Mixed_Breeds/")
famDF = read.table("Chr_Merged_BeaglePhased_SmallIndicine.fam")
colnames(famDF) = c("Breed", "ID","V3","V4","V5", "V6")
famDF$Breed = as.factor(fam_file_plus[which(as.character(fam_file_plus$SampleID) %in% as.character(famDF$ID)), "Breed"])

admixtureDF2= read.table("Chr_Merged_BeaglePhased_SmallIndicineBreeds_USV_k4Qhat.txt")
admixtureDF3= read.table("Chr_Merged_BeaglePhased_SmallIndicineBreeds_USV_k9Qhat.txt")

k4_plotDF=preparePlotData_SCOPE(famDF, admixtureDF2, n=44)
k9_plotDF=preparePlotData_SCOPE(famDF, admixtureDF3, n=44)

k4=createAdmixturePlot(k4_plotDF)
k9=createAdmixturePlot(k9_plotDF)


final_graphic_smallBreeds = ggarrange(k4, k9,  labels = LETTERS[1:2] )

print(final_graphic_smallBreeds)


ggsave(plot=final_graphic_smallBreeds,
       filename= "/Home/johanna/Dokumente/1000_bull_genomes/Graphics/For_publication/For_Publication_SmaallIndicineBreeds_k4_9.png",
       dpi=300,
       units="px", 
       height=3000,
       width= 6000
)


## Sahiwal Shaiwal: ---------
final_graphic_Sahiwal <- createBreedComparisonPlot_small(
  breedA = "Sahiwal",
  breedB = "Shaiwal",
  fam_file_plus = fam_file_plus,
  sim_by_ibs = sim_by_ibs,
  legend_decision = FALSE,
  lettering = TRUE,
  label1 = "A",
  label2 = "B"
)
print(final_graphic_Sahiwal)


## Tharparkar - TharparkerModern: ----------
final_graphic_Tharparkar <- createBreedComparisonPlot_small(
  breedA = "Tharparkar",
  breedB = "TharparkerModern",
  fam_file_plus = fam_file_plus,
  sim_by_ibs = sim_by_ibs,
  legend_decision = FALSE,
  lettering = TRUE,
  label1 = "C",
  label2 = "D"
)
print(final_graphic_Tharparkar)




final_graphic_Sahiwal_Tharparkar = ggarrange(final_graphic_Sahiwal,
                                             final_graphic_Tharparkar,
                                             ncol=1,
                                             nrow=2)

print(final_graphic_Sahiwal_Tharparkar)

ggsave(plot=final_graphic_Sahiwal_Tharparkar,
       filename= "/Home/johanna/Dokumente/1000_bull_genomes/Graphics/For_publication/For_Publication_Sahiwal_Tharparkar.png",
       dpi=600,
       units="px", 
       height=3500,
       width= 6000
)



# Correct the metadata according to the analysis: -------

## Read different metadata versions: -----

#original metadata with all missing information:
meta_original<-read.xlsx("/Home/johanna/Dokumente/1000_bull_genomes/Data/Run9-TAUINDraw-animalDistributionList-20210630.xlsx")
#including the simple NCBI query information:
meta2 = read.delim("/Home/johanna/Dokumente/1000_bull_genomes/Data/Run9-TAUIND_recovery_by_NCBI_information_v2.txt")
#including the breed splits from the NCBI information:
meta3 = read.delim("/Home/johanna/Dokumente/1000_bull_genomes/Data/Run9-TAUIND_recovery_by_NCBI_information_v3.txt")

## Meta-data after correction decisions: -----------
clean_fam=meta3[ ,c("SampleID", "Breed", "Species")]
clean_fam$Species_Johanna = clean_fam$Species

## Subspecies level: correct the subspecies according to the analysis results   -------
clean_fam[grep("BohaiBlack", clean_fam$Breed), "Species_Johanna" ]= "Taurus"
clean_fam[grep("JapaneseNative", clean_fam$Breed), "Species_Johanna" ]= "Taurus"
clean_fam[which(clean_fam$Breed == "Crossbreed" & clean_fam$Species=="other"), "Species_Johanna" ]= "Taurus"

clean_fam[grep("Ankole", clean_fam$Breed), "Species_Johanna" ] = "TaurusXIndicus"
clean_fam[grep("Benishangul",clean_fam$Breed), "Species_Johanna" ]= "TaurusXIndicus"
clean_fam[grep("Nuer",clean_fam$Breed), "Species_Johanna" ]= "TaurusXIndicus"
clean_fam[grep("Lingnan",clean_fam$Breed), "Species_Johanna" ]= "TaurusXIndicus"
clean_fam[grep("JiaxianRed",clean_fam$Breed), "Species_Johanna" ]= "TaurusXIndicus"
clean_fam[grep("Maure",clean_fam$Breed), "Species_Johanna" ]= "TaurusXIndicus"
clean_fam[grep("SichuanIndigenous",clean_fam$Breed), "Species_Johanna" ]= "TaurusXIndicus"
clean_fam[grep("Weining",clean_fam$Breed), "Species_Johanna" ]= "TaurusXIndicus"
clean_fam[grep("Qinchuan",clean_fam$Breed), "Species_Johanna" ]= "TaurusXIndicus"
clean_fam[grep("Gourounsi",clean_fam$Breed), "Species_Johanna" ]= "TaurusXIndicus"

clean_fam[grep("Kenana",clean_fam$Breed), "Species_Johanna" ]= "Indicus"
clean_fam[grep("Butana",clean_fam$Breed), "Species_Johanna" ]= "Indicus"
clean_fam[grep("Jinjiang",clean_fam$Breed), "Species_Johanna" ]= "Indicus"
clean_fam[grep("Wenshan",clean_fam$Breed), "Species_Johanna" ]= "Indicus"
clean_fam[grep("NariMaster",clean_fam$Breed), "Species_Johanna" ]= "Indicus"
clean_fam[grep("Jian",clean_fam$Breed), "Species_Johanna" ]= "Indicus"
clean_fam[grep("Leiqiong",clean_fam$Breed), "Species_Johanna" ]= "Indicus"
clean_fam[grep("Dehong",clean_fam$Breed), "Species_Johanna" ]= "Indicus"
clean_fam[grep("Fujian",clean_fam$Breed), "Species_Johanna" ]= "Indicus"

# correction of subspecies, based on quantiles, given that no breed information was available:
clean_fam[which(clean_fam$Species == "unknown"), "Species_Johanna" ]= "Taurus"
clean_fam[which(clean_fam$Species_Johanna == "PredictedTaurus"), "Species_Johanna" ]= "Taurus"
clean_fam[which(clean_fam$Species_Johanna == "ancient"), "Species_Johanna" ]= "Taurus"
clean_fam[grep("TaurusCloned", clean_fam$Species), "Species_Johanna" ]= "Taurus"

clean_fam[which(clean_fam$SampleID %in% c("Unknown-358", "Unknown-356")), "Species_Johanna" ]= "Taurus"
clean_fam[which(clean_fam$SampleID %in% c("Unknown-357", "Unknown-359")), "Species_Johanna" ]= "TaurusXIndicus"



## Visualize the changes in the subspecies: --------
df <- clean_fam %>%
  #filter(Species == "Taurus" | Species=="Indicus")%>%
  group_by(Species, Species_Johanna) %>%
  filter(Species != Species_Johanna) %>%
  summarise(freq = n(), .groups = 'drop') 


custom_colors <- c( "#56B4E9" , "#009E73", "#E69F00" )

flow = ggplot(df, aes(axis1 = Species, axis2 = Species_Johanna, y = freq)) +
  geom_alluvium(aes(fill = Species_Johanna)) +
  geom_stratum(width = 0.45) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Original Subpecies", "Corrected Subspecies"), ) +
  scale_fill_manual(values = custom_colors, name = "Final Subspecies") +
  theme_minimal()+
  labs(y="Number of changes")+
  ylim(c(0, sum(df$freq) * 1.1))


ggsave(file="/Home/johanna/Dokumente/1000_bull_genomes/Graphics/For_publication/Flow_plot_subspecies.png",
       plot=flow,
       dpi=300, 
       width=2500,
       height=3000, unit="px")


table(clean_fam$Species_Johanna)

# Breed assignments:--------------------------
clean_fam$Breed_Johanna = as.character(clean_fam$Breed)

clean_fam[grep("AngusRed", clean_fam$Breed ), "Breed_Johanna" ]= "RedAngus"
clean_fam[grep("AyrshireFinnish", clean_fam$Breed ), "Breed_Johanna" ]= "FinnishAyrshire"

clean_fam[grep("WagyuModern", clean_fam$Breed), "Breed_Johanna" ]= "Wagyu"
clean_fam[grep("TharparkerModern", clean_fam$Breed ), "Breed_Johanna" ]= "Tharparkar"

clean_fam[grep("DanishRedDairy", clean_fam$Breed), "Breed_Johanna" ]= "ModernDanishRed"
clean_fam[grep("TraditionalAngler",clean_fam$Breed), "Breed_Johanna" ]= "GermanRedAngler"
clean_fam[grep("FjällCattle",clean_fam$Breed), "Breed_Johanna" ]= "Fjäll"
clean_fam[grep("Kalmyk",clean_fam$Breed), "Breed_Johanna" ]= "Kalmykian" 
clean_fam[grep("HolsteinFriesian",clean_fam$Breed), "Breed_Johanna" ]= "Holstein" 

clean_fam[grep("TyroleanGrauvieh", clean_fam$Breed), "Breed_Johanna" ]= "TyroleanGrey"
clean_fam[grep("GreyCattle", clean_fam$Breed), "Breed_Johanna" ]= "TyroleanGrey"
clean_fam[grep("BelgiumBlue", clean_fam$Breed), "Breed_Johanna" ]= "BelgianBlue"

clean_fam[grep("Shaiwal", clean_fam$Breed ), "Breed_Johanna" ]= "Sahiwal"
clean_fam[grep("Mututu", clean_fam$Breed ), "Breed_Johanna" ]= "Muturu"
clean_fam[grep("Ndama", clean_fam$Breed ), "Breed_Johanna" ]= "NDama"


#add samples from the NCBI:
clean_fam[grep("Maine Anjou", data$ncbi_infraspecies), "Breed_Johanna"]= "MaineAnjou"
clean_fam[grep("Blonde d'Aquitaine", data$ncbi_infraspecies), "Breed_Johanna"]= "BlondedAquitaine"
clean_fam[grep("Rashoki", data$ncbi_infraspecies), "Breed_Johanna"]= "Rashoki"


nrow(table(clean_fam$Breed))
dim(table(clean_fam$Breed_Johanna))

View(table(data$Breed))




#count individuals with any kind of change: --------------
clean_fam$status = rep("unchanged", len=nrow(clean_fam))
clean_fam$original_breed=meta_original$Breed
clean_fam$original_species=meta_original$Species

for(i in 1: nrow(clean_fam)){
  if(clean_fam$original_breed[i] != clean_fam$Breed_Johanna[i]){
    clean_fam$status[i] = "changed"
  }
  
  if(clean_fam$original_species[i] != clean_fam$Species_Johanna[i]){
    clean_fam$status[i] = "changed"
  }
  
}

table(clean_fam$status)
round(588/nrow(clean_fam), digits=2)
