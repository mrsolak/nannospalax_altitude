###############DADA2 Assign Taxonomy
library(Biostrings)
library(phyloseq)
library(dada2)
rm(list = ls())
load("D:/Rwd/PhD18s/PS_diet_filtered.R")


fasta=refseq(PS.filtered)
taxa <- assignTaxonomy(as.character(fasta), "D:/18S_DATA/silva_132.18s.99_rep_set.dada2.fa.gz",multithread=8,minBoot = 50) #download here:https://zenodo.org/record/4587955

save(taxa, file="/home/hsolak/16s2/chimera_taxonomy_assign/taxa.R")


#/truba/home/hsolak/data/18sdiet/silva_132.18s.99_rep_set.dada2.fa.gz
#/truba/home/hsolak/data/18sdiet/PS_diet_filtered.R
##############################################################################



#merge ASV up to family level
###check chordata blast them

rm(list = ls())
setwd("D:/Rwd/PhD18s")
load("PHYLOSEQ_SD_DIET_RV_all_procisteny_TAXO_derep.R")
PHY=PHYLOSEQ_SD_DIET_RV_all_procisteny_TAXO_derep

PS=prune_samples(sample_data(PHY)$Dataset=="Halil", PHY)

#remove ASVs
inc=taxa_sums(PS)>0
PS=prune_taxa(inc, PS)

#modify the metadata
write.csv(sample_data(PS), "D:/Rwd/PhD18s/metaorg2.csv")
meta=read.csv("D:/Rwd/PhD18s/metaorg2.csv")
meta=sample_data(meta)
sample_names(meta)=meta$X
#_
X= otu_table(PS)
Y=tax_table(PS)
Z=refseq(PS)
PS_new=merge_phyloseq(X,Y,Z,meta)

##exclude pos and neg controls
'%!in%' <- function(x,y)!('%in%'(x,y))
PS_new=subset_samples(PS_new, alt %!in% "pos"==T) 
PS_new=subset_samples(PS_new, alt %!in% "neg"==T) 

PS=PS_new

save(PS, file="D:/Rwd/PhD18s/PS_diet.R")
#######################################################################################
###############################FILTERING##############################################
#######################################################################################
rm(list = ls())
load("D:/Rwd/PhD18s/PS_diet.R")

#exclude diet unrelated class 
exclude_class<-c(" Atractiellomycetes", " Cystobasidiomycetes"," Exobasidiomycetes", " Microbotryomycetes", " Moniliellomycetes"," Pucciniomycetes", " Tritirachiomycetes", " Ustilaginomycetes", " Wallemiomycetes", " Saccharomycetes" , " Eurotiomycetes" )

PS.filtered<-prune_taxa(!tax_table(PS)[,"Class"]%in%exclude_class,PS)


#exclude NA phyla
FILTER<-as.logical(!is.na(tax_table(PS)[,2]))
PS.filtered<-prune_taxa(FILTER,PS)

sum(otu_table(PS.filtered))/sum(otu_table(PS))
#0.85 ~15percent of the data is gone

sum(otu_table(PS.filtered))
#2.950.756

#tremadota filter
Trematoda.filter<-as.logical(tax_table(PS.filtered)[,3]==" Trematoda")
Trematoda.filter[is.na(Trematoda.filter)]<-FALSE

PS.filtered=prune_taxa(Trematoda.filter!=TRUE,PS.filtered)

#Prototheca filter
Prototheca.filter<-as.logical(tax_table(PS.filtered)[,6]==" Prototheca")
Prototheca.filter[is.na(Prototheca.filter)]<-FALSE

PS.filtered=prune_taxa(Prototheca.filter!=TRUE,PS.filtered)

#mammalia filter
Mammalia.filter<-as.logical(tax_table(PS.filtered)[,3]==" Mammalia")
Mammalia.filter[is.na(Mammalia.filter)]<-FALSE

PS.filtered
PS.filtered<-prune_taxa(Mammalia.filter!=TRUE,PS.filtered)
PS.filtered
#Mucoromycota filter
Mucoromycota.filter<-as.logical(tax_table(PS.filtered)[,2]==" Mucoromycota")
Mucoromycota.filter[is.na(Mucoromycota.filter)]<-FALSE

PS.filtered<-prune_taxa(Mucoromycota.filter!=TRUE,PS.filtered)

#exclude unrelated phyla
exclude_phylum<-c(" Acanthocephala", " Apicomplexa"," Bacillariophyta", " Blastocladiomycota", " Cercozoa"," Ciliophora", " Endomyxa", " Evosea", " Chytridiomycota", " Imbricatea" , " Nematoda", " Oomycota" , " Parabasalia", " Preaxostyla", " Tubulinea", " Zoopagomycota", " NA" )
PS.filtered<-prune_taxa(!tax_table(PS.filtered)[,"Phylum"]%in%exclude_phylum,PS.filtered)

sum(otu_table(PS.filtered))/sum(otu_table(PS))
#73% of the data left after filtering

#save phyloseq
save(PS.filtered, file="D:/Rwd/PhD18s/PS_diet_filtered.R")

#######################################################################################
##########################CHECK THE OTUs ################################
#######################################################################################
#remotes::install_github("adrientaudiere/MiscMetabar")
library(MiscMetabar)
rm(list = ls())
load("D:/Rwd/PhD18s/PS_diet_filtered.R")
PS=PS.filtered

### Basidiomycota, chordata and arthropoda excluded
PS=prune_taxa(!tax_table(PS)[,"Phylum"]%in%" Chordata",PS)
PS=prune_taxa(!tax_table(PS)[,"Phylum"]%in%" Arthropoda",PS)
PS=prune_taxa(!tax_table(PS)[,"Phylum"]%in%" Basidiomycota",PS)

## save final phyloseq ##
PS_diet=PS
sample_names(PS_diet)=c("A03_Newkit","A03","A04","A05",
                   "A06","A07","A08",
                   "A09","A10_Newkit","A10",
                   "DB-1","DB-2","DB-4", "DB-8",
                   "E9-1","E9-3", "E9-4","ER-1",
                   "ER-2", "ER-3","ER-3_Newkit_fr",
                   "ER-4", "K3","K6",               
                   "KZ-1_Newkit_fr", "KZ-1","KZ-3",
                   "KZ-4","KZ-5_Newkit_fr", "KZ-5",         
                   "KZ-6","KZ-7","KZ-8",
                   "KZ-9","U01","U03",
                   "U06","U09")



save(PS_diet, file="D:/Rwd/phyloseq_fin/diet/PS_diet_fin.R")




###############################
##MANTARLARI KONTROL ET, GEREKİRSE MAYA VE KÜFLERİ EXCLUDE ET
##BRAY MATRIXLERINI KARŞILAŞTIR DIET VE MICROBIOME
##AEROBIC VE ANAEROBIC BAKTERILERIN ORANLARINI KIYASLA
##SONUÇLARI PHD RESULTS DOSYASINA EKLE
##DOSYAYI EKIBE GÖNDER
##SON ANALIZ TAVSIYELERI VARSA AL VE MAKALEYI YAZMAYA BAŞLA

#######################################################################################
##########################NO OF SEQS PER SAMPLE ################################
#######################################################################################
rm(list = ls())
load("D:/Rwd/phyloseq_fin/diet/PS_diet_fin.R")
PS=PS_diet

SSUMS<-sample_sums(PS) #Returns the total number of individuals observed from each sample

DF<-data.frame(SSUMS,sample_data(PS)) #Create data frame with sample names
DF<-DF[order(DF$SSUMS),] #Order by sample sums, lowest 1k and highest 210k

write.table(DF,"D:/Rwd/PhD18s/tables/2SSUMS_SUMMARY.txt",row.names = F,quote = F, sep="\t")



#######################################################################################
############################# VIRSUALIZATION #################################################
#######################################################################################
rm(list = ls())
##custom functions

prepare_tax_df<-function(PHYLOSEQ,RANK,Unass.symbol=NA,Unass.repl,
                         min_prop,top_tax,bellow_top,
                         merge_categories,prop.trans,sort_abu){
  if(prop.trans==T){PHYLOSEQ=transform_sample_counts(PHYLOSEQ,function(x) x/sum(x))}
  TT<-tax_table(PHYLOSEQ)
  class(TT)<-"matrix"
  TT[is.na(tax_table(TT))]<-Unass.repl
  TT<-tax_table(TT)
  tax_table(PHYLOSEQ)<-TT
  
  PHYLOSEQ.merged<-tax_glom(PHYLOSEQ,taxrank=RANK)
  PHYLOSEQ.merged=transform_sample_counts(PHYLOSEQ.merged,function(x) x/sum(x))
  
  #select top taxa
  NAMES<-as.character(tax_table(PHYLOSEQ.merged)[,RANK])
  TS<-taxa_sums(PHYLOSEQ.merged)
  names(TS)<-NAMES
  TS<-rev(sort(TS))
  KEEP<-names(TS)[1:top_tax]
  
  #
  TT<-tax_table(PHYLOSEQ.merged)
  class(TT)<-"matrix"
  TT[,RANK][!TT[,RANK]%in%KEEP]<-bellow_top
  TT<-tax_table(TT)
  tax_table(PHYLOSEQ.merged)<-TT
  PHYLOSEQ.merged<-tax_glom(PHYLOSEQ.merged,taxrank=RANK)
  
  mdf = psmelt(PHYLOSEQ.merged)
  
  if(sort_abu==T){
    mdf[,RANK]<-as.factor(mdf[,RANK])
    ADD<-levels(mdf[,RANK])[!levels(mdf[,RANK])%in%KEEP]
    if(length(ADD)>0){mdf[,RANK]<-factor(mdf[,RANK],levels=c(KEEP,ADD))}
    if(length(ADD)==0){mdf[,RANK]<-factor(mdf[,RANK],levels=c(KEEP))}
  }
  
  mdf
}

small_fig<-function(x){
  x+theme_bw()+theme(legend.title = element_blank(),
                     axis.title = element_text(size=8),
                     axis.text =  element_text(size=8),
                     legend.text = element_text(size=8))
}

#####################

load("D:/Rwd/phyloseq_fin/diet/PS_diet_fin.R")
PS=PS_diet
total = median(sample_sums(PS))
standf = function(x, t=total) round(t * (x / sum(x)))

###alt order barplot
PS2_alt <- merge_samples(PS, "alt")
PS2_alt = transform_sample_counts(PS2_alt, standf)

#This palette is relatively fine, if you have < 25 colors:
c25 <- c("gray70","green4", # red
         "#E31A1C",
         "#6A3D9A", # purple
         "#FF7F00", # orange
         "black","dodgerblue2",
         "skyblue2","#FB9A99", # lt pink
         "palegreen2",
         "#CAB2D6", # lt purple
         "#FDBF6F", # lt orange
         "gold1", "khaki2",
         "maroon","orchid1","deeppink1","blue1","steelblue4",
         "darkturquoise","green1","yellow4","yellow3",
         "darkorange4","brown")

p3=plot_bar(PS2_alt, fill = "Order") + geom_bar(stat="identity", position="stack") +scale_fill_manual(values=c25)

# Change order of X axis by altitude group
p3=p3 + aes(x = fct_reorder(Sample, p3$data$alt2))


######Bar Plot -- Phylum by Population
PS2_pop <- merge_samples(PS, "pop")
PS2_pop = transform_sample_counts(PS2_pop, standf)


p4=plot_bar(PS2_pop, fill = "Order") + geom_bar(stat="identity", position="stack") +scale_fill_manual(values=c25)

# Change order of X axis by altitude group
library(forcats)
p4=p4 + aes(x = fct_reorder(Sample, p4$data$alt2))

#combine plots
library(ggpubr)
ggarrange(p3, p4, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")




PS = transform_sample_counts(PS, standf)

#Genus by alt
PS2 <- merge_samples(PS, "Dataset")
PS2 = transform_sample_counts(PS2, standf)


DF.genus<-prepare_tax_df(PHYLOSEQ=PS2,
                                     RANK="Genus",Unass.symbol=NA,
                                     Unass.repl="Unassigned",min_prop=NULL,top_tax=9,
                                     bellow_top="Others",merge_categories=NULL,
                                     prop.trans=TRUE,sort_abu=TRUE)


DF.genus$Genus<-gsub(" group","",DF.genus$Genus)
ta3 = ggplot(DF.genus, aes_string(x = "Dataset", y = "Abundance", fill = "Genus",order="Genus"))+theme_bw(base_size = 12)
ta3 = ta3 + geom_bar(stat = "identity", position = "stack")
ta3= ta3 + theme(axis.title.x=element_blank(),axis.text.x = element_text(size = 7, angle = 90),axis.text.y = element_text(hjust = 0,size=10))

library(forcats)
ta3 + aes(x = fct_reorder(Sample, ta3$data$alt2)) 

#Family by alt

PS_alt <- merge_samples(PS, "alt")
PS_alt = transform_sample_counts(PS_alt, standf)


DF.alt<-prepare_tax_df(PHYLOSEQ=PS_alt,
                                     RANK="Family",Unass.symbol=NA,
                                     Unass.repl="Unassigned",min_prop=NULL,top_tax=9,
                                     bellow_top="Others",merge_categories=NULL,
                                     prop.trans=TRUE,sort_abu=TRUE)


DF.alt$Phylum<-gsub(" group","",DF.alt$Phylum)
ta3 = ggplot(DF.alt, aes_string(x = "alt", y = "Abundance", fill = "Family",order="Family"))+theme_bw(base_size = 12)
ta3 = ta3 + geom_bar(stat = "identity", position = "stack")
ta3= ta3 + theme(axis.title.x=element_blank(),axis.text.x = element_text(size = 7, angle = 90),axis.text.y = element_text(hjust = 0,size=10))
library(forcats)
ta3 + aes(x = fct_reorder(Sample, ta3$data$alt2)) 


#Genus by alt

DF.alt<-prepare_tax_df(PHYLOSEQ=PS_alt,
                       RANK="Genus",Unass.symbol=NA,
                       Unass.repl="Unassigned",min_prop=NULL,top_tax=9,
                       bellow_top="Others",merge_categories=NULL,
                       prop.trans=TRUE,sort_abu=TRUE)


DF.alt$Genus<-gsub(" group","",DF.alt$Genus)
ta3 = ggplot(DF.alt, aes_string(x = "alt", y = "Abundance", fill = "Genus",order="Genus"))+theme_bw(base_size = 12)
ta3 = ta3 + geom_bar(stat = "identity", position = "stack")
ta3= ta3 + theme(axis.title.x=element_blank(),axis.text.x = element_text(size = 7, angle = 90),axis.text.y = element_text(hjust = 0,size=10))
library(forcats)
ta3 + aes(x = fct_reorder(Sample, ta3$data$alt2)) 

#######################################################################################
############################# TOP TAXA #################################################
#######################################################################################
#https://github.com/joey711/phyloseq/issues/494


############################## FAMILY(TOP 5) #############################
phy <- transform_sample_counts(PS, function(x) x/sum(x)) # get abundance in %
glom <- tax_glom(phy, taxrank = 'Family') # agglomerate taxa
dat <- psmelt(glom) # create dataframe from phyloseq object
dat$Family <- as.character(dat$Family) # convert Phylum to a character vector from a factor because R
medians <- ddply(dat, ~Family, function(x) c(median=median(x$Abundance))) # group dataframe by Phylum, calculate median rel. abundance
remainder <- medians[medians$median <= 0.01,]$Family # find Phyla whose rel. abund. is less than 1%
dat[dat$Family %in% remainder,]$Family <- 'Remainder' # change their name to "Remainder"

p= ggplot(dat,
          aes(x=reorder(Family, Abundance),
              y=Abundance)) + geom_boxplot() + coord_flip() + ggtitle ("The most abundant bacterial family")

p
get_taxa_unique(PS, "Family")
############################## GENUS(TOP 5) #############################
PS=PS2
phy <- transform_sample_counts(PS, function(x) x/sum(x)) # get abundance in %
glom <- tax_glom(phy, taxrank = 'Genus') # agglomerate taxa
dat <- psmelt(glom) # create dataframe from phyloseq object
dat$Genus <- as.character(dat$Genus) # convert Genus to a character vector from a factor because R
medians <- ddply(dat, ~Genus, function(x) c(median=median(x$Abundance))) # group dataframe by Genus, calculate median rel. abundance
remainder <- medians[medians$median <= 0.01,]$Genus # find Genus whose rel. abund. is less than 1%
dat[dat$Genus %in% remainder,]$Genus <- 'Remainder' # change their name to "Remainder"

p=ggplot(dat,
         aes(x=reorder(Genus, Abundance),
             y=Abundance)) + geom_boxplot() + coord_flip() + ggtitle ("The most abundant family genera")

p

get_taxa_unique(PS, "Genus")




#######################################################################################
############################# ALPHA D #################################################
#######################################################################################
rm(list = ls())
load("D:/Rwd/phyloseq_fin/diet/PS_diet_fin.R")
PS=PS_diet
rm(PS_diet)

#Alpha Diversity 1
a_div=estimate_richness(PS, measures=c("Observed", "Shannon", "Simpson"))
a_div$ID=sample_names(PS)


write.csv(a_div,"D:/Rwd/PhD18S/tables/alpha_table.csv",row.names = F,quote = F)    
                    
#alpha diversity plots
##ALSO CHECK THIS OUT https://rpubs.com/lconteville/713954
p=plot_richness(PS, x="alt", measures=c("Observed", "Simpson", "Shannon"))+geom_boxplot(alpha=0.6)+ theme_classic()
library(forcats)
p=p+aes(x = fct_reorder(alt, p$data$alt2))
comp=list(c("low","middle"),c("middle","high"),c("low","high"))
symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

p+stat_compare_means(method = "wilcox.test", comparisons = comp, label = "p.signif", symnum.args = symnum.args)


#### mixed modelling ALPHA
PHYLOSEQ2.rare<-rarefy_even_depth(PS)

RICH<-estimate_richness(PHYLOSEQ2.rare, measures=c("Shannon", "Simpson"))
RICH<-data.frame(RICH,sample_data(PHYLOSEQ2.rare))

#1. fitting models
model_null<-glmmTMB(Shannon~1+(1|pop),data=RICH) #null model
model_cat<-glmmTMB(Shannon~alt+(1|pop),data=RICH) #

#2. model comparission
anova(model_null,model_cat) #NS


#both wilcox test and mixed modelling not significant

#######################################################################################
############################# BETA D #################################################
#######################################################################################
rm(list = ls())
load("D:/Rwd/phyloseq_fin/diet/PS_diet_fin.R")
PS=PS_diet
rm(PS_diet)


##Beta Diversity: Plot the PCoA using Bray-Curtis distance:

bray_dist= phyloseq::distance(PS, method="bray", weighted=F)
ordination = ordinate(PS, method="PCoA", distance="bray")
p=plot_ordination(PS, ordination, color="alt") + theme(aspect.ratio=1)+ ggtitle("PCoA using the Bray-Curtis distance (PERMANOVA p=0.001)") + geom_point(size=2)

plot_ordination(PS, ordination, color="alt", shape = "pop") + 
  theme(aspect.ratio=1, legend.text = element_text(size = 12)) +
  geom_point(size=3)+scale_shape_manual(values=c(0,15,4,21,17,3))

#0,4,21,3,15,17

## PERMANOVA ##
adonis2(bray_dist ~ sample_data(PS)$alt)     ##significant

### MIXED MODELLING ###
PHYLOSEQ2.rare<-rarefy_even_depth(PS)
PHYLOSEQ2.trans<-transform_sample_counts(PHYLOSEQ2.rare,function(x) x/sum(x))
DIST<-vegdist(otu_table(PHYLOSEQ2.trans))
ORD<-ordinate(PHYLOSEQ2.trans,method="PCoA",distance = DIST)
plot_ordination(PHYLOSEQ2.trans,ORD, color="alt", shape ="pop") + theme(aspect.ratio=1) + geom_point(size=2)


ORD_df<-plot_ordination(PHYLOSEQ2.trans,ORD,justDF = T)

#Use the PCoA axes as beta diversity model
model_null<-glmmTMB(Axis.1~1+(1|pop),data=ORD_df) 
model_cat<-glmmTMB(Axis.1~alt+(1|pop),data=ORD_df)

summary(model_cat)
anova(model_null,model_cat) #cat model explains better than null p=0.02


#######################################
#Composition MDMR######################
#######################################
# MDMR calculates mixed models for distance data
# should be definetly used - models with PCoA capture only relatively low fraction of variation in the microbiota (11 and 10% for the first and the second axix, respectively)
# but information you can see in the MDMR output are limited
# you probably conduct spaparate analyses to compare high ~ middle, low ~ middle et samples. 

#altitude mdmr
mdmr.res <- MDMR::mixed.mdmr(~alt+
                               (1|pop),
                             data =ORD_df,D = DIST)
summary(mdmr.res) #significant 0.0001




######Prepare matrix for comparison microsatellite and microbiome##############
###############################################################################
#read datasets
#Microbiome
rm(list = ls())
'%!in%' <- function(x,y)!('%in%'(x,y))
load("D:/Rwd/phyloseq_fin/PS.R") #load phyloseq data
PS= subset_samples(PS, altitude %!in% "out"==T) #exclude out samples

#diet
load("D:/Rwd/phyloseq_fin/diet/PS_diet_fin.R")
PS=PS_diet
rm(PS_diet)


##Check sample IDs and get same samples for each dataset
DT=sort(sample_names(PS))
MB=sort(sample_names(PS))

PS= subset_samples(PS, sample_names(PS) %!in% c("DB-1", "A03_Newkit","A10_Newkit","ER-3_Newkit_fr","KZ-1_Newkit_fr","KZ-5_Newkit_fr")==T)
sample_names(PS)[19]="K03"
sample_names(PS)[20]="K06"


PS= subset_samples(PS, sample_ID %!in% c("DB-3","DB-5","DB-6","DB-7","E9-2",
                                                   "K01","K02","K04","K05","K08","K09","K10","K11",
                                                   "KZ-2","U02","U04","U05","U07","U08")==T)
DT=sort(sample_names(PS))
MB=sort(sample_names(PS))

DT%in%MB  ##ALL OK

######################################microbiome matrix
MBIO_PC= ordinate(PS, method="PCoA", distance="bray")


#calculate eucledian distances from the coordinates, distance matrix
library(sf)
library(geosphere)

x <- data.frame(MBIO_PC$vectors[,1])
y <- data.frame(MBIO_PC$vectors[,2])
df <- data.frame(longitude = x, latitute= y)
MBIO_MAT <- distm(df, fun= distGeo)

######################################diet matrix
DIET_PC= ordinate(PS, method="PCoA", distance="bray")

x <- data.frame(DIET_PC$vectors[,1])
y <- data.frame(DIET_PC$vectors[,2])
df <- data.frame(longitude = x, latitute= y)
DIET_MAT <- distm(df, fun= distGeo)

#Mantel Test for Similarity of Two Matrices
library(ape)

mantel.test(MBIO_MAT, DIET_MAT, nperm = 999, alternative="two.sided") #p=0.30

