library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(picante)
library(ape)
library(lavaan)
library(indicspecies)

##function
#convert phyloseq object to vegan matrix
veganotu = function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}

##import biom and tree file
bird_phy = import_biom('16S_dada2.cluster.lulu.biom', '16S_dada2.cluster.lulu.tree.phy')

##import metadata
meta.df = import_qiime_sample_data('bird_metadata_microbiota.txt')
##merge datasets
bird_phy2 = merge_phyloseq(bird_phy, meta.df)

##remove the pos.con sample
bird_phy2 = subset_samples(bird_phy2, phinchID != "pos.control"  )

#Separate out cloaca dataset
bird_phy_clo = subset_samples(bird_phy2, Sample_Type == "cloaca") 
bird_phy_clo = prune_taxa(taxa_sums(bird_phy_clo)>=1, bird_phy_clo) ##3698 taxa and 79 samples

##look at how many levels of the factor "pre_post_AH" exist in the dataset
unique(sample_data(bird_phy_clo)$pre_post_AH)

#remove the chloroplasts and the mitochondria and archea (just keep bacteria)
bird_phy_clo2 = subset_taxa(bird_phy_clo, Rank3!="c__Chloroplast") #2065 taxa and 32 samples
bird_phy_clo2 = subset_taxa(bird_phy_clo2, Rank5!="f__Mitochondria") # 1167 taxa and 32 samples
bird_phy_clo2 = subset_taxa(bird_phy_clo2, Rank1 == "d__Bacteria") # 1164 taxa and 32 samples
bird_phy_clo2 = prune_taxa(taxa_sums(bird_phy_clo2)>=1, bird_phy_clo2)

#look at the number of reads on a sample-by-sample basis to decide which sequencing depth for subsampling 
sample_sums(bird_phy_clo2)

##subsample the phyloseq object to 980 reads
#set.seed(711)
SUBSAMPLED_DATA  = rarefy_even_depth(bird_phy_clo2, sample.size = 980, rngseed = 903, replace = FALSE, trimOTUs = 1)

#check how many taxa are represented by sample
specnumber(veganotu(SUBSAMPLED_DATA))

##get rid of sample "post.53" (only has 7 taxa)<-most of this sample is made
## up of reads from one OTU (OTU2) which is Acinetobacter radioresistans)
SUBSAMPLED_DATA = subset_samples(SUBSAMPLED_DATA, X.SampleID != "post.53") 
#remove any taxa that were unique to the sample we just removed so we don't have taxa with 0 reads in our data
SUBSAMPLED_DATA = prune_taxa(taxa_sums(SUBSAMPLED_DATA)>0, SUBSAMPLED_DATA) 

##get rid of sample "post.08" (only has 15 taxa)<-most of this sample is made
## up of reads from one OTU (OTU3) which is Diplorickettsia massiliensis)
SUBSAMPLED_DATA = subset_samples(SUBSAMPLED_DATA, X.SampleID != "post.08") 
SUBSAMPLED_DATA = subset_samples(SUBSAMPLED_DATA, X.SampleID != "pre.05") 
SUBSAMPLED_DATA = subset_samples(SUBSAMPLED_DATA, X.SampleID != "pre.24") 
SUBSAMPLED_DATA = subset_samples(SUBSAMPLED_DATA, X.SampleID != "pre.42") 
SUBSAMPLED_DATA = subset_samples(SUBSAMPLED_DATA, X.SampleID != "pre.61") 
SUBSAMPLED_DATA = subset_samples(SUBSAMPLED_DATA, X.SampleID != "pre.62") 
SUBSAMPLED_DATA = subset_samples(SUBSAMPLED_DATA, X.SampleID != "pre.64") 
#remove any taxa that were unique to the sample we just removed so we don't have taxa with 0 reads in our data
SUBSAMPLED_DATA = prune_taxa(taxa_sums(SUBSAMPLED_DATA)>0, SUBSAMPLED_DATA) #1561 taxa and 65 samples


#find the number of samples and OTUs for each group
pre_AHY = subset_samples(SUBSAMPLED_DATA, pre_post_AH == "pre_Y") 
pre_AHY = prune_taxa(taxa_sums(pre_AHY)>0, pre_AHY)  #691 taxa and 11 samples

pos_AHY = subset_samples(SUBSAMPLED_DATA, pre_post_AH == "pos_Y")
pos_AHY = prune_taxa(taxa_sums(pos_AHY)>0, pos_AHY)  #611 taxa and 13 samples

pre_AHN = subset_samples(SUBSAMPLED_DATA, pre_post_AH == "pre_N") 
pre_AHN = prune_taxa(taxa_sums(pre_AHN)>0, pre_AHN)  #870 taxa and 26 samples

pos_AHN = subset_samples(SUBSAMPLED_DATA, pre_post_AH == "pos_N") 
pos_AHN = prune_taxa(taxa_sums(pos_AHN)>0, pos_AHN)  #680 taxa and 15 samples

#get datasets for post trt datapoint
post_all = subset_samples(SUBSAMPLED_DATA, pre_post == "pos")
post_all = prune_taxa(taxa_sums(post_all)>0, post_all)  # 986 taxa and 28 samples

#get datasets for all AHY (pre and post)
all_AHY = subset_samples(SUBSAMPLED_DATA, AH == "Y")
all_AHY = prune_taxa(taxa_sums(all_AHY)>0, all_AHY)  # 1008 taxa and 24 samples

#get a relative abundance otu table 
bird_rel_abun = transform_sample_counts(SUBSAMPLED_DATA, function(x) x / sum(x))

#get a relative abundance otu table for posttrt dataset
bird_rel_abun_post = transform_sample_counts(post_all, function(x) x / sum(x))

#get a relative abundance otu table for AHY dataset
bird_rel_abun_AHY = transform_sample_counts(all_AHY, function(x) x / sum(x))


#get metadata
env.data = sample_data(SUBSAMPLED_DATA)
env.data = data.frame(env.data)

env.data.post = sample_data(post_all)
env.data.post = data.frame(env.data.post)

env.data.AHY = sample_data(all_AHY)
env.data.AHY = data.frame(env.data.AHY)

#get vegan OTU table (raw abundance)
votu <- veganotu(SUBSAMPLED_DATA)
votu = as.data.frame(votu)

#get vegan OTU table (raw abundance) --post
votu_post <-veganotu(post_all)
votu_post = as.data.frame(votu_post)

#get vegan OTU table (raw abundance) --AHY
votu_AHY <-veganotu(all_AHY)
votu_AHY = as.data.frame(votu_AHY)

#get vegan OTU table (relative abundance)
rel_abun_votu <-veganotu(bird_rel_abun)
rel_abun_votu = as.data.frame(rel_abun_votu)

#get vegan OTU table (relative abundance) --post
rel_abun_votu_post <-veganotu(bird_rel_abun_post)
rel_abun_votu_post = as.data.frame(rel_abun_votu_post)

#get vegan OTU table (relative abundance) --AHY
rel_abun_votu_AHY <-veganotu(bird_rel_abun_AHY)
rel_abun_votu_AHY = as.data.frame(rel_abun_votu_AHY)

#get vegan OTU table (presence-absence)
PA_votu= votu
PA_votu[PA_votu>0] <-1
PA_votu = as.data.frame(PA_votu)

###############################################
##Find indicator species ##RELATIVE ABUNDANCE -- between Post AHY and Post AHN##
###############################################
rel_abun_sosp_indic_post=multipatt(rel_abun_votu_post, env.data.post[["AH"]], func = "IndVal.g", duleg = TRUE, set.seed(5508), control = how(nperm = 999))
summary(rel_abun_sosp_indic_post, indvalcomp = TRUE)
# List of species associated to each combination: 
#   
# Group N  #sps.  3 
# A      B  stat p.value   
# OTU22 0.7096 1.0000 0.842   0.008 **
#   OTU49 0.7540 0.8000 0.777   0.015 * 
#   OTU74 1.0000 0.4667 0.683   0.007 **
#   
#   Group Y  #sps.  1 
# A      B  stat p.value  
# OTU169 0.9872 0.4615 0.675   0.015 *

tax_table(bird_rel_abun_post)[c("OTU22","OTU49","OTU74","OTU169")]
# Taxonomy Table:     [4 taxa by 7 taxonomic ranks]:
# 
# Rank5                    Rank6                 Rank7                             
# OTU22  "f__Microbacteriaceae"   "g__Salinibacterium"  "s__Salinibacterium_xinjiangense" 
# OTU49  "f__Cryptosporangiaceae" "g__Jatrophihabitans" "s__Jatrophihabitans_endophyticus"
# OTU74  "f__Conexibacteraceae"   "g__Conexibacter"     "s__Conexibacter_woesei"          
# OTU169 "f__Moraxellaceae"       "g__Acinetobacter"    "s__Acinetobacter_johnsonii"     

###############################################
##Find indicator species ##RELATIVE ABUNDANCE -- between pre and post AHY birds
###############################################
rel_abun_sosp_indic_AHY=multipatt(rel_abun_votu_AHY, env.data.AHY[["pre_post"]], func = "IndVal.g", duleg = TRUE,set.seed(5508), control = how(nperm = 999))
summary(rel_abun_sosp_indic_AHY, indvalcomp = TRUE)
# List of species associated to each combination: 
#   
#   Group pre  #sps.  4 
# A      B  stat p.value  
# OTU297  0.8640 0.4545 0.627   0.043 *
#   OTU91   1.0000 0.3636 0.603   0.033 *
#   OTU211  1.0000 0.3636 0.603   0.046 *
#   OTU5984 1.0000 0.3636 0.603   0.035 *
tax_table(bird_rel_abun_post)[c("OTU297","OTU91","OTU211","OTU5984")]
# Rank5                              Rank6                  Rank7                         
                          
# OTU297  "f__Nocardioidaceae"               "g__Aeromicrobium"     "s__Aeromicrobium_fastidiosum"
# OTU91   "f__Erysipelotrichaceae"           "g__Clostridium_XVIII" NA                            
# OTU211  "f__Bradyrhizobiaceae"             "g__Bosea"             "s__Bosea_lathyri"            
# OTU5984 "f__Verrucomicrobiaceae"           "g__Luteolibacter"     NA                            

##BETA DIVERSITY--ORDINATION
###############################################
#NMDS#
###############################################

#NMDS with Unifracs

#calculate distance matrix 
#weighted unifrac
wuni_dist<-phyloseq::distance(SUBSAMPLED_DATA, method = "wunifrac")

#ordinate unsing abundance data
wunifrac.nmds<-metaMDS(wuni_dist, dist = "wunifrac", k=2, trymax = 500)


#calculate distance matrix using PRESENCE ABSENCE data
#Raupcrick
raup_dist<-raupcrick(PA_votu,  null="r1", nsimul=999, chase=FALSE)

#ordinate using relative abundance data
raup.nmds<-metaMDS(raup_dist, dist = "raup", k=3, trymax = 300)

#plot

# New facet label names for anthelminthic
AH.labs <- c("Control (water)", "Anthelminthic")
names(AH.labs) <- c("N", "Y")

#ggplot version 
wuni_plot<- plot_ordination(SUBSAMPLED_DATA, wunifrac.nmds,  type = "samples", shape = 'pre_post',color = 'pre_post') +
  scale_color_manual(values = c("black","black"))+
  stat_ellipse(geom = "polygon", alpha = 0.3, aes(fill = pre_post))+
  geom_point(aes(colour=factor(pre_post), 
                 fill = factor(pre_post)), shape=21, size = 6) + 
  scale_fill_manual(values=c("white", "black"),breaks=c("pos","pre"),
                    labels=c("Pre-treatment", "Post-treatment")) +
  guides(shape = FALSE)+
  guides(color = FALSE)+
  theme(panel.background = element_rect(fill="white", colour="black", size=0.5,
                                        linetype="solid"), panel.grid = element_blank(),
        text = element_text(size=10),strip.background = element_blank(),
        legend.title=element_blank(),legend.position="none",legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 20),
        strip.text = element_text(size = 14))+
  geom_path(aes(group=phinchID), color = 'black',size = 0.9)+
  # geom_text(aes(label=phinchID))+ #just a check to make sure these paths are actually connecting the right points
  facet_grid(~AH,labeller = labeller(AH = AH.labs))+
  labs(title="Weighted Unifrac")

wuni_plot


raup_plot <- plot_ordination(SUBSAMPLED_DATA, raup.nmds,  type = "samples", shape = 'pre_post',color = 'pre_post') +
  scale_color_manual(values = c("black","black"))+
  stat_ellipse(geom = "polygon", alpha = 0.3, aes(fill = pre_post))+
  geom_point(aes(colour=factor(pre_post), 
                 fill = factor(pre_post)), shape=21, size = 6) + 
  scale_fill_manual(values=c("white", "black"),breaks=c("pos","pre"),
                    labels=c("Pre-treatment", "Post-treatment")) +
  guides(shape = FALSE)+
  guides(color = FALSE)+
  theme(panel.background = element_rect(fill="white", colour="black", size=0.5,
                                        linetype="solid"), panel.grid = element_blank(),
        text = element_text(size=10),strip.background = element_blank(),
        legend.title=element_blank(),legend.position="bottom",legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 20),
        strip.text = element_text(size = 14))+
  geom_path(aes(group=phinchID), color = 'black',size = 0.9)+
  # geom_text(aes(label=phinchID))+ #just a check to make sure these paths are actually connecting the right points
  facet_grid(~AH,labeller = labeller(AH = AH.labs))+
  labs(title="Raup-Crick")
raup_plot

if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
#Make figure 4
#layout(matrix(c(1,2), 24, 2, byrow = TRUE))
# tiff("Figure4.tiff", units="in", width=6.5, height=10, res=300)
# 
# 
# plot_grid(wuni_plot, raup_plot , 
#           labels = c("A", "B"),
#           ncol = 1, nrow = 2)
# dev.off()



##PERMANOVA
#############################################################################################
#get ordinations for sp, genus and phylum levels for post-AH-treatment birds
wuni_dist_sp_post <-phyloseq::distance(post_all, method = "wunifrac")
wunifrac.nmds.sp.post<-metaMDS(wuni_dist_sp_post, dist = "wunifrac", k=2, trymax = 300)

PA_votu_post_sp= votu_post
PA_votu_post_sp[PA_votu_post_sp>0] <-1
PA_votu_post_sp = as.data.frame(PA_votu_post_sp)

raup_dist_post_sp<-raupcrick(PA_votu_post_sp,  null="r1", nsimul=999, chase=FALSE)
raup.nmds.post.sp<-metaMDS(raup_dist_post_sp, dist = "raup", k=3, trymax = 300)

post_all_genus = tax_glom(post_all, taxrank = rank_names(post_all)[6])
post_all_phylum = tax_glom(post_all, taxrank = rank_names(post_all)[2])

wuni_dist_genus_post<-phyloseq::distance(post_all_genus, method = "wunifrac")
wunifrac.nmds.genus.post<-metaMDS(wuni_dist_genus_post, dist = "wunifrac", k=2, trymax = 300)

wuni_dist_phylum_post<-phyloseq::distance(post_all_phylum, method = "wunifrac")
wunifrac.nmds.phylum.post<-metaMDS(wuni_dist_phylum_post, dist = "wunifrac", k=2, trymax = 300)

#get vegan OTU table (raw abundance) --post
votu_post_genus <-veganotu(post_all_genus)
votu_post_genus = as.data.frame(votu_post_genus)
votu_post_phylum <-veganotu(post_all_phylum)
votu_post_phylum = as.data.frame(votu_post_phylum)

#get vegan OTU table (presence-absence)
PA_votu_genus= votu_post_genus
PA_votu_genus[PA_votu_genus>0] <-1
PA_votu_genus = as.data.frame(PA_votu_genus)

PA_votu_phylum= votu_post_phylum
PA_votu_phylum[PA_votu_phylum>0] <-1
PA_votu_phylum = as.data.frame(PA_votu_phylum)

raup_dist_post_genus<-raupcrick(PA_votu_genus,  null="r1", nsimul=999, chase=FALSE)
raup.nmds.post.genus<-metaMDS(raup_dist_post_genus, dist = "raup", k=3, trymax = 300)

raup_dist_post_phylum<-raupcrick(PA_votu_phylum,  null="r1", nsimul=999, chase=FALSE)
raup.nmds.post.phylum<-metaMDS(raup_dist_post_phylum, dist = "raup", k=3, trymax = 300)

###########adonis for weighted unifrac --species level
set.seed(111)
w_uni_adonis_post_sp= adonis(wuni_dist_sp_post ~ AH, data = env.data.post, permutations = 999)
w_uni_betadisper_post_sp = betadisper(wuni_dist_sp_post, env.data.post$AH,type = "centroid")
set.seed(111)
w_permutest_post_sp = permutest(w_uni_betadisper_post_sp, permutations = 9999) 
###########adonis for weighted unifrac --genus level
set.seed(111)
w_uni_adonis_post_genus= adonis(wuni_dist_genus_post ~ AH, data = env.data.post, permutations = 99)
w_uni_betadisper_post_genus = betadisper(wuni_dist_genus_post, env.data.post$AH,type = "centroid")
set.seed(111)
w_permutest_post_genus = permutest(w_uni_betadisper_post_genus, permutations = 9999) 
###########adonis for weighted unifrac --phylum level
set.seed(111)
w_uni_adonis_post_phylum= adonis(wuni_dist_phylum_post ~ AH, data = env.data.post, permutations = 99)
w_uni_betadisper_post_phylum = betadisper(wuni_dist_phylum_post, env.data.post$AH,type = "centroid")
set.seed(111)
w_permutest_post_phylum = permutest(w_uni_betadisper_post_phylum, permutations = 9999) 


###########adonis for raup --species level
set.seed(111)
raup_adonis_post_sp = adonis (raup_dist_post_sp ~ AH, data = env.data.post, permutations = 99)
raup_betadisper_post_sp = betadisper(raup_dist_post_sp, env.data.post$AH, type = "centroid")
set.seed(111)
raup_permutest_post_sp = permutest(raup_betadisper_post_sp, permutations = 9999)

###########adonis for raup --genus level
set.seed(111)
raup_adonis_post_genus = adonis (raup_dist_post_genus ~ AH, data = env.data.post, permutations = 99)
raup_betadisper_post_genus = betadisper(raup_dist_post_genus, env.data.post$AH, type = "centroid")
set.seed(111)
raup_permutest_post_genus = permutest(raup_betadisper_post_genus, permutations = 9999)

###########adonis for raup --phylum level
set.seed(111)
raup_adonis_post_phylum = adonis (raup_dist_post_phylum ~ AH, data = env.data.post, permutations = 99)
raup_betadisper_post_phylum = betadisper(raup_dist_post_phylum, env.data.post$AH, type = "centroid")
set.seed(111)
raup_permutest_post_phylum = permutest(raup_betadisper_post_phylum, permutations = 9999)
#############################################################################################
#PCOA
#weighted unifrac
#ordinate weighted unifrac
wunifrac.pcoa<-cmdscale(wuni_dist_sp_post, eig = TRUE, k=4)
PCOA.loadings.4axes.wunifrac= wunifrac.pcoa$points[,c(1:4)]
colnames(PCOA.loadings.4axes.wunifrac) <- c("Wuni.PCoA.1", "Wuni.PCoA.2","Wuni.PCoA.3","Wuni.PCoA.4")
PCOA.loadings.4axes.wunifrac = data.frame(PCOA.loadings.4axes.wunifrac)


#Raupcrick
#ordinate 
raup.pcoa<-cmdscale(raup_dist_post_sp, eig = TRUE, k=4)
PCOA.loadings.4axes.raup= raup.pcoa$points[,c(1:4)]
colnames(PCOA.loadings.4axes.raup) <- c("Raup.PCoA.1", "Raup.PCoA.2","Raup.PCoA.3","Raup.PCoA.4")
PCOA.loadings.4axes.raup = data.frame(PCOA.loadings.4axes.raup)

#############################################################################################
##Alpha diversity
#############################################################################################
#OTU richness
OTU_richness = specnumber(votu_post)
OTU.richness = data.frame(OTU_richness)
OTU.richness$AH = env.data$AH[match(row.names(OTU.richness), row.names(env.data))]
range(OTU.richness$OTU_richness) #16 242
AHY = OTU.richness %>% filter(AH == "Y")
range(AHY$OTU_richness) #16 231
AHN = OTU.richness %>% filter(AH =="N")
range(AHN$OTU_richness) #18 242

#ShannonIndex
Shannon_Index = diversity(votu_post, index = "shannon")
#SimpsonIndex 
Simpson = diversity(votu_post, index = "simpson")
#Faith's phylogenetic diversity
Faith_PD = ses.pd(votu_post, phy_tree(post_all), null.model = "richness", runs = 999, iterations = 1000)
obs_faith_pd_mean = Faith_PD$pd.rand.mean

##add the alpha diversity information to other sample data 
bird_data = sample_data(post_all)
bird_data$OTU_richness = OTU_richness
bird_data$Shannon_Index = Shannon_Index
bird_data$Simpson_Index = Simpson
bird_data$Faith_PD = obs_faith_pd_mean

###compare alpha diversity with student t test
t.test(OTU_richness ~ AH, alternative = "greater", var.equal=TRUE,data = data.frame(bird_data))
# data:  OTU_richness by AH
# t = 0.61023, df = 26, p-value = 0.2735
t.test(Shannon_Index ~ AH, alternative = "greater", var.equal=TRUE,data = data.frame(bird_data))
# data:  Shannon_Index by AH
# t = 0.32179, df = 26, p-value = 0.3751
t.test(Simpson_Index ~ AH, alternative = "greater",var.equal=TRUE, data = data.frame(bird_data))
# data:  Simpson_Index by AH
# t = 0.14439, df = 26, p-value = 0.4432
t.test(Faith_PD ~ AH, alternative = "greater", var.equal=TRUE, data = data.frame(bird_data))
# data:  Faith_PD by AH
# t = 0.80481, df = 26, p-value = 0.2141

#add the PCOA axis loadings from the weighted unifrac and raupcrip pcoas 
bird_data$Raup_PCOA1 = PCOA.loadings.4axes.raup$Raup.PCoA.1
bird_data$Raup_PCOA2 = PCOA.loadings.4axes.raup$Raup.PCoA.2
bird_data$Raup_PCOA3 = PCOA.loadings.4axes.raup$Raup.PCoA.3
bird_data$Raup_PCOA4 = PCOA.loadings.4axes.raup$Raup.PCoA.4

bird_data$Wuni_PCOA1 = PCOA.loadings.4axes.wunifrac$Wuni.PCoA.1
bird_data$Wuni_PCOA2 = PCOA.loadings.4axes.wunifrac$Wuni.PCoA.2
bird_data$Wuni_PCOA3 = PCOA.loadings.4axes.wunifrac$Wuni.PCoA.3
bird_data$Wuni_PCOA4 = PCOA.loadings.4axes.wunifrac$Wuni.PCoA.4

##clean up the dataframe so we just see the variables we want to model
clean_bird_data = select(bird_data, X.SampleID, AH, LPS, Condition, 
                         Fever, Activity, OTU_richness, Shannon_Index, Simpson_Index, 
                         Faith_PD, Raup_PCOA1,Raup_PCOA2,Raup_PCOA3,Raup_PCOA4,Wuni_PCOA1,Wuni_PCOA2,
                         Wuni_PCOA3,Wuni_PCOA4)

clean_bird_data$AH = factor(clean_bird_data$AH)



#####################################################################################################################
#SEM
#####################################################################################################################
##FUNCTIONS
#AIC TABLE
library(AICcmodavg)

AICc.lavaan<-function(object, second.ord=TRUE, c.hat = 1, return.K = FALSE){
  object <- as.list(fitMeasures(object))
  npar<-object$baseline.df - object$df
  if(return.K==T) return(object$npar)
  if(second.ord==F && c.hat>1) return(-2*object$logl/c.hat+2*npar)
  if(second.ord==F) return(object$aic)
  if(c.hat>1) return( -2*object$logl/c.hat+2*npar + 2*( npar*(object$npar+1))/(object$ntotal-npar-1))
  object$aic + 2*( npar*(npar+1))/(object$ntotal-npar-1)
}

aictab.lavaan<-function(cand.set, modnames, sort = TRUE, c.hat = 1, second.ord = TRUE, nobs = NULL){
  if(is.null(modnames)) modnames<-1:length(cand.set)
  # check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
  # if (length(unique(check.resp)) > 1) 
  #     stop("You must use the same response variable for all models\n")
  Results <- NULL
  Results <- data.frame(Modnames = modnames)
  Results$K <- unlist(lapply(X = cand.set, FUN = AICc.lavaan, 
                             return.K = TRUE, c.hat = c.hat,second.ord = second.ord))
  Results$AICc <- unlist(lapply(X = cand.set, FUN = AICc.lavaan, 
                                return.K = FALSE, c.hat = c.hat,second.ord = second.ord))
  Results$Delta_AICc <- Results$AICc - min(Results$AICc)
  Results$ModelLik <- exp(-0.5 * Results$Delta_AICc)
  Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)
  if (length(unique(Results$AICc)) != length(cand.set)) 
    warning("\nCheck model structure carefully as some models may be redundant\n")
  if (second.ord == TRUE && c.hat == 1) {
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))
  }
  if (second.ord == TRUE && c.hat > 1) {
    colnames(Results) <- c("Modnames", "K", "QAICc", "Delta QAICc", 
                           "ModelLik", "QAICcWt")
    LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))
    Results$Quasi.LL <- LL/c.hat
    Results$c_hat <- c.hat
  }
  if (second.ord == FALSE && c.hat == 1) {
    colnames(Results) <- c("Modnames", "K", "AIC", "Delta AIC", 
                           "ModelLik", "AICWt")
    Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))
  }
  if (second.ord == FALSE && c.hat > 1) {
    colnames(Results) <- c("Modnames", "K", "QAIC", "Delta QAIC", 
                           "ModelLik", "QAICWt")
    LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))
    Results$Quasi.LL <- LL/c.hat
    Results$c_hat <- c.hat
  }
  if (sort) {
    Results <- Results[rev(order(Results[, 6])), ]
    Results$Cum.Wt <- cumsum(Results[, 6])
  }
  else {
    Results$Cum.Wt <- NULL
  }
  class(Results) <- c("aictab", "data.frame")
  return(Results)
  
}

clean_bird_data$OTU74_rel_abun = rel_abun_votu_post[,"OTU74"]
clean_bird_data$OTU135_rel_abun = rel_abun_votu_post[,"OTU135"]
clean_bird_data$OTU22_rel_abun = rel_abun_votu_post[,"OTU22"]
clean_bird_data$OTU49_rel_abun = rel_abun_votu_post[,"OTU49"]
clean_bird_data$OTU169_rel_abun = rel_abun_votu_post[,"OTU169"]

clean_bird_data$OTU74_PA = PA_votu_post_sp[,"OTU74"]
clean_bird_data$OTU135_PA = PA_votu_post_sp[,"OTU135"]
clean_bird_data$OTU22_PA = PA_votu_post_sp[,"OTU22"]
clean_bird_data$OTU49_PA = PA_votu_post_sp[,"OTU49"]
clean_bird_data$OTU169_PA = PA_votu_post_sp[,"OTU169"]

#remove the rows that have missing data
row.names(clean_bird_data)
lavaan_data = clean_bird_data[-c(13,18,21,25),]#this left df with 24 samples
lavaan_data$Condition = as.character(lavaan_data$Condition)
lavaan_data$Condition = as.numeric(lavaan_data$Condition)
lavaan_data$Fever = as.character(lavaan_data$Fever)
lavaan_data$Fever = as.numeric(lavaan_data$Fever)
lavaan_data$Activity = as.character(lavaan_data$Activity)
lavaan_data$Activity = as.numeric(lavaan_data$Activity)
lavaan_data$OTU_richness = as.character(lavaan_data$OTU_richness)
lavaan_data$OTU_richness = as.numeric(lavaan_data$OTU_richness)

#standardize the data so that the variance isn't so whacky
std_lavaan_data = as.data.frame(scale(lavaan_data[,c(4:23)]))
head(std_lavaan_data)
##put the treatment data and ID data back in 
std_lavaan_data$ID = lavaan_data$X.SampleID
std_lavaan_data$AH = lavaan_data$AH
std_lavaan_data$LPS = lavaan_data$LPS
std_lavaan_data$OTU74_PA = lavaan_data$OTU74_PA
std_lavaan_data$OTU135_PA = lavaan_data$OTU135_PA
std_lavaan_data$OTU22_PA = lavaan_data$OTU22_PA
std_lavaan_data$OTU49_PA = lavaan_data$OTU49_PA
std_lavaan_data$OTU169_PA = lavaan_data$OTU169_PA


#convert the factors AH and LPS into 1s and 0s (Y = 1 N = 0)
std_lavaan_data$AH<-as.factor(ifelse(std_lavaan_data$AH=="Y",1,0))
std_lavaan_data$LPS<-as.factor(ifelse(std_lavaan_data$LPS=="Y",1,0))
std_lavaan_data$OTU74_PA = as.numeric(std_lavaan_data$OTU74_PA)
std_lavaan_data$OTU135_PA =as.numeric(std_lavaan_data$OTU135_PA)
std_lavaan_data$OTU22_PA =  as.numeric(std_lavaan_data$OTU22_PA)
std_lavaan_data$OTU49_PA =  as.numeric(std_lavaan_data$OTU49_PA)
std_lavaan_data$OTU169_PA =  as.numeric(std_lavaan_data$OTU169_PA)

##############################################################################################
#                             FEVER 
##############################################################################################

model1<-'
Fever ~ AH + LPS + Condition + Wuni_PCOA1
'
fit1<-sem(model1,data= std_lavaan_data)
summary(fit1)
AICc.lavaan(fit1)

model2<-'
Fever ~ AH + LPS + Condition + Faith_PD'
fit2<-sem(model2,data= std_lavaan_data)
summary(fit2)
AICc.lavaan(fit2) 


model3<-'Fever ~ AH +LPS+ Condition + OTU_richness'
fit3<-sem(model3,data= std_lavaan_data)
summary(fit3)
AICc.lavaan(fit3) 


model4<-'
Fever~ AH + LPS + Condition + OTU74_rel_abun'
fit4<-sem(model4,data= std_lavaan_data)
summary(fit4)
AICc.lavaan(fit4) 


model5<-'
Wuni_PCOA1 ~ AH
Fever ~ Condition + AH + LPS + Wuni_PCOA1'
fit5<-sem(model5,data= std_lavaan_data)
summary(fit5)
AICc.lavaan(fit5) 

model6<-'
Faith_PD ~ AH
Fever ~ Condition + AH + LPS + Faith_PD'
fit6<-sem(model6,data= std_lavaan_data)
summary(fit6)
AICc.lavaan(fit6) 

model7 <-'
OTU_richness ~ AH
Fever ~ Condition + AH + LPS + OTU_richness'
fit7<-sem(model7,data = std_lavaan_data)
summary(fit7)
AICc.lavaan(fit7) 

model8 <-'
OTU74_rel_abun ~ AH
Fever ~ Condition + AH + LPS + OTU74_rel_abun'

fit8<-sem(model8, data = std_lavaan_data)
summary(fit8)
AICc.lavaan(fit8) 

model9 <-'
Condition ~ AH
Fever ~ AH + LPS + Condition
'
fit9<-sem(model9, data = std_lavaan_data)
summary(fit9)
AICc.lavaan(fit9)

model10 <-'
Fever ~ AH + LPS + Condition
'
fit10<-sem(model10, data = std_lavaan_data)
summary(fit10)
# Regressions:
#   Estimate  Std.Err  z-value  P(>|z|)
# Fever ~                                             
#   AH                0.904    0.353    2.560    0.010
# LPS               0.447    0.344    1.301    0.193
# Condition         0.395    0.180    2.197    0.028
AICc.lavaan(fit10)


aictab.lavaan(list(fit1, fit2,  fit3,fit4, fit5, fit6, fit7,fit8, fit9, fit10), 
              c('fit1','fit2','fit3', 'fit4', 'fit5', 'fit6','fit7', 'fit8', 'fit9','fit10'))

# Model selection based on AICc:
#   
# K   AICc Delta_AICc AICcWt Cum.Wt     LL
# fit10 4  67.40       0.00   0.45   0.45 -29.10
# fit4  5  69.25       1.86   0.18   0.63 -28.57
# fit1  5  69.52       2.12   0.16   0.79 -28.71
# fit3  5  70.29       2.90   0.11   0.89 -29.09
# fit2  5  70.30       2.90   0.11   1.00 -29.10
# fit8  7 135.40      68.01   0.00   1.00 -59.03
# fit9  6 136.65      69.26   0.00   1.00 -61.27
# fit5  7 141.47      74.07   0.00   1.00 -62.07
# fit6  7 142.01      74.62   0.00   1.00 -62.34
# fit7  7 142.30      74.91   0.00   1.00 -62.48

##############################################################################################
#                            ACTIVITY 
##############################################################################################

modelA<-'

Activity ~ AH + LPS + Condition + Wuni_PCOA1
'
fitA<-sem(modelA,data= std_lavaan_data)
summary(fitA)
AICc.lavaan(fitA)

modelB<-'
Activity ~ AH + LPS + Condition + Faith_PD'
fitB<-sem(modelB,data= std_lavaan_data)
summary(fitB)
AICc.lavaan(fitB) 


modelC<-'
Activity ~ AH + LPS + Condition + OTU74_rel_abun '

fitC<-sem(modelC,data= std_lavaan_data)
summary(fitC)
AICc.lavaan(fitC) 


modelD<-'
Activity ~ AH +LPS + Condition+OTU_richness
'
fitD<-sem(modelD,data= std_lavaan_data)
summary(fitD)
AICc.lavaan(fitD) 


modelE<-'
Wuni_PCOA1 ~ AH
Activity ~ Condition + AH + LPS +Wuni_PCOA1
'
fitE<-sem(modelE,data= std_lavaan_data)
summary(fitE)
AICc.lavaan(fitE) 

modelF<-'
Faith_PD ~ AH
Activity ~ Condition + AH + LPS + Faith_PD
'
fitF<-sem(modelF,data= std_lavaan_data)
summary(fitF)
AICc.lavaan(fitF) 
#63.44109

modelG <-'
OTU74_rel_abun ~ AH
Activity ~  AH + LPS + Condition +OTU74_rel_abun '
fitG<-sem(modelG,data = std_lavaan_data)
summary(fitG)
AICc.lavaan(fitG) 


modelH <-'
OTU_richness ~ AH 
Activity ~ Condition + LPS + AH + OTU_richness
'
fitH<-sem(modelH, data = std_lavaan_data)
summary(fitH)
AICc.lavaan(fitH) 

modelI <-'
Condition ~ AH
Activity ~ AH + LPS + Condition'
fitI<-sem(modelI, data = std_lavaan_data)
summary(fitI)
AICc.lavaan(fitI)

modelJ <-'
Activity ~ AH + LPS + Condition'
fitJ<-sem(modelJ, data = std_lavaan_data)
summary(fitJ)
# Regressions:
#   Estimate  Std.Err  z-value  P(>|z|)
# Activity ~                                          
# AH                0.038    0.314    0.120    0.904
# LPS              -1.260    0.306   -4.121    0.000
# Condition        -0.223    0.160   -1.393    0.164
AICc.lavaan(fitJ)


aictab.lavaan(list(fitA, fitB, fitC,fitD,  fitE, fitF, fitG,fitH,  fitI, fitJ), 
              c('fitA','fitB', 'fitC','fitD',  'fitE', 'fitF', 'fitG', 'fitH', 'fitI','fitJ'))


# 
# Model selection based on AICc:
#     K   AICc Delta_AICc AICcWt Cum.Wt     LL
# fitJ 4  61.79       0.00   0.36   0.36 -26.30
# fitC 5  62.18       0.38   0.30   0.66 -25.04
# fitD 5  63.42       1.62   0.16   0.82 -25.66
# fitA 5  64.37       2.58   0.10   0.92 -26.13
# fitB 5  64.68       2.89   0.08   1.00 -26.29
# fitG 7 128.33      66.53   0.00   1.00 -55.50
# fitI 6 131.05      69.26   0.00   1.00 -58.47
# fitH 7 135.43      73.63   0.00   1.00 -59.05
# fitE 7 136.32      74.53   0.00   1.00 -59.50
# fitF 7 136.40      74.60   0.00   1.00 -59.53





################
#supplemental analyses including egg detection data
################

##Get a phyloseq object that just has the samples with data about egg presence/absence in it
##this means getting rid of samples that DON'T have this information, which we will want to un-do for
##other analyses that have nothing to do with whether an egg was detected
##so just remember, we want to use the dataset that include samples '06','12','35','57','61','65',and '69' later
all_egg_data = subset_samples(post_all, parasite_egg_ever != "#N/A" )

#remove any taxa that were unique to the sample we just removed so we don't have taxa with 0 reads in our data
all_egg_data = prune_taxa(taxa_sums(all_egg_data)>0, all_egg_data) #809 taxa and 21 samples

env.data_all_egg = data.frame(sample_data(all_egg_data))

#get vegan OTU table (raw abundance)
votu_egg <- veganotu(all_egg_data)
votu_egg = as.data.frame(votu_egg)

#get vegan OTU table (presence-absence)
PA_votu_egg= votu_egg
PA_votu_egg[PA_votu_egg>0] <-1
PA_votu_egg = as.data.frame(PA_votu_egg)

###zoom out to the level of genera
egg_genera = tax_glom(all_egg_data, taxrank = rank_names(all_egg_data)[6])
tax_table(egg_genera)
egg_genera_votu<- veganotu(egg_genera)
egg_genera_votu = as.data.frame(egg_genera_votu)

#get PA vegan table for genera
PA_votu_genera_egg = egg_genera_votu
PA_votu_genera_egg[PA_votu_genera_egg>0]<-1
PA_votu_genera_egg = as.data.frame(PA_votu_genera_egg)


###zoom out to the level of phylum
egg_phylum = tax_glom(all_egg_data, taxrank = rank_names(all_egg_data)[2])
egg_phylum_votu<- veganotu(egg_phylum)
egg_phylum_votu = as.data.frame(egg_phylum_votu)

#get PA vegan table for phylum
PA_votu_phyla_egg = egg_phylum_votu
PA_votu_phyla_egg[PA_votu_phyla_egg>0]<-1
PA_votu_phyla_egg = as.data.frame(PA_votu_phyla_egg)
#########
#Indicator species
#get a relative abundance otu table 
#for genus
egg_rel_abun_genus = transform_sample_counts(egg_genera, function(x) x / sum(x))
rel_abun_egg_genera_votu<-veganotu(egg_rel_abun_genus)
rel_abun_egg_genera_votu = as.data.frame(rel_abun_egg_genera_votu)

rel_egg_indic_genus=multipatt(rel_abun_egg_genera_votu, env.data_all_egg[["parasite_egg_ever"]], func = "IndVal.g", duleg = TRUE, set.seed(1),control = how(nperm = 999))
summary(rel_egg_indic_genus)
tax_table(all_egg_data)[c("OTU402","OTU1577","OTU834", "OTU3413")]

######
##BETA DIVERSITY--ORDINATION
###############################################
#weighted unifrac
egg_wuni_dist<-phyloseq::distance(all_egg_data, method = "wunifrac")

#ordinate weighted unifrac
egg.wunifrac.nmds<-metaMDS(egg_wuni_dist, distance = "wunifrac",autotransform = FALSE, k=2, trymax=50)

#weighted unifrac at genus level
egg_genera_wuni_dist<-phyloseq::distance(egg_genera, method = "wunifrac")
#ordinate weighted unifrac
genera_grace.egg.wunifrac.nmds<-metaMDS(egg_genera_wuni_dist, distance = "wunifrac", autotransform = FALSE,k=2, trymax=50)


#weighted unifrac at phyla level
egg_phyla_wuni_dist<-phyloseq::distance(egg_phylum, method = "wunifrac")
#ordinate weighted unifrac
phyla_grace.egg.wunifrac.nmds<-metaMDS(egg_phyla_wuni_dist, distance = "wunifrac", autotransform = FALSE,k=2, trymax=50)

#############################################################################################
#Raupcrick
egg_raup_dist<-vegdist(PA_votu_egg, method = "raup")
#ordinate using presence absence data
egg.raup.nmds<-metaMDS(egg_raup_dist, dist = "raup",autotransform = FALSE, k=2, trymax = 200)

#############################################################################################
#raup at genus level
egg_raup_genera_dist<-raupcrick(PA_votu_genera_egg, null = "r1", nsimul = 999, chase = FALSE)
#ordinate using presence absence data
egg.raup.genera.nmds<-metaMDS(egg_raup_genera_dist, dist = "raup", autotransform = FALSE, k=3, trymax = 200)
#############################################################################################
#raup at phylum level
egg_raup_phyla_dist<-raupcrick(PA_votu_phyla_egg, null = "r1", nsimul = 999, chase = FALSE)
#ordinate using presence absence data
egg.raup.phyla.nmds<-metaMDS(egg_raup_phyla_dist, dist = "raup", autotransform = FALSE, k=3, trymax = 200)

###PERMANOVA and Betadispersion tests
###########adonis for weighted unifrac --species level
set.seed(111)
w_uni_adonis_egg_sp= adonis(egg_wuni_dist ~ parasite_egg_ever, data = env.data_all_egg, permutations = 99)
w_uni_betadisper_egg_sp = betadisper(egg_wuni_dist, env.data_all_egg$parasite_egg_ever,type = "centroid")
set.seed(111)
w_permutest_egg_sp = permutest(w_uni_betadisper_egg_sp, permutations = 9999) 
###########adonis for weighted unifrac --genus level
set.seed(111)
w_uni_adonis_egg_genus= adonis(egg_genera_wuni_dist ~ parasite_egg_ever, data = env.data_all_egg, permutations = 99)
w_uni_betadisper_egg_genus = betadisper(egg_genera_wuni_dist, env.data_all_egg$parasite_egg_ever,type = "centroid")
set.seed(111)
w_permutest_egg_genus = permutest(w_uni_betadisper_egg_genus, permutations = 9999) 
###########adonis for weighted unifrac --phylum level
set.seed(111)
w_uni_adonis_egg_phylum= adonis(egg_phyla_wuni_dist ~ parasite_egg_ever, data = env.data_all_egg, permutations = 99)
w_uni_betadisper_egg_phylum = betadisper(egg_phyla_wuni_dist, env.data_all_egg$parasite_egg_ever,type = "centroid")
set.seed(111)
w_permutest_egg_phylum = permutest(w_uni_betadisper_egg_phylum, permutations = 9999) 


###########adonis for raup --species level
set.seed(111)
raup_adonis_egg_sp = adonis(egg_raup_dist ~ parasite_egg_ever, data = env.data_all_egg, permutations = 99)
raup_betadisper_egg_sp = betadisper(egg_raup_dist, env.data_all_egg$parasite_egg_ever, type = "centroid")
set.seed(111)
raup_permutest_egg_sp = permutest(raup_betadisper_egg_sp, permutations = 9999)

###########adonis for raup --genus level
set.seed(111)
raup_adonis_egg_genus = adonis (egg_raup_genera_dist ~ parasite_egg_ever, data = env.data_all_egg, permutations = 99)
raup_betadisper_egg_genus = betadisper(egg_raup_genera_dist, env.data_all_egg$parasite_egg_ever, type = "centroid")
set.seed(111)
raup_permutest_egg_genus = permutest(raup_betadisper_egg_genus, permutations = 9999)

###########adonis for raup --phylum level
set.seed(111)
raup_adonis_egg_phylum = adonis (egg_raup_phyla_dist ~ parasite_egg_ever, data = env.data_all_egg, permutations = 99)
raup_betadisper_egg_phylum = betadisper(egg_raup_phyla_dist, env.data_all_egg$parasite_egg_ever, type = "centroid")
set.seed(111)
raup_permutest_egg_phylum = permutest(raup_betadisper_egg_phylum, permutations = 9999)


#############################################################################################
#plot

tiff("FigureS2.tiff", units="in", width=8.5, height=3, res=300)
par(mfrow= c(1,2), las = 1)

plot(egg.wunifrac.nmds$points[,1:2], type = "n",xlab="NMDS Axis 1",ylab="NMDS Axis 2")
#text(egg.wunifrac.nmds$points[,1:2], labels=env.data_all_egg$X.SampleID)
points(egg.wunifrac.nmds$points[env.data_all_egg$parasite_egg_ever=="N",1],egg.wunifrac.nmds$points[env.data_all_egg$parasite_egg_ever=="N",2], pch= 21, cex = 1.25, col = 'black',bg='grey60', lwd = 2)
points(egg.wunifrac.nmds$points[env.data_all_egg$parasite_egg_ever=="Y",1],egg.wunifrac.nmds$points[env.data_all_egg$parasite_egg_ever=="Y",2], pch= 21, cex = 1.25, col = 'black',bg='purple', lwd = 2)
title(main = "Weighted Unifrac")
legend("bottomright", c("Capillariid-eggs\n detected","No eggs"), xpd = TRUE, inset = c(-0.45,0),
       bty = "n", pch = c(21, 21), pt.bg = c("purple","grey60"),col = c("black","black"),pt.lwd = 2,pt.cex = 1.4, cex = 0.75)


#plot
plot(egg.raup.nmds$points[,1:2], type = "n", xlab = "NMDS Axis 1", ylab = "NMDS Axis 2")
points(egg.raup.nmds$points[env.data_all_egg$parasite_egg_ever=="Y",1],egg.raup.nmds$points[env.data_all_egg$parasite_egg_ever=="Y",2], pch= 21, cex = 1.25, col = 'black',bg='purple', lwd = 2)
points(egg.raup.nmds$points[env.data_all_egg$parasite_egg_ever=="N",1],egg.raup.nmds$points[env.data_all_egg$parasite_egg_ever=="N",2], pch= 21, cex = 1.25, col = 'black',bg='grey60', lwd = 2)
title(main = "Raup-Crick")

dev.off()
##################
###alpha diversity
#match egg detection data into the bird_data sample dataframe (where all the alpha div metrics are)
egg_sample_data = sample_data(all_egg_data)
bird_data$egg_det = egg_sample_data$parasite_egg_ever[match(row.names(bird_data),row.names(egg_sample_data))]
#cut this dataframe down to only include rows with data for parasite egg detection
bird_data_w_egg_data = bird_data %>%
  filter(egg_det != "<NA>")

###compare alpha diversity with student t test
t.test(OTU_richness ~ egg_det, alternative = "greater", var.equal=TRUE,data = data.frame(bird_data_w_egg_data))
# data:  OTU_richness by egg_det
# t = 0.49517, df = 19, p-value = 0.3131
t.test(Shannon_Index ~ AH, alternative = "greater", var.equal=TRUE,data = data.frame(bird_data_w_egg_data))
# data:  Shannon_Index by AH
# t = 0.18197, df = 19, p-value = 0.4288
t.test(Simpson_Index ~ AH, alternative = "greater",var.equal=TRUE, data = data.frame(bird_data_w_egg_data))
# data:  Simpson_Index by AH
# t = 0.045034, df = 19, p-value = 0.4823
t.test(Faith_PD ~ AH, alternative = "greater", var.equal=TRUE, data = data.frame(bird_data_w_egg_data))
# data:  Faith_PD by AH
# t = 0.71568, df = 19, p-value = 0.2414


##############
#look at fecal float data to see whether egg detection changed after anthelmintic treatment

float<-read.csv("float_data.csv")
head(float)
float$parasite_det<-ifelse(float$parasite == "npd", 0,1)
float$pre_detect<-ifelse(float$pre.post == "pre" & float$parasite_det == 1, 1, 0)
sum(float$pre_detect)
length(float$pre_detect) #14
float$pre = ifelse(float$pre.post == "pre", 1,0)
sum(float$pre) #52

#proportion of pre-treatment samples in which parasite eggs detected = 14/52 = 0.2692308

float$post_AH <-ifelse(float$pre.post == "post" & float$Anti.helm..== "Y", 1, 0)
float$post_ah_detect <- ifelse(float$post_AH == 1 & float$parasite_det == 1, 1, 0)

library(binom)
binom.confint(x = c(2, 4), n = 100, tol = 1e-8)

head(float)  
length(float$bird_id)  
binom.confint(22, 80,conf.level = 0.95,methods = "wilson")  

sum(float$parasite_det)
length(float$bird_id)
length(float$parasite_det[which(float$pre.post=="pre" & float$Anti.helm..=="N")])
length(float$parasite_det[which(float$pre.post=="pre" & float$Anti.helm..=="Y")])
length(float$parasite_det[which(float$pre.post=="post" & float$Anti.helm..=="N")])
length(float$parasite_det[which(float$pre.post=="post" & float$Anti.helm..=="Y")])

length(float$parasite_det[which(float$pre.post=="pre")])

binom.confint(sum(float$parasite_det[which(float$pre.post=="pre")]),52,conf.level = 0.95,methods = "wilson")


binom.confint(sum(float$parasite_det[which(float$pre.post=="pre" & float$Anti.helm..=="N")]), length(float$parasite_det[which(float$pre.post=="pre" & float$Anti.helm..=="N")]),
              conf.level = 0.95,methods = "wilson")  
binom.confint(sum(float$parasite_det[which(float$pre.post=="pre" & float$Anti.helm..=="Y")]), length(float$parasite_det[which(float$pre.post=="pre" & float$Anti.helm..=="Y")]),
              conf.level = 0.95,methods = "wilson")  
binom.confint(sum(float$parasite_det[which(float$pre.post=="post" & float$Anti.helm..=="N")]), length(float$parasite_det[which(float$pre.post=="post" & float$Anti.helm..=="N")]),
              conf.level = 0.95,methods = "wilson")  
binom.confint(sum(float$parasite_det[which(float$pre.post=="post" & float$Anti.helm..=="Y")]), length(float$parasite_det[which(float$pre.post=="post" & float$Anti.helm..=="Y")]),
              conf.level = 0.95,methods = "wilson")  
pp<-as.vector(c("pre", "post","post"))
trt<-as.vector(c("NA","N","Y"))
mean<-as.vector(c(0.2692308, 0.2857143,0.2857143))
lower<-as.vector(c(0.1676896, 0.1172138, 0.1172138))
upper<-as.vector(c(0.4025222, 0.5464908, 0.5464908))
n<-as.vector(c("52","14","14"))
fecal_egg<-data.frame(pp,trt,mean,lower,upper,n)
head(fecal_egg)
str(fecal_egg)

fecal_egg$pp<-relevel(fecal_egg$pp, "pre")
fecal_egg$pp<-factor(fecal_egg$pp,labels=c("Pre-treatment", "Post-treatment"))
fecal_egg$trt<-factor(fecal_egg$trt,labels=c("Control", "Pre-treatment","Anthelminthic"))
fecal_egg$trt<-relevel(fecal_egg$trt, "Pre-treatment")

limits<-aes(ymax = fecal_egg$upper, ymin = fecal_egg$lower)

length(unique(float$bird_id))

n_fun<-function(x){
  return(data.frame(y=mean(x),label=paste0("n=",length(x))))}

p<-ggplot(data = fecal_egg, aes(x=factor(pp),y=mean, group = trt, shape = factor(trt)))


p + geom_point (position = position_dodge(0.9),size = 4)+
  geom_errorbar(limits,position=position_dodge(0.9), 
                width = 0.25)+
  labs(x = "Pre or Post treatment", y = "Proportion of samples with parasite eggs detected")+
  theme(panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = c(0.15,0.85),
        legend.background = element_blank(),
        axis.title =  element_text(size=14))+
  annotate("text", x = "Pre-treatment",y=0.413, label = "52")+
  annotate("text", x = 1.771,y=0.556, label = "14")+
  annotate("text", x = 2.221,y=0.556, label = "14")

