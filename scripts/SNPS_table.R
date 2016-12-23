
####R version: 3.2.2 (2015-08-14) -- "Fire Safety"####
####Copyright (C) 2015 The R Foundation for Statistical Computing###
####Platform: x86_64-pc-linux-gnu (64-bit)###
####Authors: Jan Hattendorf & Daniela Brites###
####Date:1.12.2016####
####content: This script performs extracts allele frequencies and performs some statistics based on the SFS for Mtb/HIVpos and Mtb/HIVneg####

#read snps table; this table contains all position after filtering the annotation for IS elements, phages,transposases,PE/PPE/PGRS.
table_snps <- read.table("~/data/table.HIV_UGII.genotypes",na.strings="",stringsAsFactors=F,sep="\t",header=T)

dim(table_snps)

#read annotation
annotation <- read.table("~/data/HIV_UgandaII.annovar_in.finalout_updated",na.strings="",stringsAsFactors=F,sep="\t",header=F,fill=T)

dim(annotation)

colnames(annotation) <- c("position","ref","mut","NS/S/I","Rv","in_vitro","macrophage","tuberculist_cat","size_aa","DR")

#put together table of snps plus annotation
table_snps_annotation <-cbind(table_snps,annotation)

###get table with drug resistance mutations based on Walker 2015
DR_mutations <- table_snps_annotation$DR=="DR_position"
to.select_DR <- which(DR_mutations)
table_DR <- table_snps_annotation[c(to.select_DR),]
dim(table_DR)

##remove drug resistance positions from the table
table_snps_annotation_filtered<- table_snps_annotation[-c(to.select_DR),]
dim(table_snps_annotation_filtered)


# Remove recent duplication within the ESX family (J (Rv1038c), W (Rv3620c), K (Rv1197), P (Rv2347c), M (Rv1792), I (Rv1037c), V (Rv3619c), N (Rv1793), L (Rv1198), O (Rv2346c), G (Rv0287), S (Rv3020c), H (Rv0288), R (Rv3019c), Q (Rv3017c) )

esx_duplicated <- c("Rv1038c", "Rv3620c", "Rv1197", "Rv2347c","Rv1792","Rv1037c", "Rv3619c", "Rv1793", "Rv1198", "Rv2346c", "Rv0287", "Rv3020c", "Rv0288","Rv3019c", "Rv3017c")
no_esx_duplicated <- apply(table_snps_annotation_filtered, 1, function(x) !any(x %in% esx_duplicated ))
table_snps_annotation_filtered <- table_snps_annotation_filtered[c(no_esx_duplicated),]
dim(table_snps_annotation_filtered)


##Remove genes wich have mutations associated with drug resistances (the position of rrs have been removed alredy) 
DR_genes <- c("Rv2428", "Rv1483", "Rv1484", "Rv1908c","Rv0667","Rv3795","Rv2043c", "Rv1630", "Rv0006", "Rv0682","Rv3919c","Rv1694","Rv2416c")
no_DR_genes <- apply(table_snps_annotation_filtered, 1, function(x) !any(x %in% DR_genes))
table_snps_annotation_filtered <- table_snps_annotation_filtered[c(no_DR_genes),]
dim(table_snps_annotation_filtered)


##Remove also the genes for which there is identical duplications of streches of 50 bp. 
genes_repetions <- c("Rv1572c","Rv1574","Rv1575","Rv0336","Rv0515","IG1195_Rv1174c-Rv1175c","IG127_Rv0126-Rv0127","IG1711_Rv1682-Rv1683","IG18_Rv0018c-Rv0019c","IG3012_Rv2965c-Rv2966c","IG3013_Rv2966c-Rv2967c","IG533_Rv0525-Rv0526","IG559_Rv0551c-Rv0552","IG622_Rv0612-Rv0613c","IG71_Rv0071-Rv0072","IG784_Rv0769-Rv0770","IG877_Rv0861c-Rv0862c","Rv0031","Rv0094c","Rv0095c","Rv0096","Rv0109","Rv0124","Rv0151c","Rv0152c","Rv0159c","Rv0160c","Rv0256c","Rv0257","Rv0277c","Rv0278c","Rv0279c","Rv0280","Rv0285","Rv0286","Rv0297","Rv0304c","Rv0305c","Rv0335c","Rv0353","Rv0354c","Rv0355c","Rv0387c","Rv0388c","Rv0393","Rv0397","Rv0442c","Rv0453","Rv0487","Rv0490","Rv0532","Rv0538","Rv0578c","Rv0605","Rv0605","Rv0740","Rv0741","Rv0742","Rv0746","Rv0747","Rv0750","Rv0754","Rv0755A","Rv0755c","Rv0795","Rv0796","Rv0797","Rv0814c","Rv0823c","Rv0829","Rv0832","Rv0833","Rv0834c","Rv0850","Rv0867c","Rv0872c","Rv0878c","Rv0915c","Rv0916c","Rv0920c","Rv0921","Rv0922","Rv0977","Rv0978c","Rv0980c","Rv1034c","Rv1035c","Rv1036c","Rv1037c","Rv1038c","Rv1039c","Rv1040c","Rv1041c","Rv1042c","Rv1047","Rv1067c","Rv1068c","Rv1087","Rv1088","Rv1089","Rv1091","Rv1128c","Rv1135c","Rv1148c","Rv1149","Rv1150","Rv1168c","Rv1169c","Rv1172c","Rv1173","Rv1195","Rv1196","Rv1197","Rv1198","Rv1199c","Rv1214c","Rv1243c","Rv1288","Rv1295","Rv1313c","Rv1318c","Rv1319c","Rv1325c","Rv1361c","Rv1369c","Rv1370c","Rv1386","Rv1387","Rv1396c","Rv1430","Rv1441c","Rv1450c","Rv1452c","Rv1458c","Rv1468c","Rv1489A","Rv1493","Rv1548c","Rv1557","Rv1558","Rv1573","Rv1574","Rv1575","Rv1576c","Rv1577c","Rv1578c","Rv1579c","Rv1580c","Rv1581c","Rv1582c","Rv1583c","Rv1584c","Rv1585c","Rv1586c","Rv1587c","Rv1588c","Rv1646","Rv1651c","Rv1702c","Rv1705c","Rv1706c","Rv1753c","Rv1756c","Rv1757","Rv1758","Rv1759c","Rv1763","Rv1764","Rv1765A","Rv1765c","Rv1768","Rv1787","Rv1788","Rv1789","Rv1790","Rv1791","Rv1793","Rv1800","Rv1801","Rv1802","Rv1803c","Rv1806","Rv1807","Rv1808","Rv1809","Rv1818c","Rv1829","Rv1840c","Rv1910c","Rv1911c","Rv1917c","Rv1918c","Rv1945","Rv1983","Rv2013","Rv2014","Rv2015c","Rv2048c","Rv2082","Rv2085","Rv2090","Rv2105","Rv2106","Rv2107","Rv2108","Rv2112c","Rv2123","Rv2126c","Rv2162c","Rv2167c","Rv2168c","Rv2177c","Rv2196","Rv2258c","Rv2277c","Rv2278","Rv2279","Rv2328","Rv2340c","Rv2346c","Rv2347c","Rv2352c","Rv2353c","Rv2354","Rv2355","Rv2356c","Rv2371","Rv2396","Rv2408","Rv2424c","Rv2430c","Rv2431c","Rv2460c","Rv2461c","Rv2479c","Rv2480c","Rv2487c","Rv2489","Rv2490c","Rv2512c","Rv2519","Rv2543","Rv2544","Rv2591","Rv2608","Rv2615c","Rv2634c","Rv2648","Rv2649","Rv2650c","Rv2651c","Rv2652c","Rv2653c","Rv2654c","Rv2655c","Rv2656c","Rv2657c","Rv2659c","Rv2665","Rv2666","Rv2673","Rv2680","Rv2689c","Rv2690c","Rv2741","Rv2768c","Rv2769c","Rv2770c","Rv2774c","Rv2791c","Rv2792c","Rv2805","Rv2807","Rv2810c","Rv2812","Rv2814c","Rv2815c","Rv2825c","Rv2828c","Rv2853","Rv2859c","Rv2882c","Rv2885c","Rv2886c","Rv2892c","Rv2931","Rv2932","Rv2943","Rv2943A","Rv2944","Rv2961","Rv2977c","Rv2978c","Rv2979c","Rv2980","Rv3018A","Rv3018c","Rv3021c","Rv3022A","Rv3022c","Rv3023c","Rv3097c","Rv3115","Rv3125c","Rv3135","Rv3136","Rv3144c","Rv3159c","Rv3184","Rv3185","Rv3186","Rv3187","Rv3191c","Rv3281","Rv3325","Rv3326","Rv3327","Rv3343c","Rv3344c","Rv3345c","Rv3346c","Rv3347c","Rv3348","Rv3349c","Rv3350c","Rv3355c","Rv3367","Rv3380c","Rv3381c","Rv3386","Rv3387","Rv3388","Rv3424c","Rv3425","Rv3426","Rv3427c","Rv3428c","Rv3429","Rv3430c","Rv3431c","Rv3466","Rv3467","Rv3474","Rv3475","Rv3477","Rv3478","Rv3507","Rv3508","Rv3511","Rv3512","Rv3513c","Rv3514","Rv3515c","Rv3532","Rv3533c","Rv3539","Rv3558","Rv3590c","Rv3595c","Rv3611","Rv3619c","Rv3620c","Rv3621c","Rv3622c","Rv3636","Rv3637","Rv3638","Rv3639c","Rv3640c","Rv3650","Rv3680","Rv3710","Rv3738c","Rv3739c","Rv3746c","Rv3798","Rv3812","Rv3826","Rv3827c","Rv3828c","Rv3844","Rv3873","Rv3876","Rv3892c","Rv3893c")
no_genes_repetions <- apply(table_snps_annotation_filtered, 1, function(x) !any(x %in% genes_repetions))
table_snps_annotation_filtered <- table_snps_annotation_filtered[ c(no_genes_repetions),]
dim(table_snps_annotation_filtered)

write.table(table_snps_annotation_filtered,file="table_snps_annotation_filtered",sep="\t",col.names=T)

##The total number of mutations respect to the reference after filtering is 5312 for 180 genomes.


hiv_status <- read.table("~/data/list_name_status",na.strings="",stringsAsFactors=F,sep="\t",header=F,fill=T)
 

#######Analysis of Mtb variants intra-host#####

#Obtain the number of mutations per genome

mutations_per_genome <- ifelse (table_snps_annotation_filtered==table_snps_annotation_filtered$mut,1,0)

mutations_per_genome_table <- colSums(mutations_per_genome[,1:180])
mutations_per_genome_table <- as.data.frame (mutations_per_genome_table, row.names = NULL)

mutations_per_genome_table <- cbind(mutations_per_genome_table,hiv_status)
colnames(mutations_per_genome_table) <- c("nHomoMut","name","HIV","year")


#plot

boxplot(nHomoMut~HIV,data=mutations_per_genome_table,ylab="N fixed mutations per isolate", names=c("Mtb/HIV-","Mtb/HIV+"),border=c("blue","red"))

wilcox.test(nHomoMut~HIV,data=mutations_per_genome_table)

#Transforms all ambiguities in indeterminations
R_per_genome <- ifelse (table_snps_annotation_filtered=="R",1,0)
R_per_genome_table <- colSums(R_per_genome[,1:180])

S_per_genome <- ifelse (table_snps_annotation_filtered=="S",1,0)
S_per_genome_table <- colSums(S_per_genome[,1:180])

W_per_genome <- ifelse (table_snps_annotation_filtered=="W",1,0)
W_per_genome_table <- colSums(W_per_genome[,1:180])

K_per_genome <- ifelse (table_snps_annotation_filtered=="K",1,0)
K_per_genome_table <- colSums(K_per_genome[,1:180])

M_per_genome <- ifelse (table_snps_annotation_filtered=="M",1,0)
M_per_genome_table <- colSums(M_per_genome[,1:180])

indertermination_table <-cbind( R_per_genome_table,S_per_genome_table, W_per_genome_table,K_per_genome_table,M_per_genome_table)
indertermination_table <- as.data.frame (indertermination_table, row.names = NULL)

indertermination_table$nIndeter <- indertermination_table$R_per_genome_table+indertermination_table$S_per_genome_table+indertermination_table$W_per_genome_table+indertermination_table$K_per_genome_table+indertermination_table$M_per_genome_table
mutations_per_genome_table<- cbind(mutations_per_genome_table,indertermination_table$nIndeter)
colnames(mutations_per_genome_table) <- c("nHomoMut","name","HIV","year","nIndeter")

hiv <- factor(mutations_per_genome_table$HIV==1, labels=c("nonHIV","HIV"))
table(mutations_per_genome_table$nIndeter,hiv)

myfun <- function(x){c(mean = mean(x), median = median(x))}
tapply(mutations_per_genome_table$nIndeter,mutations_per_genome_table$HIV,myfun)

boxplot(nIndeter~HIV,data=mutations_per_genome_table,ylab="ivSNPs per genome", names=c("Mtb/HIV-","Mtb/HIV+"),border=c("blue","red"))

wilcox.test(nIndeter~HIV,data=mutations_per_genome_table)


####calculations of SFS

#d_hiv contains all strains (rows) and positons (5312) plus the two last collums are hiv status and year

dim(table_snps_annotation_filtered)

dim(hiv_status)

##remove annoatation from the table of SNPs: d has only SNPs
d <- table_snps_annotation_filtered[,1:180]
dim(d)

#create a table with the information on SNPs and HIV status. 
d_hiv <- cbind(as.data.frame(t(d), stringsAsFactors = F),hiv_status[,2:3])
dim(d_hiv)

colnames(d_hiv)[5313] <- "hiv"
colnames(d_hiv)[5314] <- "year"

d_hiv$hiv <- as.numeric(as.character(d_hiv$hiv))
d_hiv <- d_hiv[order(d_hiv$hiv),]

#How many genomes per HIV status and year. 
table(d_hiv$year,d_hiv$hiv)

#plot the previous table
barplot(table(d_hiv$hiv,d_hiv$year),ylab= "Number of genomes", xlab="Isolation year",beside=T,	legend = table(d_hiv$hiv),names.arg=c("1995","1996","1997","1998","1999","2000","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012"),las=2)

#annoatation per variable position
annotation_filtered <-  as.data.frame(t(table_snps_annotation_filtered[,181:190]))
dim(annotation_filtered)

#define mutation (derived allele) 
mut <- annotation_filtered[3,]

#define HIV.
hiv <- d_hiv$hiv

#remove HIV and year from d_hiv. d will be again just the table of snps but it will be ordered by hiv wich is important for some of the steps bellow.
d <- d_hiv[,-5313:-5314]

##define function to detect invariable sites. In this way invariable sites will be those where the mutation is fixed in all strains, i.e if there is an intermination, or gap, that site will not be consired as fixed.
is.not.variable <- function(x) length(table(x[x %in% c("A","G","C","T","N")]))==1
not.variable <- apply(d,2,is.not.variable)
table(not.variable)



#non-variable 593 sites, remove them from the data set
to.remove0 <- which(not.variable)
d <- d[, -c(to.remove0)]
dim(d)

#remove those from mut  and annotatopn as well.
mut <-  mut[,-c(to.remove0)]

annotation_filtered <- annotation_filtered[, -c(to.remove0)]

dim(annotation_filtered)


## set ambiguities to NA
table(unlist(d))
d[d=="K"] <- NA
d[d=="M"] <- NA
d[d=="R"] <- NA
d[d=="S"] <- NA
d[d=="W"] <- NA
d[d=="Y"] <- NA

####dataset prepared, start with the analysis of the segregating sites in the population####


#create a dataset the same size as d but with the information of mut.
dmut <- mut[rep(1,dim(d)[1]),]
dim(dmut)


# identifying the differences (missing values are not considered as difference). mutg has the count of derived alleles per position
ismut <- (d == dmut)
mutg <- colSums(ismut, na.rm=T)
table(mutg, useNA="always")

#number of singletons
table(mutg==1)

#define factor singleton and non-singleton 
glob.singl <- factor(mutg==1, labels = c("non-singleton","singleton"))
table(glob.singl, useNA="always")

##counts for HIV+ and HIV-
ismuthiv <- (d[hiv==1,] == dmut[hiv==1,])
mut1 <- colSums(ismuthiv, na.rm=T)
table(mut1, useNA="always")

#Singletons in Mtb/HIVpos
table(mut1==1)

#attribute factor to singletons and nonsingletons
hiv.factor <- cut(mut1, breaks=c(0,1,2,10000), labels=c("no-mut","singel","non-sing"), right = F)
table(hiv.factor, useNA="ifany")

#How many Mtb/HIVpos singletons are also "Population singletons".
table(hiv.factor[glob.singl=="singleton"], useNA="always")

#How many Mtb/HIVpos singletons are "Group singletons" 
table(hiv.factor[glob.singl=="non-singleton"], useNA="always")

#Singletons in Mtb/HIVneg 
ismutneg <- (d[hiv==0,] == dmut[hiv==0,])
mut0 <- colSums(ismutneg, na.rm=T)
table(mut0, useNA="always")
table(mut0==1)

#attribute factor to singletons and nonsingletons
neg.factor <- cut(mut0, breaks=c(0,1,2,10000), labels=c("no-mut","singel","non-sing"), right = F)
table(neg.factor, useNA="always")


#How many Mtb/HIVneg singletons are also "Population singletons".
table(neg.factor[glob.singl=="singleton"], useNA="always")  

#How many Mtb/HIVneg singletons are "Group singletons" 
table(neg.factor[glob.singl=="non-singleton"], useNA="always")  

#Summarizing Total Singletons, Population and group singletons for Mtb/HIVpos and for Mtb/HIVneg
table(hiv.factor, glob.singl)
table(neg.factor, glob.singl)

dhiv <- d[hiv==1,]
dim(dhiv)
#[1]   86 4719

dneg <- d[hiv==0,]
dim(dneg)
#[1]   94 4719

dim(d)
dim(dmut)

##remove stop gain and stop lost fom the annotation data.
annotation_filtered[annotation_filtered=="stoplost"]<-NA
annotation_filtered[annotation_filtered=="stopgain"]<-NA


####start with bootstrap for obtaining 95% confidence intervals 

set.seed(110303)

# sample Mtb/hIVneg with n = 94 and Mtb/HiVpos with n= 86
 
loops <- 1000
n0 <- sum(hiv %in% 0)
n1 <- sum(hiv %in% 1)

bor <-  matrix(NA,loops,33+38)
colnames(bor) <-c(
c("prop1","prop0","bn1","bn0","nonsing.i.1","nonsing.ns.1","nonsing.s.1","sing.i.1","sing.ns.1","sing.s.1","ns.i.0","ns.ns.0","ns.s.0","s.i.0","s.ns.0","s.s.0","nonsing.e.1","nonsing.ne.1","sing.e.1","sing.ne.1","ns.e.0","ns.ne.0","s.e.0","s.ne.0","dif")
,
paste(c(paste(c("nonsing","sing"),".ns.",sort(rep(c("e","ne"),2)),sep="") ,paste(c("nonsing","sing"),".s.",sort(rep(c("e","ne"),2)),sep="")),".1",sep="")
,
paste(c(paste(c("nonsing","sing"),".ns.",sort(rep(c("e","ne"),2)),sep="") ,paste(c("nonsing","sing"),".s.",sort(rep(c("e","ne"),2)),sep="")),".0",sep="")
,
paste(c("nseg.sites.1","pr.nonsing.ns.1","pr.sing.ns.1","pr.nonsing.ns.0","pr.sing.ns.0","pr.nonsing.s.1","pr.sing.s.1","pr.nonsing.s.0","pr.sing.s.0","pr.nonsing.i.1","pr.sing.i.1","pr.nonsing.i.0","pr.sing.i.0","pr.e.nonsing.ns.1","pr.e.sing.ns.1","pr.e.nonsing.ns.0","pr.e.sing.ns.0","pr.e.nonsing.s.1","pr.e.sing.s.1","pr.e.nonsing.s.0","pr.e.sing.s.0","pr.ne.nonsing.ns.1","pr.ne.sing.ns.1","pr.ne.nonsing.ns.0","pr.ne.sing.ns.0","pr.ne.nonsing.s.1","pr.ne.sing.s.1","pr.ne.nonsing.s.0","pr.ne.sing.s.0","nseg.sites.1")))

#diff contains difference in singletons in HIV+ and HIV-
#bor1 has Population and group singletons with n=86 vs n=94 and sampling sites instead of strains

bor1  <- matrix(NA,loops,14)
colnames(bor1) <-c("n.singletons","prop.sing","pop.sing.1.86","prop.pop.sing.1.86","group.sing.1.86","prop.group.sing.1.86","nonsing.1.86","prop.nonsing.1.86","pop.sing.0.94","prop.pop.sing.0.94","group.sing.0.94","prop.group.sing.0.94","nonsing.0.94","prop.nonsing.0.94")

for(i in 1:loops){
  #sample the positions in the original data
  genboot  <- sort(sample.int(dim(d)[2], dim(d)[2], replace=T))
  annoboot <- annotation_filtered[,genboot]
  if(i %% 50 == 0) print(paste(i,"--", Sys.time()))
 
  mutboot.94 <- colSums(d[,genboot]==dmut[,genboot])
  glob.sing.94 <- factor(mutboot.94==1, labels = c("non-singleton","singleton"))
  
  #number of singletons in the 180 strains and their proportion
  bor1[i,1] <- table(mutboot.94==1,useNA="always") [2]
  bor1[i,2] <- table(mutboot.94==1,useNA="always") [2]/(table(mutboot.94==1,useNA="always") [2]+table(mutboot.94==1,useNA="always") [1])
  
  #From all the mutations sampled, which ones are in hiv+
  muthivboot <- colSums(d[hiv==1,genboot] == dmut[hiv==1,genboot])
  hiv.factorboot <- cut(muthivboot, breaks=c(0,1,2,10000), labels=c("no-mut","singel","non-sing"), right = F)
  table(hiv.factor[glob.sing.94=="singleton"], useNA="always")

  # The following commands fill the bor matrix,1 receives the proportion of singletons;3 receives the absolute number.
  bor[i,1] <- table(hiv.factorboot)[2]/(table(hiv.factorboot)[2]+table(hiv.factorboot)[3])
  bor[i,3] <- table(muthivboot==1)[2]
  #number of Population singletons
  bor1[i,3] <-table(hiv.factorboot[glob.sing.94=="singleton"], useNA="always")[2]
  # proportion of population singletons
  bor1[i,4] <-table(hiv.factorboot[glob.sing.94=="singleton"], useNA="always")[2]/(table(hiv.factorboot)[2]+table(hiv.factorboot)[3])
  #number of group singletons
  bor1[i,5] <-table(hiv.factorboot[glob.sing.94=="non-singleton"], useNA="always")[2]
  #proportion of group singletons
  bor1[i,6] <- table(hiv.factorboot[glob.sing.94=="non-singleton"], useNA="always")[2]/(table(hiv.factorboot)[2]+table(hiv.factorboot)[3])
  #number of nonsingletons
  bor1[i,7] <- table(hiv.factorboot)[3]
  #proportion of nonsingletons
  bor1[i,8] <-table(hiv.factorboot)[3]/(table(hiv.factorboot)[2]+table(hiv.factorboot)[3])
  
  #Classify bootstrap positions as NS/SYN/I
  t1 <- table(hiv.factorboot, as.character(t(annoboot[4,])), useNA="always")
  
  #nonsingletons I/NS/SYN,"nonsing.i.1","nonsing.ns.1","nonsing.s.1"
  bor[i,5:7] <- t1[3,c(1,2,3)]
  
  #singletons I/NS/SYN,"sing.i.1","sing.ns.1","sing.s.1"
  bor[i,8:10] <- t1[2,c(1,2,3)]
  
  #Classify bootstrap positions as Essential/Non-essential in vitro
  t2 <- table(hiv.factorboot, as.character(t(annoboot[6,])), useNA="always")
  
  #nonsingletons essentials/non-essentials in vitro
  bor[i,17:18] <- t2[3,c(2:3)]
  bor[i,19:20] <- t2[2,c(2:3)]
  bor[i,26:33] <-  as.numeric(table(hiv.factorboot, as.character(t(annoboot[6,])),as.character(t(annoboot[4,])))[c(3:2),2:3,c(2,3)])
  
  #store the total number of variable sites for mtb/hivpos
  bor[i,42] <- (table(hiv.factorboot)[2]+table(hiv.factorboot)[3])
  
  #43:""pr.nonsing.ns.1",44:pr.sing.ns.1"
  bor[i,43:44] <- c(bor[i,6]/bor[i,42],bor[i,9]/bor[i,42]) 
  
  #47:"pr.nonsing.s.1",48:"pr.sing.s.1"
  bor[i,47:48] <- c(bor[i,7]/bor[i,42],bor[i,10]/bor[i,42]) 
  
  #51:"pr.nonsing.i.1",52:"pr.sing.i.1"
  bor[i,51:52] <- c(bor[i,5]/bor[i,42],bor[i,8]/bor[i,42])
  
  #55:"pr.e.nonsing.ns.1"  56:"pr.e.sing.ns.1" 
  bor[i,55:56] <-c(bor[i,26]/bor[i,42],bor[i,27]/bor[i,42])
  
  #59:"pr.e.nonsing.s.1"   60:"pr.e.sing.s.1"
  bor[i,59:60] <-c(bor[i,30]/bor[i,42],bor[i,31]/bor[i,42])
  
  #63:"pr.ne.nonsing.ns.1" 64:"pr.ne.sing.ns.1
  bor[i,63:64] <-c(bor[i,28]/bor[i,42],bor[i,29]/bor[i,42])
  
  #67:pr.ne.nonsing.s.1"68:"pr.ne.sing.s.1"
  bor[i,67:68] <-c(bor[i,32]/bor[i,42],bor[i,33]/bor[i,42])
  
  #The same but for Mtb/HIVneg
  
  mutnegboot <- colSums(d[hiv==0,genboot] == dmut[hiv==0,genboot])
  neg.factorboot <- cut(mutnegboot, breaks=c(0,1,2,10000), labels=c("no-mut","singel","non-sing"), right = F)
  bor[i,2] <- table(neg.factorboot)[2]/(table(neg.factorboot)[2]+table(neg.factorboot)[3])
  bor[i,4] <- (table(mutnegboot==1))[2]
  
  #population and group singletons
  #number of Population singletons
  bor1[i,9] <-table(neg.factorboot[glob.sing.94=="singleton"], useNA="always")[2]
  
  # proportion of population singletons
  bor1[i,10] <-table(neg.factorboot[glob.sing.94=="singleton"], useNA="always")[2]/(table(neg.factorboot)[2]+table(neg.factorboot)[3])
  
  #number of group singletons
  bor1[i,11] <-table(neg.factorboot[glob.sing.94=="non-singleton"], useNA="always")[2]
  
  #proportion of group singletons
  bor1[i,12] <- table(neg.factorboot[glob.sing.94=="non-singleton"], useNA="always")[2]/(table(neg.factorboot)[2]+table(neg.factorboot)[3])
  
  #number of nonsingletons
  bor1[i,13] <- table(neg.factorboot)[3]
  
  #proportion of nonsingletons
  bor1[i,14] <-table(neg.factorboot)[3]/(table(neg.factorboot)[2]+table(neg.factorboot)[3])
  
  t1 <- table(neg.factorboot, as.character(t(annoboot[4,])), useNA="always")
  bor[i,11:13] <- t1[3,c(1,2,3)]
  bor[i,14:16] <- t1[2,c(1,2,3)]
  t2 <- table(neg.factorboot, as.character(t(annoboot[6,])), useNA="always")
  bor[i,21:22] <- t2[3,c(2:3)]
  bor[i,23:24] <- t2[2,c(2:3)]
  bor[i,34:41] <-  as.numeric(table(neg.factorboot, as.character(t(annoboot[6,])),as.character(t(annoboot[4,])))[c(3:2),2:3,c(2,3)])
  
  #number of segregating sites
  bor[i,71] <-(table(neg.factorboot)[2]+table(neg.factorboot)[3])
  
  #43:""pr.nonsing.ns.0",44:pr.sing.ns.0"
  bor[i,45:46] <- c(bor[i,12]/bor[i,71],bor[i,15]/bor[i,71]) 
  
  #47:"pr.nonsing.s.0",48:"pr.sing.s.0"
  bor[i,49:50] <- c(bor[i,13]/bor[i,71],bor[i,16]/bor[i,71]) 
  
  #51:"pr.nonsing.i.0",52:"pr.sing.i.0"
  bor[i,53:54] <- c(bor[i,11]/bor[i,71],bor[i,14]/bor[i,71])
  
  #55:"pr.e.nonsing.ns.0"  56:"pr.e.sing.ns.0" 
  bor[i,57:58] <-c(bor[i,34]/bor[i,71],bor[i,35]/bor[i,71])
  
  #59:"pr.e.nonsing.s.0"   60:"pr.e.sing.s.0"
  bor[i,61:62] <-c(bor[i,38]/bor[i,71],bor[i,39]/bor[i,71])
  
  #63:"pr.ne.nonsing.ns.0" 64:"pr.ne.sing.ns.0
  bor[i,65:66] <-c(bor[i,36]/bor[i,71],bor[i,37]/bor[i,71])
  
  #67:pr.ne.nonsing.s.0"68:"pr.ne.sing.s.0"
  bor[i,69:70] <-c(bor[i,40]/bor[i,71],bor[i,41]/bor[i,71])
  bor[i,25] <- bor[i,3] - bor[i,4]
  }
  
  dbor <- as.data.frame(list(var94=colnames(bor),
          median94=round(apply(bor, 2, median, na.rm=T),digits=2)))
  dbor$lci94 <- round(apply(bor, 2, quantile, probs=0.025, na.rm=T),digits=2)
  dbor$uci94 <- round(apply(bor, 2, quantile, probs=0.975, na.rm=T),digits=2)

  dbor1 <- as.data.frame(list(varGLO94=colnames(bor1),medianGLO94=round(apply(bor1,2,median,na.rm=T),digits=2)))
  dbor1$lciGLO94 <- round(apply(bor1, 2, quantile, probs=0.025, na.rm=T),digits=2)
  dbor1$uciGLO94 <- round(apply(bor1, 2, quantile, probs=0.975, na.rm=T),digits=2)
  
  write.table(dbor,file="Table_bootstrap94", sep="\t", row.names=FALSE)
  write.table(dbor1,file="Table_boot_global_local_singletonsSITES.94",sep="\t",row.names=F)

  
#### As Mtb/HIV- has  n = 94 vs Mtb/HIV+ has n=86 let see if sampling Mtb/HIV- n=86 has impact on the confidence intervals


set.seed(110303)

#create vector with 86 "0" and 86 "1" 

hiv86<- c(rep(0,86),rep(1,86))
    
loops <- 1000

#"glob.sing"=total singletons n=180, "hiv.sing" = total singletons Mtb/HIV+=86, "hiv.loc.sing"=singletons only in HIV+ and not singletons in Mtb/HIV-, "hiv.glob.sing" =singletons in n=180 but exclusively in Mtb/HIV+
#bor2 has Population and group singletons with n=86 for both Mtb/HIVpos and Mtb/HIVneg
bor2 <- matrix(NA, loops, 9)
colnames(bor2) <- c("loci","not.var","glob.sing","hiv.sing","hiv.loc.sing","hiv.glob.sing","neg.sing","neg.loc.sing","neg.glob.sing")

bor3 <-  matrix(NA,loops,33+38)
colnames(bor3) <-c(c("prop1.86","prop0.86","bn1.86","bn0.86","nonsing.i.1.86","nonsing.ns.1.86","nonsing.s.1.86","sing.i.1.86","sing.ns.1.86","sing.s.1.86","ns.i.0.86","ns.ns.0.86","ns.s.0.86","s.i.0.86","s.ns.0.86","s.s.0.86","nonsing.e.1.86","nonsing.ne.1.86","sing.e.1.86","sing.ne.1.86","ns.e.0.86","ns.ne.0.86","s.e.0.86","s.ne.0.86","dif.86")
,
paste(c(paste(c("nonsing","sing"),".ns.",sort(rep(c("e","ne"),2)),sep="") ,paste(c("nonsing","sing"),".s.",sort(rep(c("e","ne"),2)),sep="")),".1.86",sep="")
,
paste(c(paste(c("nonsing","sing"),".ns.",sort(rep(c("e","ne"),2)),sep="") ,paste(c("nonsing","sing"),".s.",sort(rep(c("e","ne"),2)),sep="")),".0.86",sep="")
,
paste(c("nseg.sites.1.86","pr.nonsing.ns.1.86","pr.sing.ns.1.86","pr.nonsing.ns.0.86","pr.sing.ns.0.86","pr.nonsing.s.1.86","pr.sing.s.1.86","pr.nonsing.s.0.86","pr.sing.s.0.86","pr.nonsing.i.1.86","pr.sing.i.1.86","pr.nonsing.i.0.86","pr.sing.i.0.86","pr.e.nonsing.ns.1.86","pr.e.sing.ns.1.86","pr.e.nonsing.ns.0.86","pr.e.sing.ns.0.86","pr.e.nonsing.s.1.86","pr.e.sing.s.1.86","pr.e.nonsing.s.0.86","pr.e.sing.s.0.86","pr.ne.nonsing.ns.1.86","pr.ne.sing.ns.1.86","pr.ne.nonsing.ns.0.86","pr.ne.sing.ns.0.86","pr.ne.nonsing.s.1.86","pr.ne.sing.s.1.86","pr.ne.nonsing.s.0.86","pr.ne.sing.s.0.86","nseg.sites.1.86")))

bor4 <- matrix(NA,loops,12)
colnames(bor4) <-c("pop.sing.1.86","prop.pop.sing.1.86","group.sing.1.86","prop.group.sing.1.86","nonsing.1.86","prop.nonsing.1.86","pop.sing.0.86","prop.pop.sing.0.86","group.sing.0.86","prop.group.sing.0.86","nonsing.0.86","prop.nonsing.0.86")

  
#Create a random sample of 86 numbers from 94 and extract those from the MTB/HIVneg isolates. This is now d86 and it has now 86 HIV+ and 86 HIV-.	
#d needs to be sorted according to HIV status

  for(i in 1:loops){
  if(i %% 50 == 0) print(paste(i,"--", Sys.time()))
  d86 <- d[c(sort(sample.int(94,86,replace=F)),95:180),]
  dim(d86) 
  not.variable <- apply(d86,2,is.not.variable)
  bor2[i,1:2] <- table(not.variable)
  to.removeb <- which(not.variable)
  d86 <- d86[, -c(to.removeb)]
  mut86 <-  mut[,-c(to.removeb)]
  annotation_filtered86 <-  annotation_filtered[,-c(to.removeb)]

  #Create a table with the information on the mutations, 172 rows (86*2) and 4581 collumns 
  dmut86 <- mut86[rep(1,86*2),]
  dhiv86 <- d86[hiv86==1,]
  dneg86 <- d86[hiv86==0,]
  
  #Count the mutation in d86, assign the number of singletons to the 3rd collumn of bor2
  ismut86 <- (d86 == dmut86)
  mutg86 <- colSums(ismut86, na.rm=T)
  bor2[i,3] <- table(mutg86==1)[2]
 
  #Count the mutation in MTB/HIV+
  ismuthiv86 <- (d86[hiv86==1,] == dmut86[hiv86==1,])
  mut186 <- colSums(ismuthiv86, na.rm=T)

  #in column 4 store the the number of singletons in Mtb/HIV+ called "hiv.sing"
  bor2[i,4] <- table(mut186==1)[2]
  
  ismutneg86 <- (d86[hiv86==0,] == dmut86[hiv86==0,])
  mut086 <- colSums(ismutneg86, na.rm=T)
  
  #in column 7 store the the number of singletons in Mtb/HIV- called "neg.sing"
  bor2[i,7] <- table(mut086==1)[2]

  #in collumn 5 is the number of positions which are singletons in Mtb/HIV+ but non-singleton in Mtb/HIV- "hiv.loc.sing";in column 6 are the n° of positions wich are singletons in HIV+ and have 0 counts in HIV-, .
  #in collumn 8 is the number of positions which are singletons in Mtb/HIV- but non-singleton in Mtb/HIV+ "neg.loc.sing";in column 9 are the n° of positions wich are singletons in HIV- and have 0 counts in HIV+,"neg.glob.sing" 
  bor2[i,5:6] <- table(mutg86==1,mut186==1)[,2]
  bor2[i,8:9] <- table(mutg86==1,mut086==1)[,2]

  #sample ns/s/i in essential and non-essential in vitro.  
  genboot86  <- sort(sample.int(dim(d86)[2], dim(d86)[2], replace=T))
  annoboot86 <- annotation_filtered86[,genboot86]
  
  mutboot.86 <- colSums(d86[,genboot86]==dmut86[,genboot86])
  glob.sing.86 <- factor(mutboot.86==1, labels = c("non-singleton","singleton"))
  
  
  #From all the mutations sampled, which ones are in Mtb/hiv+. This has been calculated already before as the number of Mtb/HIV+ does not change. 
  muthivboot86 <- colSums(dhiv86[,genboot86] == dmut86[hiv86==1,genboot86])
  hiv.factorboot86 <- cut(muthivboot86, breaks=c(0,1,2,10000), labels=c("no-mut","singel","non-sing"), right = F)
 
  #number of Population singletons
  bor4[i,1] <-table(hiv.factorboot86[glob.sing.86=="singleton"], useNA="always")[2]
  
  # proportion of population singletons
  bor4[i,2] <-table(hiv.factorboot86[glob.sing.86=="singleton"], useNA="always")[2]/(table(hiv.factorboot86)[2]+table(hiv.factorboot86)[3])
  
  #number of group singletons
  bor4[i,3] <-table(hiv.factorboot86[glob.sing.86=="non-singleton"], useNA="always")[2]
  
  #proportion of group singletons
  bor4[i,4] <- table(hiv.factorboot86[glob.sing.86=="non-singleton"], useNA="always")[2]/(table(hiv.factorboot86)[2]+table(hiv.factorboot86)[3])
  
  #number of nonsingletons
  bor4[i,5] <- table(hiv.factorboot86)[3]
  
  #proportion of nonsingletons
  bor4[i,6] <-table(hiv.factorboot86)[3]/(table(hiv.factorboot86)[2]+table(hiv.factorboot86)[3])
  
  #daniela: proportion of singletons. 
  bor3[i,1] <- table(hiv.factorboot86)[2]/(table(hiv.factorboot86)[2]+table(hiv.factorboot86)[3])
  
  #daniela:number of singletons 
  bor3[i,3] <- table(muthivboot86==1)[2]
  
  t1 <- table(hiv.factorboot86, as.character(t(annoboot86[4,])), useNA="always")
  bor3[i,5:7] <- t1[3,c(1,2,3)]
  bor3[i,8:10] <- t1[2,c(1,2,3)]
  
  t2 <- table(hiv.factorboot86, as.character(t(annoboot86[6,])), useNA="always")
  bor3[i,17:18] <- t2[3,c(2:3)]
  bor3[i,19:20] <- t2[2,c(2:3)]
  
  #daniela:nonsingletons essentials/non-essentials I/NS/SYN
  bor3[i,26:33] <-  as.numeric(table(hiv.factorboot86, as.character(t(annoboot86[6,])),as.character(t(annoboot86[4,])))[c(3:2),2:3,c(2,3)])
  bor3[i,42] <- (table(hiv.factorboot86)[2]+table(hiv.factorboot86)[3])
  
  #43:""pr.nonsing.ns.1.86",44:pr.sing.ns.1.86"
  bor3[i,43:44] <- c(bor3[i,6]/bor3[i,42],bor3[i,9]/bor3[i,42]) 
  
  #47:"pr.nonsing.s.1.86",48:"pr.sing.s.1.86"
  bor3[i,47:48] <- c(bor3[i,7]/bor3[i,42],bor3[i,10]/bor3[i,42]) 
  
  #51:"pr.nonsing.i.1.86",52:"pr.sing.i.1.86"
  bor3[i,51:52] <- c(bor3[i,5]/bor3[i,42],bor3[i,8]/bor3[i,42])
  
  #55:"pr.e.nonsing.ns.1.86"  56:"pr.e.sing.ns.1.86" 
  bor3[i,55:56] <-c(bor3[i,26]/bor3[i,42],bor3[i,27]/bor3[i,42])
  
  #59:"pr.e.nonsing.s.1.86"   60:"pr.e.sing.s.1.86"
  bor3[i,59:60] <-c(bor3[i,30]/bor3[i,42],bor3[i,31]/bor3[i,42])
  
  #63:"pr.ne.nonsing.ns.1.86" 64:"pr.ne.sing.ns.1.86
  bor3[i,63:64] <-c(bor3[i,28]/bor3[i,42],bor3[i,29]/bor3[i,42])
  
  #67:pr.ne.nonsing.s.1.86"68:"pr.ne.sing.s.1.86"
  bor3[i,67:68] <-c(bor3[i,32]/bor3[i,42],bor3[i,33]/bor3[i,42])
  
  ####Mtb/HIVneg
  
  mutnegboot86 <- colSums(dneg86[,genboot86] == dmut86[hiv86==0,genboot86])
  neg.factorboot86 <- cut(mutnegboot86, breaks=c(0,1,2,10000), labels=c("no-mut","singel","non-sing"), right = F)
  
  #number of Population singletons
  bor4[i,7] <-table(neg.factorboot86[glob.sing.86=="singleton"], useNA="always")[2]
  
  # proportion of population singletons
  bor4[i,8] <-table(neg.factorboot86[glob.sing.86=="singleton"], useNA="always")[2]/(table(neg.factorboot86)[2]+table(neg.factorboot86)[3])
  
  #number of group singletons
  bor4[i,9] <-table(neg.factorboot86[glob.sing.86=="non-singleton"], useNA="always")[2]
  
  #proportion of group singletons
  bor4[i,10] <- table(neg.factorboot86[glob.sing.86=="non-singleton"], useNA="always")[2]/(table(neg.factorboot86)[2]+table(neg.factorboot86)[3])
  
  #number of nonsingletons
  bor4[i,11] <- table(neg.factorboot86)[3]
  #proportion of nonsingletons
  bor4[i,12] <-table(neg.factorboot86)[3]/(table(neg.factorboot86)[2]+table(neg.factorboot86)[3])
  bor3[i,2] <- table(neg.factorboot86)[2]/(table(neg.factorboot86)[2]+table(neg.factorboot86)[3])
  bor3[i,4] <- table(mutnegboot86==1) [2]
  
  t1 <- table(neg.factorboot86, as.character(t(annoboot86[4,])), useNA="always")
  bor3[i,11:13] <- t1[3,c(1,2,3)]
  bor3[i,14:16] <- t1[2,c(1,2,3)]
  t2 <- table(neg.factorboot86, as.character(t(annoboot86[6,])), useNA="always")
  bor3[i,21:22] <- t2[3,c(2:3)]
  bor3[i,23:24] <- t2[2,c(2:3)]
  bor3[i,34:41] <- as.numeric(table(neg.factorboot86, as.character(t(annoboot86[6,])),as.character(t(annoboot86[4,])))[c(3:2),2:3,c(2,3)])
  
  #number if segregating sites
  bor3[i,71] <-(table(neg.factorboot86)[2]+table(neg.factorboot86)[3])
  #43:""pr.nonsing.ns.0.86",44:pr.sing.ns.0.86"
  bor3[i,45:46] <- c(bor3[i,12]/bor3[i,71],bor3[i,15]/bor3[i,71]) 
  #47:"pr.nonsing.s.0.86",48:"pr.sing.s.0.86"
  bor3[i,49:50] <- c(bor3[i,13]/bor3[i,71],bor3[i,16]/bor3[i,71]) 
  #51:"pr.nonsing.i.0.86",52:"pr.sing.i.0.86"
  bor3[i,53:54] <- c(bor3[i,11]/bor3[i,71],bor3[i,14]/bor3[i,71])
  #55:"pr.e.nonsing.ns.0.86"  56:"pr.e.sing.ns.0.86" 
  bor3[i,57:58] <-c(bor3[i,34]/bor3[i,71],bor3[i,35]/bor3[i,71])
  #59:"pr.e.nonsing.s.0.86"   60:"pr.e.sing.s.0.86"
  bor3[i,61:62] <-c(bor3[i,38]/bor3[i,71],bor3[i,39]/bor3[i,71])
  #63:"pr.ne.nonsing.ns.0.86" 64:"pr.ne.sing.ns.0.86
  bor3[i,65:66] <-c(bor3[i,36]/bor3[i,71],bor3[i,37]/bor3[i,71])
  #67:pr.ne.nonsing.s.0.86"68:"pr.ne.sing.s.0.86"
  bor3[i,69:70] <-c(bor3[i,40]/bor3[i,71],bor3[i,41]/bor3[i,71])
  bor3[i,25] <- bor3[i,3]-bor3[i,4]
                }
 
dbor2 <- as.data.frame(list(varGLO=colnames(bor2),
          medianGLO=round(apply(bor2, 2, median, na.rm=T),digits=2)))
dbor2$lciGLO <- round(apply(bor2, 2, quantile, probs=0.025, na.rm=T),digits=2)
dbor2$uciGLO <- round(apply(bor2, 2, quantile, probs=0.975, na.rm=T),digits=2)

write.table(dbor2, file="Table_boot_global_local_singletons", sep="\t", row.names=FALSE)
 
 
dbor3 <- as.data.frame(list(var86=colnames(bor3),
          median86=round(apply(bor3, 2, median, na.rm=T),digits=2)))
dbor3$lci86 <- round(apply(bor3, 2, quantile, probs=0.025, na.rm=T),digits=2)
dbor3$uci86 <- round(apply(bor3, 2, quantile, probs=0.975, na.rm=T),digits=2)

write.table(dbor3, file="Table_bootstrap86", row.names=FALSE)
 
dbor4 <- as.data.frame(list(varGLO86=colnames(bor4),medianGLO86=round(apply(bor4,2,median,na.rm=T),digits=2)))
dbor4$lciGLO86 <-  round(apply(bor4, 2, quantile, probs=0.025, na.rm=T),digits=2)
dbor4$uciGLO86 <- round(apply(bor4, 2, quantile, probs=0.975, na.rm=T),digits=2)
 
write.table(dbor4, file="Table_boot_global_local_singletonsSITES.86", row.names=FALSE)
 
 
######permutation tests for Singletons, Population singletons, Group singletons, Nonsingletons
 
 set.seed(777)
 loops <-1000
 per <-  matrix(NA,loops,12)  
 colnames(per) <-c("pop1sing","pop2sing","diffSing","PopSingpop1","PopSingpop2","diffPopsing","GroupSingpop1","GroupSingpop2","diffGroupsing","nonSingpop1","nonSingpop2","diffnonSin")
 
 
 for(i in 1:loops){
  if(i%%10==0)cat(paste(i/10, ".", sep=""))
  perm <- sample.int(dim(d)[1],dim(d)[1])
  pop1 <- sort(perm[1:86])
  pop2 <- sort(perm[-(1:86)])
 
  #number of singletons
  ismutpop1 <- colSums(d[pop1,]==dmut[pop1,])
  per[i,1]<- table(ismutpop1==1)[2] 
 
  ismutpop2 <- colSums(d[pop2,]==dmut[pop2,])
  per[i,2] <- table(ismutpop2==1)[2] 
 
  per[i,3] <- table(ismutpop1==1)[2]-table(ismutpop2==1)[2]

  #number of population singletons
  per [i,4] <- table(glob.singl=="singleton",ismutpop1==1)[2,2] 
  per [i,5] <- table(glob.singl=="singleton",ismutpop2==1)[2,2] 
  per [i,6] <- table(glob.singl=="singleton",ismutpop1==1)[2,2]-table(glob.singl=="singleton",ismutpop2==1)[2,2]
 
  #number of group singletons
  per [i,7] <- table(glob.singl=="singleton",ismutpop1==1)[1,2] 
  per [i,8] <- table(glob.singl=="singleton",ismutpop2==1)[1,2] 
  per [i,9] <- table(glob.singl=="singleton",ismutpop1==1)[1,2]-table(glob.singl=="singleton",ismutpop2==1)[1,2]
 
  #nonsingletons
  per[i,10]<- table(ismutpop1>1)[2] 
  per[i,11] <- table(ismutpop2>1)[2]
  per[i,12] <- table(ismutpop1>1)[2] -table(ismutpop2>1)[2]
   }
 
##Observed differences and p values for the permutation
##p values- multiply by 2 go get asymptotic 2 sided p value. 
 
#singletons Mtb/HIV+ -Mtb/HIV-
2121-1697 

colMeans(per[,1:2])
quantile(per[,3])
#p value for permutation 
table(per[,3]>424)/loops*2

#Population singletons
1279-1141
 
colMeans(per[,4:5])
quantile(per[,6])
#p value for permutation
table(per[,6]>138)/loops*2

#Group singletons
842-546

colMeans(per[,7:8])
quantile(per[,9])
table(per[,9]>296)/loops*2

#Nonsingletons
1181-1433

colMeans(per[,10:11])
table(per[,12]<=-252)[2]/1000
quantile(per[,12]) 
table(per[,12]<(-252))/loops*2
