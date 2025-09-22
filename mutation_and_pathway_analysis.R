source('/Users/cfungtammasan/Analysis/setup_slim_R.R')
source('/Users/cfungtammasan/gitlab/tmp_code/scripts/functions_ML.R')
library(maftools)
library(dplyr)
library(data.table)
library(readxl)
library(ggsignif)


#####################
## Helper function ##
#####################

variant_type_cal <- function(ref,alt){
  length_diff=nchar(alt)-nchar(ref)
  case_when(
    length_diff==0 ~ 'SNP',
    length_diff>0 ~ 'INS',
    length_diff<0 ~ 'DEL'
  )
}

map_SnpEff2MAF_classification <- function(effect,Variant_Type,length_ref,length_alt) {
  case_when(
    is.na(effect) | trimws(effect)=="" ~ "Targeted_Region",
    grepl("splice_acceptor_variant|splice_donor_variant", effect) ~ "Splice_Site",
    grepl("stop_gained", effect) ~ "Nonsense_Mutation",
    grepl("frameshift_variant", effect) & Variant_Type=="DEL" ~ "Frame_Shift_Del",
    grepl("frameshift_variant", effect) & Variant_Type=="INS" ~ "Frame_Shift_Ins",
    grepl("protein_protein_contact|structural_interaction_variant", effect) & Variant_Type=="DEL" & (length_ref-length_alt)%%3!=0 ~ "Frame_Shift_Del",
    grepl("protein_protein_contact|structural_interaction_variant", effect) & Variant_Type=="INS" & (length_alt-length_ref)%%3!=0 ~ "Frame_Shift_Ins",
    grepl("stop_lost", effect) ~ "Nonstop_Mutation",
    grepl("start_lost|initiator_codon_variant", effect) ~ "Translation_Start_Site",
    grepl("inframe_insertion", effect) ~ "In_Frame_Ins",
    grepl("inframe_deletion", effect) ~ "In_Frame_Del",
    grepl("protein_protein_contact|structural_interaction_variant", effect) & Variant_Type=="INS" & (length_alt-length_ref)%%3==0~ "In_Frame_Ins",
    grepl("protein_protein_contact|structural_interaction_variant", effect) & Variant_Type=="DEL" & (length_ref-length_alt)%%3==0 ~ "In_Frame_Del",
    grepl("protein_protein_contact", effect) ~ "Missense_Mutation", 
    grepl("structural_interaction_variant", effect) ~ "Missense_Mutation",  
    grepl("missense_variant|rare_amino_acid_variant", effect) ~ "Missense_Mutation",
    grepl("intron_variant", effect) ~ "Intron",
    grepl("splice_region_variant", effect) ~ "Splice_Region",
    grepl("synonymous_variant|stop_retained_variant", effect) ~ "Silent",
    grepl("mature_miRNA_variant|exon_variant|non_coding_exon_variant|non_coding_transcript_exon_variant|non_coding_transcript_variant|nc_transcript_variant", effect) ~ "RNA",
    grepl("5_prime_UTR_variant|5_prime_UTR_premature_start_codon_gain_variant", effect) ~ "5'UTR",
    grepl("3_prime_UTR_variant", effect) ~ "3'UTR",
    grepl("TF_binding_site_variant|regulatory_region_variant|regulatory_region|intergenic_variant|intergenic_region", effect) ~ "IGR",
    grepl("upstream_gene_variant", effect) ~ "5'Flank",
    grepl("downstream_gene_variant", effect) ~ "3'Flank",
    TRUE ~ "Unknown"
  )
}

map_impact <- function(impact) {
  case_when(
    is.na(impact) | trimws(impact)=="" ~ NA,
    grepl("HIGH", impact) ~ "HIGH",
    grepl("MODERATE", impact) ~ "MODERATE",
    grepl("LOW", impact) ~ "LOW",
    grepl("MODIFIER", impact) ~ "MODIFIER"
  )
}

unique_element<-function(input_string,max_unq_count=1){
  unq_count=ifelse(
    grepl(',',input_string),
    length(unique(trimws(unlist(strsplit(input_string,','))))),
    1
  )
  case_when(
    unq_count == 1 ~ input_string,
    unq_count <= max_unq_count ~ paste(unique(trimws(unlist(strsplit(input_string,',')))),collapse=','),
    unq_count > max_unq_count ~ NA
  )
}

##############################################
## Read annotated data and transform to maf ##
##############################################

all_data=fread('all_data_WESdata20250725.txt')
stopifnot(length(unique(all_data[,`Registration number`]))==5549)
length(unique(all_data[WGD=='positive',`Registration number`])) #2716
length(unique(all_data[WGD=='negative',`Registration number`])) #2833
all_data[,WGD:=toupper(WGD)]

subset_data=all_data[,.(Hugo_Symbol=unique_element(snpEff_ANN_Gene_Name),
                        Entrez_Gene_Id=unique_element(snpEff_ID),
                        Center='',
                        NCBI_Build='GRCh37',
                        Chromosome=snpEff_CHROM,
                        Start_Position=snpEff_POS,
                        End_Position='',
                        Strand='',
                        Variant_Classification="",
                        Variant_Type=variant_type_cal(snpEff_REF...229,snpEff_ALT),
                        Reference_Allele=snpEff_REF...229,
                        Tumor_Seq_Allele1=snpEff_REF...229,
                        Tumor_Seq_Allele2=snpEff_ALT,
                        dbSNP_RS='',
                        dbSNP_Val_Status='',
                        Tumor_Sample_Barcode=`Registration number`,
                        IMPACT='HIGH',#map_impact(snpEff_ANN_Annotation_Impact),
                        copynumber,
                        WGD,
                        ploidy,
                        `Registration number`,
                        Age,
                        Sex,
                        snpEff_ANN_Annotation,
                        snpEff_TYPE,
                        snpEff_QUAL,
                        snpEff_CIGAR,
                        snpEff_LEN,
                        snpEff_ANN_Allele,
                        snpEff_ANN_Transcript_BioType,
                        snpEff_ANN_HGVS_c,
                        snpEff_ANN_HGVS_p)]

length(unique(subset_data$`Registration number`))
subset_data[,Variant_Classification:=map_SnpEff2MAF_classification(snpEff_ANN_Annotation,Variant_Type,nchar(Reference_Allele),nchar(Tumor_Seq_Allele2))]
subset_data[,Var.ID.DNA:=paste(Chromosome,Start_Position, Reference_Allele,Tumor_Seq_Allele2, sep="-")] # found that we don't need snpEff_TYPE because same variant could be classified as mnp or complex
subset_data[, Var.ID.Prot:=paste(Hugo_Symbol, snpEff_ANN_HGVS_p, sep="_")]
subset_data[snpEff_ANN_HGVS_p == "", Var.ID.Prot:= paste(Hugo_Symbol, Var.ID.DNA, sep="_")]
subset_data[,Freq.DNA:=.N, by=Var.ID.DNA ]
subset_data[,Freq.DNA.All:=Freq.DNA]
subset_data[!Variant_Type =='SNP', PolyX_ID:=paste(Chromosome,as.numeric(Start_Position)+1,sep='_')]
subset_data %<>% mergeR(dt.polyX[,.(PolyX_ID, Rep, Len)], by='PolyX_ID')
length(unique(subset_data$`Registration number`))

subset_data[,Freq.Prot:=.N, by=Var.ID.Prot]
subset_data %<>% mergeR(unique(subset_data[,.(Hugo_Symbol, Tumor_Sample_Barcode)])[,.(Freq.Gene=.N), by="Hugo_Symbol"], by.x="Hugo_Symbol", by.y="Hugo_Symbol") 
length(unique(subset_data$`Registration number`))

subset_data[Variant_Classification %in% c('Missense_Mutation','Frame_Shift_Ins','In_Frame_Del',
'Frame_Shift_Del','Nonsense_Mutation','In_Frame_Ins','Nonstop_Mutation','Splice_Site',
'Translation_Start_Site') , Variant_Classification2 := "Non-synonymous mutation"]
length(unique(subset_data$`Registration number`))

subset_data[Variant_Classification == "Silent" , Variant_Classification2 := "Synonymous mutation"]
subset_data[is.na(subset_data$Variant_Classification2), Variant_Classification2 := "Other"]
subset_data[Hugo_Symbol %in% Cancer.genes$Hugo.Symbol, Is.Cancer.gene:="Yes"]
length(unique(subset_data$`Registration number`))

Total.tumors <- as.numeric(unqN(subset_data$Tumor_Sample_Barcode))
subset_data[,Prev.DNA:=100*Freq.DNA/Total.tumors]
subset_data[,Prev.Gene:=100*Freq.Gene/Total.tumors]
subset_data=subset_data[!is.na(Hugo_Symbol)]
length(unique(subset_data$`Registration number`)) # lost 2 individual

subset_data[,FLAGS:=Hugo_Symbol %in% top100flags[1:20]]
dim(subset_data[WGD=="POSITIVE"])
dim(subset_data[WGD=="NEGATIVE"])
length(unique(subset_data$`Registration number`))

write.csv(subset_data,'subset_data_20250725.csv',row.names=F)

########################
## Merge MSI TMB info ##
########################

subset_data=fread('subset_data_20250725.csv')
length(unique(subset_data$Tumor_Sample_Barcode)) # 5547

TMB_data=fread("George_WGD20240402.tsv")
table(TMB_data$TMB,useNA='always')
TMB_data=TMB_data[,.(Tumor_Sample_Barcode=PtsID,TMB)]
MSI_data=fread('MSIstatus20250210.xlsx - Sheet1.csv')
MSI_data=MSI_data[,.(Tumor_Sample_Barcode=PtsID,MSI)]
table(MSI_data$MSI,useNA='always')

subset_data=MSI_data[subset_data,on='Tumor_Sample_Barcode']
subset_data=TMB_data[subset_data,on='Tumor_Sample_Barcode']
table(subset_data[,.(WGD,MSI)],useNA='always')
table(subset_data[,.(WGD,TMB)],useNA='always')

write.csv(subset_data,'subset_data_20250725_withMSITMB.csv',row.names=F)

#####################
## Load final data ##
#####################

subset_data=fread('subset_data_20250725_withMSITMB.csv')

##############
## Oncoplot ##
##############

positive_subset_data=subset_data[WGD=="POSITIVE" & Prev.Gene > 5 & Freq.Gene > 1]
negative_subset_data=subset_data[WGD=="NEGATIVE" & Prev.Gene > 5 & Freq.Gene > 1]


# pdf('oncoplot_pathway_WGD.pdf',width=40,height=8)
png('oncoplot_pathway_WGD.png',width=1800,height=600)
# oncoplot(maf = laml_subset, pathways = "sigpw", gene_mar = 8, fontSize = 0.4, topPathways = 10)
laml_subset <- read.maf(maf = positive_subset_data)
oncoplot(maf = laml_subset, pathways = "sigpw", gene_mar = 8, fontSize = 0.4, topPathways = 5, path_order=c('TP53','WNT','RTK-RAS','NOTCH','Hippo'),titleText='WGD+')
laml_subset <- read.maf(maf = negative_subset_data)
oncoplot(maf = laml_subset, pathways = "sigpw", gene_mar = 8, fontSize = 0.4, topPathways = 5,path_order=c('TP53','WNT','RTK-RAS','NOTCH','Hippo'),titleText='WGD-')
# oncoplot(maf = laml_subset, pathways = "sigpw", gene_mar = 8, fontSize = 0.4, selectedPathways='TP53')
dev.off()

png('oncoplot_pathway_WGD_positive.png',width=3200,height=2000)
# pdf('oncoplot_pathway_WGD_positive.pdf',width=40,height=20)
laml_subset <- read.maf(maf = positive_subset_data)
oncoplot(maf = laml_subset, pathways = "sigpw", gene_mar = 8, fontSize = 1, legendFontSize=3, topPathways = 5, path_order=c('TP53','WNT','RTK-RAS','NOTCH','Hippo'),titleText='WGD+',titleFontSize = 2)
dev.off()
png('oncoplot_pathway_WGD_negative.png',width=3200,height=2000)
# pdf('oncoplot_pathway_WGD_negative.pdf',width=40,height=20)
laml_subset <- read.maf(maf = negative_subset_data)
oncoplot(maf = laml_subset, pathways = "sigpw", gene_mar = 8, fontSize = 1, legendFontSize=3, topPathways = 5, path_order=c('TP53','WNT','RTK-RAS','NOTCH','Hippo'),titleText='WGD-',titleFontSize = 2)
dev.off()


#########################
## mutation enrichment ##
#########################
subset_data2=subset_data[is.na(Len) | Len<=7]
subset_data3=subset_data2[IMPACT %in% c("HIGH","MODERATE")]
subset_data3[MSI=='MSI-High',MSI.status:='MSI-High']
subset_data3[MSI=='MSS',MSI.status:='MSS']
subset_data3[,TMB.status:=TMB]


subset_data3=subset_data3[Prev.Gene > 2 & Freq.Gene > 5] 
# Prev.Gene = % of people that carry it; Freq.Gene number of people that carry it
list_A <- unique(subset_data3[WGD=='POSITIVE', Tumor_Sample_Barcode])
list_B <- unique(subset_data3[WGD=='NEGATIVE', Tumor_Sample_Barcode])
name_A <- "WGD+"
name_B <- "WGD-"
compare_genomic_groups_genes(subset_data3, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.02,min_gene_Freq=1, Plot_n_genes=30 )


subset_data3=subset_data3[Prev.Gene > 2 & Freq.Gene > 5] 
list_A <- unique(subset_data3[MSI.status=='MSI-High' & WGD=='POSITIVE', Tumor_Sample_Barcode])
list_B <- unique(subset_data3[MSI.status=='MSI-High' & WGD=='NEGATIVE', Tumor_Sample_Barcode])
name_A <- "MSI_WGD+"
name_B <- "MSI_WGD-"
compare_genomic_groups_genes(subset_data3, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.02,min_gene_Freq=1, Plot_n_genes=30 )


subset_data3=subset_data3[Prev.Gene > 2 & Freq.Gene > 5] 
list_A <- unique(subset_data3[TMB.status=='TMB-High' & MSI.status=='MSI-High' & WGD=='POSITIVE', Tumor_Sample_Barcode])
list_B <- unique(subset_data3[TMB.status=='TMB-High' & MSI.status=='MSI-High' & WGD=='NEGATIVE', Tumor_Sample_Barcode])
name_A <- "TMB_High_MSI_WGD+"
name_B <- "TMB_High_MSI_WGD-"
compare_genomic_groups_genes(subset_data3, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.02,min_gene_Freq=1, Plot_n_genes=30 )


subset_data3=subset_data3[Prev.Gene > 2 & Freq.Gene > 5] 
list_A <- unique(subset_data3[MSI.status=='MSS' & WGD=='POSITIVE', Tumor_Sample_Barcode])
list_B <- unique(subset_data3[MSI.status=='MSS' & WGD=='NEGATIVE', Tumor_Sample_Barcode])
name_A <- "MSS_WGD+"
name_B <- "MSS_WGD-"
compare_genomic_groups_genes(subset_data3, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.02,min_gene_Freq=1, Plot_n_genes=30 )


subset_data3=subset_data3[Prev.Gene > 2 & Freq.Gene > 5] 
list_A <- unique(subset_data3[TMB.status=='TMB-Low' & MSI.status=='MSS' & WGD=='POSITIVE', Tumor_Sample_Barcode])
list_B <- unique(subset_data3[TMB.status=='TMB-Low' & MSI.status=='MSS' & WGD=='NEGATIVE', Tumor_Sample_Barcode])
name_A <- "TMB_Low_MSS_WGD+"
name_B <- "TMB_Low_MSS_WGD-"
compare_genomic_groups_genes(subset_data3, min_Freq=5, list_A=list_A, list_B=list_B, name_A=name_A, name_B=name_B,fraction = 0.02,min_gene_Freq=1, Plot_n_genes=30 )

