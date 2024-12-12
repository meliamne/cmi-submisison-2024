library(tidyverse)
library(ggplot2)
library(glue)
library(ComplexHeatmap)
library(pcaMethods)
library(logger)
library(circlize)


filelocs=list.files("/user/leuven/342/vsc34263/myvibdata/challenge/data_pp",full.names=T,recursive=T,pattern="pca_scores.*.tsv")

# filelocs=list(
# 	cellfreq="/user/leuven/342/vsc34263/myvibdata/challenge/data_pp/cellfreq/pca_scores.tsv",
# 	anti="/user/leuven/342/vsc34263/myvibdata/challenge/data_pp/anti/pca_scores.tsv",
# 	# metadata="data/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/training_subject_specimen.tsv",
# 	cyto="/user/leuven/342/vsc34263/myvibdata/challenge/data_pp/cyto/pca_scores.tsv",
# 	tact="/user/leuven/342/vsc34263/myvibdata/challenge/data/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/training_t_cell_activation_wide.tsv",
# 	tpol="/user/leuven/342/vsc34263/myvibdata/challenge/data_pp/tpol/pca_scores.tsv"
# )

dataset="merged"
data_train=read_tsv("data/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/training_subject_specimen.tsv")
data_challenge=read_tsv("data/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/challenge_subject_specimen.tsv")
data=rbind(data_train,data_challenge)
# data=data %>% filter(timepoint>=0) #filter out early time points due to high missingness in training data # actually dont do this because the challenge data has timepoints <0
data$age=(date("2024-01-01") - data$year_of_birth)/365
#get rid of rarely seen or missing race ingo
data$race_clean=case_when( 
	data$race=="White" ~ "White",
	data$race=="Asian" ~ "Asian",
	T ~ NA
)


all_subjects = unique(data$subject_id)
cormatsize=1200
ndim=15
keeppca=10

colnames(data)

dir.create(glue("./plots/{dataset}"), showWarnings = FALSE, recursive = TRUE)
dir.create(glue("./data_pp/{dataset}"), showWarnings = FALSE, recursive = TRUE)


for(i in 1:length(filelocs)){
	file=filelocs[[i]]
	imputed=ifelse(str_detect(file,"imputed"),"_imputed","")
	corrected=ifelse(str_detect(file,"corrected"),"_corrected","")
	data_modality=str_extract(file,"cellfreq|anti|cyto|tact|tpol")
	data_add=read_tsv(file)
	colnames(data_add)[-1]=glue("{colnames(data_add)[-1]}_{data_modality}{imputed}{corrected}")
	colnames(data_add)[1]="specimen_id"
	newcols=colnames(data_add)[-1]
	overlapcols=intersect(newcols,colnames(data))
	if(length(overlapcols)>0){
		log_warn(glue_collapse("Columns {overlapcols} already present in {dataset}"),sep="\n")
	}
	data=full_join(data,data_add,by="specimen_id")
}


# for(i in 1:length(filelocs)){
# 	file=filelocs[[i]]
# 	data_add=read_tsv(file)
# 	colnames(data_add)[-1]=paste(colnames(data_add)[-1],names(filelocs)[i],sep="_")
# 	newcols=colnames(data_add)[-1]
# 	overlapcols=intersect(newcols,colnames(data))
# 	if(length(overlapcols)>0){
# 		log_warn(glue_collapse("Columns {overlapcols} already present in {dataset}"),sep="\n")
# 	}
# 	data=full_join(data,data_add,by="specimen_id")
# }
# head(data)
# colnames(data)


data_sub=select(data,matches("PC\\d+_"))
numberofnotna=data_sub %>% apply(2,function(x) sum(!is.na(x))) %>% as.character()
annot = HeatmapAnnotation(obsnum=anno_text(numberofnotna))
mat=data_sub %>%  as.matrix()
cormat=cor(mat,use="pairwise.complete.obs")
png(file=glue("plots/{dataset}/correlation_heatmap.png"),width=cormatsize,height=cormatsize)
ht=Heatmap(cormat, name="correlation_raw", show_row_names=T, show_column_names=FALSE, cluster_rows=T, cluster_columns=T, top_annotation=annot, col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))
print(ht)
dev.off()
png(file=glue("plots/{dataset}/correlation_heatmap_lab.png"),width=cormatsize,height=cormatsize)
ht=Heatmap(cormat, name="correlation_raw", show_row_names=T, show_column_names=FALSE, cluster_rows=T, cluster_columns=T, top_annotation=annot, col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), 
			cell_fun = function(j, i, x, y, width, height, fill) {
			  grid.text(round(cormat[i, j], 2), x, y, gp = gpar(fontsize = 10, col = "black"))
			})
print(ht)
dev.off()
# looks like there is limited correlation between the modalities
data$time_num= data$timepoint %>% as.numeric()
data$pre_post=ifelse(data$time_num<=0, "pre", "post")
table(data$pre_post,data$timepoint)
meanfunc=function(x){
	if(all(is.na(x))){
		return(NA)
	} else {
	   return(mean(x,na.rm=TRUE))
	}
}

#take average of pre timepoints as better representation of baseline
modalitystr="tact|anti|cyto|tpol|cellfreq"
nrow(data)
data_train=data %>%
	filter(pre_post=="pre") %>%
	group_by(subject_id) %>%
	select(matches(modalitystr)) %>%
	summarise_all(meanfunc)
colnames(data_train)
nrow(data_train)

# make outcome measures and add to data
data_outcome_anti_IgG_PT_anti_post=data %>%
	filter(time_num>=7) %>%
	group_by(subject_id) %>%
	select(IgG_PT_anti,IgG_PT_anti_corrected) %>%
	summarise_all(meanfunc)
data_outcome_anti_IgG_PT_anti_post$IgG_PT_anti_post=data_outcome_anti_IgG_PT_anti_post$IgG_PT_anti
data_outcome_anti_IgG_PT_anti_post$IgG_PT_anti_corrected_post=data_outcome_anti_IgG_PT_anti_post$IgG_PT_anti_corrected
data_outcome_anti_IgG_PT_anti_post = data_outcome_anti_IgG_PT_anti_post %>% select(subject_id,IgG_PT_anti_post,IgG_PT_anti_corrected_post)

data_outcome_Monocytes_cellfreq_post=data %>%
	filter(timepoint==1) %>%
	select(subject_id,Monocytes_cellfreq,Monocytes_cellfreq_corrected)
data_outcome_Monocytes_cellfreq_post$Monocytes_cellfreq_post=data_outcome_Monocytes_cellfreq_post$Monocytes_cellfreq
data_outcome_Monocytes_cellfreq_post$Monocytes_cellfreq_corrected_post=data_outcome_Monocytes_cellfreq_post$Monocytes_cellfreq_corrected
data_outcome_Monocytes_cellfreq_post = data_outcome_Monocytes_cellfreq_post %>% select(subject_id,Monocytes_cellfreq_post,Monocytes_cellfreq_corrected_post)

data_meta=data %>%
	select(subject_id,age,infancy_vac,biological_sex,ethnicity,race,year_of_birth,dataset) %>%
	distinct()
if(!all(data_train$subject_id %in% data_meta$subject_id)){
	log_warn("Not all subjects in data_train are in data_meta")
}
data_train=full_join(data_train,data_meta,by="subject_id") 
data_train=full_join(data_train,data_outcome_anti_IgG_PT_anti_post,by="subject_id")
data_train=full_join(data_train,data_outcome_Monocytes_cellfreq_post,by="subject_id")
# load and add gex data, this includes post for CCL3
data_gex=read_tsv("/user/leuven/342/vsc34263/myvibdata/challenge/data_pp/gex/gex_data.tsv")
data_train=full_join(data_train,data_gex,by="subject_id")
missingsubjects=setdiff(all_subjects,data_train$subject_id)
if(length(missingsubjects)>0){
	log_warn(glue_collapse("Missing subjects: {missingsubjects}"),sep="\n")
}
#writing data
data_train %>% write_tsv(glue("data_pp/{dataset}/training_data.tsv"))
data_train$IgG_PT_anti_corrected

# table(data_train$timepoint)
data %>% write_tsv(glue("data_pp/{dataset}/merged_data.tsv"))

# data_train %>% group_by(timepoint) %>% select(IgG_PT_anti) %>% summarise_all(mean,na.rm=T) %>% print(n=Inf)

#get output measures

dim(data_train)
dim(data)
