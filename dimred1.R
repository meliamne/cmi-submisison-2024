library(tidyverse)
library(ggplot2)
library(glue)
library(ComplexHeatmap)
library(pcaMethods)
library(logger)
library(circlize)
#challengefiles contain relatively more cyto and tact and tpol data
challengefilelocs=list(
	cellfreq="/user/leuven/342/vsc34263/myvibdata/challenge/data/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/challenge_pbmc_cell_frequency_wide.tsv",
	anti="/user/leuven/342/vsc34263/myvibdata/challenge/data/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/challenge_plasma_antibody_levels_wide.tsv",
	# metadata="data/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/challenge_subject_specimen.tsv",
	cyto="/user/leuven/342/vsc34263/myvibdata/challenge/data/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/challenge_plasma_cytokine_concentrations_by_legendplex_wide.tsv",
	tact="/user/leuven/342/vsc34263/myvibdata/challenge/data/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/challenge_t_cell_activation_wide.tsv",
	tpol="/user/leuven/342/vsc34263/myvibdata/challenge/data/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/challenge_t_cell_polarization_wide.tsv"
)
# lapply(challengefilelocs,function(x){print(dim(read_tsv(x)))})


filelocs=list(
	cellfreq="/user/leuven/342/vsc34263/myvibdata/challenge/data/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/training_pbmc_cell_frequency_wide.tsv",
	anti="/user/leuven/342/vsc34263/myvibdata/challenge/data/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/training_plasma_antibody_levels_wide.tsv",
	# metadata="data/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/training_subject_specimen.tsv",
	cyto="/user/leuven/342/vsc34263/myvibdata/challenge/data/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/training_plasma_cytokine_concentrations_by_legendplex_wide.tsv",
	tact="/user/leuven/342/vsc34263/myvibdata/challenge/data/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/training_t_cell_activation_wide.tsv",
	tpol="/user/leuven/342/vsc34263/myvibdata/challenge/data/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/training_t_cell_polarization_wide.tsv"
)

# lapply(filelocs,function(x){print(dim(read_tsv(x)))})
data_meta_challenge=read_tsv("data/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/challenge_subject_specimen.tsv")
data_meta_challenge$challenge=TRUE
data_meta=read_tsv("data/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/training_subject_specimen.tsv")
data_meta$challenge=FALSE


data_meta=rbind(data_meta,data_meta_challenge)
data_meta$dataset_timepoint=paste(data_meta$dataset,data_meta$timepoint,sep="_")
table(data_meta$dataset_timepoint,useNA="always")

table(data_meta$timepoint,data_meta$dataset)
table(data_meta$timepoint,data_meta$challenge)
all_subjects = unique(data_meta$subject_id)
all_specimens = unique(data_meta$specimen_id)


cormatsize=1200
ndim=15
keeppca=10
maxna_col=0.8


metacols=c("specimen_id","subject_id","dataset","timepoint","challenge","dataset_timepoint")
metacols_drop=metacols[-1]



pptdata=function(dataset,dataloc,dataloc_challenge,meta,keeppca){
	print(glue("Processing {dataset}"))
	dir.create(glue("./plots/{dataset}"), showWarnings = FALSE, recursive = TRUE)
	dir.create(glue("./data_pp/{dataset}"), showWarnings = FALSE, recursive = TRUE)
	data=read_tsv(dataloc)
	if(!all(data$specimen_id %in% meta$specimen_id)){
		log_warn(glue("Not all specimen_ids in {dataset} are in the metadata"))
	}
	data_challenge=read_tsv(dataloc_challenge)
	if(!all(data_challenge$specimen_id %in% meta$specimen_id)){
		log_warn(glue("Not all specimen_ids in {dataset} are in the metadata"))
	}
	datacols=colnames(data)[-1]
	meta_sub=meta %>% select(all_of(metacols))
	data=data %>%
		inner_join(meta_sub,by=c("specimen_id"="specimen_id"))
	data_challenge=data_challenge %>%
		inner_join(meta_sub,by=c("specimen_id"="specimen_id"))
	if(!(all(colnames(data)==colnames(data_challenge)))){log_warn(glue("Column names in {dataset} and challenge data do not match"))}
	data=full_join(data,data_challenge)
	# data=data %>% filter(timepoint>=0) #filter out early time points due to high missingness in training data
	sdzero= data %>% select(!c(metacols)) %>% summarise_all(sd,na.rm=T) %>% unlist()
	sdzero=sdzero==0
	if(any(sdzero)){
		badcols=names(sdzero)[sdzero]
		log_warn(glue_collapse("Columns with zero standard deviation: {badcols} in {dataset} being removed"),sep="\n")
		data=data %>% select(!all_of(badcols))
	}
	nacols=lapply(data,function(x){mean(is.na(x))})
	nadf=tibble(colnames=names(unlist(nacols)),nanfrac=unlist(nacols)) %>% arrange(desc(nanfrac))
	nadf %>% filter(nanfrac>0) %>% print(n=Inf)
	narows=apply(data,1,function(x){mean(is.na(x))})
	narowsdf=tibble(rownames=rownames(data),nanfrac=narows) %>% arrange(desc(nanfrac))
	narowsdf %>% filter(nanfrac>0) %>% print(n=Inf)
	dropnacols=names(nadf)[nadf$nanfrac>maxna_col]
	data=data %>% select(!all_of(dropnacols))
	numberofnotna=data %>% select(!c(metacols)) %>% apply(2,function(x) sum(!is.na(x))) %>% as.character()
	annot = HeatmapAnnotation(obsnum=anno_text(numberofnotna))
	mat=data %>% select(!c(metacols)) %>% as.matrix()
	cormat=cor(mat,use="pairwise.complete.obs")
	png(file=glue("plots/{dataset}/correlation_heatmap.png"),width=cormatsize,height=cormatsize)
	ht=Heatmap(cormat, name="correlation_raw", show_row_names=T, show_column_names=FALSE, cluster_rows=T, cluster_columns=T, top_annotation=annot, col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))
	print(ht)
	dev.off()
	dimred=pca(mat,nPcs=ncol(mat)-1,method="bpca",center=TRUE,completeObs = T,scale="vector",maxSteps =10000) #-1 because cellfreqs are not independent
	plot=ggplot(data.frame(x=1:length(dimred@R2),y=dimred@R2),aes(x=x,y=y))+geom_point()+geom_line()
	ggsave(glue("plots/{dataset}/pca.png"),plot)
	dimreddata=dimred@scores[,1:min(keeppca,ncol(dimred@scores))]
	dimreddata = cbind(data$specimen_id,dimreddata) %>% as_tibble() 
	colnames(dimreddata)=c("specimen_id",glue("PC{1:(ncol(dimreddata)-1)}"))
	data_sub=select(data,!any_of(metacols_drop))
	outdata=left_join(dimreddata,data_sub,by="specimen_id")
	outdata %>% write_tsv(glue("data_pp/{dataset}/pca_scores.tsv"))
	imputed = completeObs(dimred)
	pca=princomp(imputed)
	png(file=glue("plots/{dataset}/pca_imputed.png"),width=800,height=800)
	plot(pca,type="l")
	dev.off()
	dimreddata=dimred@scores[,1:min(keeppca,ncol(dimred@scores))]
	dimreddata = cbind(data$specimen_id,imputed,dimreddata) %>% as_tibble() %>% write_tsv(glue("data_pp/{dataset}/pca_scores_imputed.tsv"))
	cormat=cbind(mat,dimreddata) %>% cor(use="pairwise.complete.obs")
	png(file=glue("plots/{dataset}/correlation_heatmap_pca.png"),width=cormatsize,height=cormatsize)
	ht=Heatmap(cormat, name="correlation_pca", show_row_names=T, show_column_names=FALSE, cluster_rows=T, cluster_columns=T, col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))
	print(ht)
	dev.off()

	# now make batch corrected data
	data_corrected=matrix(NA,nrow=nrow(mat),ncol=ncol(mat))
	for(set_cor in unique(data$dataset_timepoint)){
		idx=which(data$dataset_timepoint==set_cor)
		data_corrected[idx,]=apply(mat[idx,],2,scale)
	}
	data_corrected= data_corrected %>% as_tibble()
	dropnarows=which(apply(data_corrected,1,function(x){mean(is.na(x))})==1)
	data_corrected=cbind(data$specimen_id,data_corrected)
	colnames(data_corrected)=c("specimen_id",colnames(mat))
	table(data$dataset_timepoint,useNA="always")
	nacols=lapply(data_corrected,function(x){mean(is.na(x),na.rm=T)})
	dropnacols=names(nacols)[nacols>maxna_col]
	data_corrected=data_corrected %>% select(!all_of(dropnacols))
	if(length(dropnarows)>0){data_corrected=data_corrected[-dropnarows,]}
	# data_corrected %>% head
	#
	mat=data_corrected %>% select(!c("specimen_id")) %>% as.matrix()
	mat[is.nan(mat)]=NA
	# apply(mat,1,function(x){mean(is.na(x))}) %>% sort()
	cormat=cor(mat,use="pairwise.complete.obs")
	numberofnotna=data_corrected %>% select(!c("specimen_id")) %>% apply(2,function(x) sum(!is.na(x))) %>% as.character()
	annot = HeatmapAnnotation(obsnum=anno_text(numberofnotna))
	png(file=glue("plots/{dataset}/correlation_heatmap_corrected.png"),width=cormatsize,height=cormatsize)
	ht=Heatmap(cormat, name="correlation_raw", show_row_names=T, show_column_names=FALSE, cluster_rows=T, cluster_columns=T, top_annotation=annot, col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))
	print(ht)
	dev.off()
	dimred=pca(mat,nPcs=ncol(mat)-1,method="bpca",center=TRUE,completeObs = T,scale="vector",maxSteps =10000) #-1 because cellfreqs are not independent
	plot=ggplot(data.frame(x=1:length(dimred@R2),y=dimred@R2),aes(x=x,y=y))+geom_point()+geom_line()
	ggsave(glue("plots/{dataset}/pca_corrected.png"),plot)
	dimreddata=dimred@scores[,1:min(keeppca,ncol(dimred@scores))]
	dimreddata = cbind(data_corrected$specimen_id,dimreddata) %>% as_tibble() 
	colnames(dimreddata)=c("specimen_id",glue("PC{1:(ncol(dimreddata)-1)}"))
	data_corrected_sub=select(data_corrected,!any_of(metacols_drop))
	outdata=left_join(dimreddata,data_corrected_sub,by="specimen_id")
	outdata %>% write_tsv(glue("data_pp/{dataset}/pca_scores_corrected.tsv"))
	imputed = completeObs(dimred)
	pca=princomp(imputed)
	png(file=glue("plots/{dataset}/pca_imputed_corrected.png"),width=800,height=800)
	plot(pca,type="l")
	dev.off()
	dimreddata=dimred@scores[,1:min(keeppca,ncol(dimred@scores))]
	dimreddata = cbind(data$specimen_id,imputed,dimreddata) %>% as_tibble() %>% write_tsv(glue("data_pp/{dataset}/pca_scores_imputed_corrected.tsv"))
	cormat=cbind(mat,dimreddata) %>% cor(use="pairwise.complete.obs")
	png(file=glue("plots/{dataset}/correlation_heatmap_pca_corrected.png"),width=cormatsize,height=cormatsize)
	ht=Heatmap(cormat, name="correlation_pca", show_row_names=T, show_column_names=FALSE, cluster_rows=T, cluster_columns=T, col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))
	print(ht)
	dev.off()
}


#old version
# pptdata=function(dataset,dataloc,meta,keeppca){
# 	print(glue("Processing {dataset}"))
# 	dir.create(glue("./plots/{dataset}"), showWarnings = FALSE, recursive = TRUE)
# 	dir.create(glue("./data_pp/{dataset}"), showWarnings = FALSE, recursive = TRUE)
# 	data=read_tsv(dataloc)
# 	datacols=colnames(data)[-1]
# 	meta_sub=meta %>% select(metacols)
# 	data=data %>%
# 		left_join(meta_sub,by=c("specimen_id"="specimen_id"))
# 	data=data %>% filter(timepoint>=0) #filter out early time points due to high missingness in training data
# 	sdzero= data %>% select(!c(metacols)) %>% summarise_all(sd,na.rm=T) %>% unlist()
# 	sdzero=sdzero==0
# 	if(any(sdzero)){
# 		badcols=names(sdzero)[sdzero]
# 		log_warn(glue_collapse("Columns with zero standard deviation: {badcols} in {dataset} being removed"),sep="\n")
# 		data=data %>% select(!all_of(badcols))
# 	}
# 	nacols=lapply(data,function(x){mean(is.na(x))})
# 	nadf=tibble(colnames=names(unlist(nacols)),nanfrac=unlist(nacols)) %>% arrange(desc(nanfrac))
# 	nadf %>% filter(nanfrac>0) %>% print(n=Inf)
# 	narows=apply(data,1,function(x){mean(is.na(x))})
# 	narowsdf=tibble(rownames=rownames(data),nanfrac=narows) %>% arrange(desc(nanfrac))
# 	narowsdf %>% filter(nanfrac>0) %>% print(n=Inf)
# 	table(data$dataset)
# 	numberofnotna=data %>% select(!c(metacols)) %>% apply(2,function(x) sum(!is.na(x))) %>% as.character()
# 	annot = HeatmapAnnotation(obsnum=anno_text(numberofnotna))
# 	mat=data %>% select(!c(metacols)) %>% as.matrix()
# 	cormat=cor(mat,use="pairwise.complete.obs")
# 	png(file=glue("plots/{dataset}/correlation_heatmap.png"),width=cormatsize,height=cormatsize)
# 	ht=Heatmap(cormat, name="correlation_raw", show_row_names=T, show_column_names=FALSE, cluster_rows=T, cluster_columns=T, top_annotation=annot, col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))
# 	print(ht)
# 	dev.off()
# 	dimred=pca(mat,nPcs=ncol(mat)-1,method="bpca",center=TRUE,completeObs = T,scale="vector",maxSteps =10000) #-1 because cellfreqs are not independent
# 	plot=ggplot(data.frame(x=1:length(dimred@R2),y=dimred@R2),aes(x=x,y=y))+geom_point()+geom_line()
# 	ggsave(glue("plots/{dataset}/pca.png"),plot)
# 	dimreddata=dimred@scores[,1:min(keeppca,ncol(dimred@scores))]
# 	dimreddata = cbind(data$specimen_id,dimreddata) %>% as_tibble() 
# 	colnames(dimreddata)=c("specimen_id",glue("PC{1:(ncol(dimreddata)-1)}"))
# 	data_sub=select(data,!c(subject_id,dataset,timepoint))
# 	outdata=left_join(dimreddata,data_sub,by="specimen_id")
# 	outdata %>% write_tsv(glue("data_pp/{dataset}/pca_scores.tsv"))
# 	cormat=cbind(mat,dimreddata) %>% cor(use="pairwise.complete.obs")
# 	png(file=glue("plots/{dataset}/correlation_heatmap_pca.png"),width=cormatsize,height=cormatsize)
# 	ht=Heatmap(cormat, name="correlation_pca", show_row_names=T, show_column_names=FALSE, cluster_rows=T, cluster_columns=T, col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))
# 	print(ht)
# 	dev.off()
# }


keeppcalist=list(
	cellfreq=13,
	anti=10,
	cyto=3,
	tact=10,
	tpol=3
)
for(i in 1:length(filelocs)){
	dataset=names(filelocs)[[i]]
	pptdata(dataset,filelocs[[dataset]],challengefilelocs[[dataset]],data_meta,keeppcalist[[dataset]])
}

print("Done")