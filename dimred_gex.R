library(tidyverse)
library(DESeq2)
library(org.Hs.eg.db)
chosengenes=c("CCL3")
chosengenes_id=c('ENSG00000277632.1')
# gex=read_tsv("/user/leuven/342/vsc34263/myvibdata/challenge/data/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/training_pbmc_gene_expression_wide_raw_count.tsv")
data_meta_challenge=read_tsv("data/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/challenge_subject_specimen.tsv")
data_meta_challenge$challenge=TRUE
data_meta=read_tsv("data/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/training_subject_specimen.tsv")
data_meta$challenge=FALSE
meta=rbind(data_meta,data_meta_challenge)
metacols=colnames(meta)
rownames(meta)=meta$specimen_id %>% as.character
gex=read.csv("/user/leuven/342/vsc34263/myvibdata/challenge/data/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/training_pbmc_gene_expression_wide_raw_count.tsv", sep="\t", header=TRUE)
gex_challenge=read.csv("/user/leuven/342/vsc34263/myvibdata/challenge/data/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/challenge_pbmc_gene_expression_wide_raw_count.tsv", sep="\t", header=TRUE)
# gex=rbind(gex,gex_challenge)

subjects=gex$specimen_id
gex=gex[,-1]
gex=t(gex)
colnames(gex) = subjects

meta_full=meta[as.character(subjects),]
nrow(meta_full)
ncol(gex)

meta_full$prepost=case_when(
  meta_full$timepoint == 0 ~ "pre",
  meta_full$timepoint == 3 ~ "post",
  T ~ NA
)
meta_full$prepost=as.factor(meta_full$prepost) %>% relevel(ref="pre")
meta_full$subject_id = meta_full$subject_id %>% as.character %>% as.factor
keeppbool=which(!is.na(meta_full$prepost))
meta_full=meta_full[keeppbool,]
gex=gex[,keeppbool]
#fit model
dds <- DESeqDataSetFromMatrix(
  countData = gex,
  colData = meta_full,
  design = ~ 1
)
dds <- estimateSizeFactors(dds)
gex_normalized <- counts(dds, normalized = TRUE)
if(!all(rownames(gex_normalized) == rownames(gex))){
  stop("rownames of normalized counts do not match input")
}
variances=apply(gex_normalized, 1, var)
varbool=order(variances, decreasing=T)[1:5000]
gex = gex[varbool,]
meta_full$CCL3=gex_normalized['ENSG00000277632.1',] %>% scale()
dds <- DESeqDataSetFromMatrix(
  countData = gex,
  colData = meta_full,
  design = ~ prepost + subject_id
)
dds <- DESeq(dds)
res <- results(dds, contrast=c("prepost", "post", "pre"))
res <- res[order(res$padj),]
res$gene = rownames(res)
res$gene_id= rownames(res)
res$gene = res$gene %>% str_split("\\.") %>% map_chr(1)
res$gene = mapIds(org.Hs.eg.db, keys = res$gene, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$gene[is.na(res$gene)] = paste0("Unknown_", 1:sum(is.na(res$gene)))
res$gene = factor(res$gene)
res_tib=res %>% as_tibble()
res_tib %>% print(n=100)
res_tib %>% filter(str_detect(gene, "CCL")) %>% print(n=100)
chosengenes=c(chosengenes,res_tib %>% head(20) %>% pull(gene))
chosengenes_id=c(chosengenes_id,res_tib %>% head(20) %>% pull(gene_id))
dds <- DESeqDataSetFromMatrix(
  countData = gex,
  colData = meta_full,
  design = ~ CCL3 + dataset
)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name = "CCL3")
res <- res[order(res$padj),]
res$gene = rownames(res)
res$gene_id= rownames(res)
res$gene = res$gene %>% str_split("\\.") %>% map_chr(1)
res$gene = mapIds(org.Hs.eg.db, keys = res$gene, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$gene[is.na(res$gene)] = paste0("Unknown_", 1:sum(is.na(res$gene)))
res$gene = factor(res$gene)
res_tib=res %>% as_tibble()
res_tib %>% print(n=100)
chosengenes=c(chosengenes,res_tib %>% head(11) %>% pull(gene))
chosengenes_id=c(chosengenes_id,res_tib %>% head(11) %>% pull(gene_id))
res_tib %>% filter(str_detect(gene, "CCL")) %>% print(n=100)

gex=read.csv("/user/leuven/342/vsc34263/myvibdata/challenge/data/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/training_pbmc_gene_expression_wide_raw_count.tsv", sep="\t", header=TRUE)
gex_challenge=read.csv("/user/leuven/342/vsc34263/myvibdata/challenge/data/cmipb_challenge_datasets/current/3rd_challenge/harmonized_and_processed_data/master_harmonized_data_TSV/challenge_pbmc_gene_expression_wide_raw_count.tsv", sep="\t", header=TRUE)
gex=rbind(gex,gex_challenge)

subjects=gex$specimen_id
gex=gex[,-1]
gex=t(gex)
colnames(gex) = subjects

# chosengenes_id=c('ENSG00000277632.1',rownames(gex)[1:10])

meta_full=meta[as.character(subjects),]
nrow(meta_full)
ncol(gex)

meta_full$prepost=case_when(
  meta_full$timepoint == 0 ~ "pre",
  meta_full$timepoint > 0 ~ "post",
  T ~ NA
)
meta_full$prepost=as.factor(meta_full$prepost) %>% relevel(ref="pre")
meta_full$subject_id = meta_full$subject_id %>% as.character %>% as.factor
keeppbool=which(!is.na(meta_full$prepost))
meta_full=meta_full[keeppbool,]
gex=gex[,keeppbool]
#fit model
dds <- DESeqDataSetFromMatrix(
  countData = gex,
  colData = meta_full,
  design = ~ 1
)
dds <- estimateSizeFactors(dds)
gex_normalized <- counts(dds, normalized = TRUE)
gex_normalized = gex_normalized[chosengenes_id,]

geneids=rownames(gex_normalized) %>% str_split("\\.") %>% map_chr(1)
genanames=mapIds(org.Hs.eg.db, keys = geneids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
genanames[is.na(genanames)] = paste0("Unknown_", 1:sum(is.na(genanames)))
rownames(gex_normalized)=paste0(genanames,"_gex")

meta_full$dataset_timepoint=paste(meta_full$dataset,meta_full$timepoint,sep="_")
table(meta_full$dataset_timepoint)

data_corrected=matrix(NA,nrow=nrow(gex_normalized),ncol=ncol(gex_normalized))
for(set_cor in unique(meta_full$dataset_timepoint)){
	idx=which(meta_full$dataset_timepoint==set_cor)
	data_corrected[,idx]=apply(gex_normalized[,idx],1,scale)
}
rownames(data_corrected)=paste0(rownames(gex_normalized),"_corrected")
colnames(data_corrected)=colnames(gex_normalized)

gex_merge=rbind(gex_normalized,data_corrected)
gex_merge_df= gex_merge %>% t %>% as_tibble()
metacols=colnames(meta_full)

gex_merge_df=bind_cols(meta_full,gex_merge_df)

gex_merge_df_pre=gex_merge_df %>% filter(prepost=="pre")
gex_merge_df_post=gex_merge_df %>% filter(prepost=="post")

data_outcome_CCL3_post=gex_merge_df_post %>%
	dplyr::select(all_of(c("subject_id","CCL3_gex","CCL3_gex_corrected")))
data_outcome_CCL3_post$CCL3_gex_post=data_outcome_CCL3_post[['CCL3_gex']]
data_outcome_CCL3_post$CCL3_gex_corrected_post=data_outcome_CCL3_post[['CCL3_gex_corrected']]
data_outcome_CCL3_post = data_outcome_CCL3_post %>% dplyr::select(subject_id,CCL3_gex_post,CCL3_gex_corrected_post)


outdata=full_join(gex_merge_df_pre,data_outcome_CCL3_post,by="subject_id")
subjectids=outdata$subject_id
outdata <- outdata %>% dplyr::select(-all_of(metacols))
outdata$subject_id=subjectids

outdata %>% write_tsv("data_pp/gex/gex_data.tsv")