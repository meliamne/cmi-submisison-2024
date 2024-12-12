library(tidymodels)
library(glue)
library(modeltime)
library(tidyverse)
library(parsnip)
library(ggplot2)
library(recipes)
library(yardstick)
library(xgboost)
library(glmnet)
library(brulee)
# library(cubist)
# library(bartMachine)
library(kernlab)
library(rules)
library(baguette)
library(plsmod)
library(brulee)
library(shapviz)
library(vip)
library(logger)
library(vip)
source("metrics.R")
library(future)
# plan(multisession,workers=5)
plan(sequential)

# library(rsample)
# library(workflows)
set.seed(42)

bayes_iters=1000

data=read_tsv("/user/leuven/342/vsc34263/myvibdata/challenge/data_pp/merged/training_data.tsv")
challengedata=data %>% filter(dataset=="2023_dataset")

#remove uninformative columns
dropcols=c('specimen_id','dataset','timepoint',"actual_day_relative_to_boost","planned_day_relative_to_boost","specimen_type","visit","year_of_birth","date_of_boost","time_num","pre_post","race")
data$ethnicity=ifelse(data$ethnicity=="Unknown",NA,data$ethnicity)
data=data %>% select(!any_of(dropcols))
colnames(data)

nan_counts <- lapply(data, function(x) sum(is.nan(x)))
any(nan_counts > 0)
inf_counts <- lapply(data, function(x) sum(is.infinite(x)))
any(inf_counts > 0)
maxvals=lapply(data, function(x){max(x,na.rm=T)})
maxvals=maxvals[unlist(lapply(maxvals,is.numeric))]
sort(unlist(maxvals),decreasing = T) %>% head

outcomecols=c('IgG_PT_anti_corrected','Monocytes_cellfreq_corrected',"CCL3_gex_corrected")

findbest=function(data_train,predvar,challengedata){
	log_info(glue("Fitting model for {predvar}"))
	postpoint=glue("{predvar}_post")
	changevar=glue("{predvar}_change")
	data_train %>% print()
	colnames(data_train)
	data_train[[changevar]]=data_train[[postpoint]]-data_train[[predvar]]
	data_train=data_train %>% drop_na(all_of(c(changevar))) #drop rows with missing values, also gets rid of challenge data
	#get rid of post columns and subject_id to prevent data leakage and overfitting
	outcomevals=data_train[[postpoint]]
	data_train= data_train %>% select(!any_of(c('subject_id')))
	data_train=data_train %>% select(!matches("post"))
	log_info(glue("data has {nrow(data_train)} rows"))
	folds <- vfold_cv(data_train, v = 5)
	#xgboost workflow
	form=as.formula(glue("{predvar}_change ~ ."))
	normalize_recipe=recipe(form, data = data_train) %>%
	step_normalize(all_numeric_predictors()) %>%
	step_dummy(all_nominal_predictors()) %>%
	step_novel(all_nominal_predictors())
	summary(normalize_recipe) %>% print(n=Inf)
	workflow_tb=workflow() %>% add_model(boost_tree(mode="regression",engine = "xgboost")) %>% add_recipe(normalize_recipe)
	tuneflow = workflow() %>%
		add_model(boost_tree(mode="regression",engine = "xgboost",learn_rate = tune(),trees = tune(),tree_depth = tune())) %>%
			add_recipe(normalize_recipe)
	bayes_controls=control_bayes(no_improve = 10,parallel_over = "resamples",verbose = TRUE,save_pred = TRUE,allow_par = F) #parallization fails
	res=tune_bayes(tuneflow,resamples=folds,metrics=metric_set(rmse,spearman_cor_abs),iter=bayes_iters,initial=10,control=bayes_controls)
	res_met=collect_metrics(res)
	res_met %>% filter(.metric=="spearman_cor_abs") %>% arrange(desc(mean)) %>% print(n=Inf)
	res_met %>% filter(.metric=="rmse") %>% arrange((mean)) %>% print(n=Inf)
	bestmodel=select_best(res,metric="rmse")
	bestmodel$predvar=glue("{predvar}_change")
	bestflow = workflow() %>%
		add_model(boost_tree(mode="regression",engine = "xgboost",learn_rate = bestmodel$learn_rate,trees = bestmodel$trees,tree_depth = bestmodel$tree_depth)) %>%
			add_recipe(normalize_recipe)
	cv_control=control_grid(save_pred = TRUE)
	cv_fit=fit(bestflow,data=data_train) %>% fit_resamples(resamples=folds,metrics=metric_set(spearman_cor_abs,rmse),control=cv_control)
	bestpred=bind_rows(cv_fit$.predictions) %>% arrange(.row)
	predictions_vec=bestpred[['.pred']]+data_train[[predvar]]
	predictionsdf=data.frame(pred=predictions_vec,truth=outcomevals) %>% as_tibble()
	final_fit = fit(bestflow, data=data_train)
	importance = vi(extract_fit_parsnip(final_fit))
	#get top 20 features
	top20features=importance %>% arrange(desc(Importance)) %>% head(20) %>% pull(Variable)
	data_train=data_train %>% select(all_of(c(top20features,changevar)))
	#refit model with top 20 features
	normalize_recipe=recipe(form, data = data_train) %>%
	step_normalize(all_numeric_predictors()) %>%
	step_dummy(all_nominal_predictors()) %>%
	step_novel(all_nominal_predictors())
	summary(normalize_recipe) %>% print(n=Inf)
	workflow_tb=workflow() %>% add_model(boost_tree(mode="regression",engine = "xgboost")) %>% add_recipe(normalize_recipe)
	tuneflow = workflow() %>%
		add_model(boost_tree(mode="regression",engine = "xgboost",learn_rate = tune(),trees = tune(),tree_depth = tune())) %>%
			add_recipe(normalize_recipe)
	bayes_controls=control_bayes(no_improve = 20,parallel_over = "resamples",verbose = TRUE,save_pred = TRUE,allow_par = F) #parallization fails
	res=tune_bayes(tuneflow,resamples=folds,metrics=metric_set(rmse,spearman_cor_abs),iter=bayes_iters,initial=10,control=bayes_controls)
	res_met=collect_metrics(res)
	res_met %>% filter(.metric=="spearman_cor_abs") %>% arrange(desc(mean)) %>% print(n=Inf)
	res_met %>% filter(.metric=="rmse") %>% arrange((mean)) %>% print(n=Inf)
	bestmodel=select_best(res,metric="rmse")
	bestmodel$predvar=glue("{predvar}_change")
	bestflow = workflow() %>%
		add_model(boost_tree(mode="regression",engine = "xgboost",learn_rate = bestmodel$learn_rate,trees = bestmodel$trees,tree_depth = bestmodel$tree_depth)) %>%
			add_recipe(normalize_recipe)
	cv_control=control_grid(save_pred = TRUE)
	cv_fit=fit(bestflow,data=data_train) %>% fit_resamples(resamples=folds,metrics=metric_set(spearman_cor_abs,rmse),control=cv_control)
	bestpred=bind_rows(cv_fit$.predictions) %>% arrange(.row)
	predictions_vec=bestpred[['.pred']]+data_train[[predvar]]
	predictionsdf=data.frame(pred=predictions_vec,truth=outcomevals) %>% as_tibble()
	final_fit = fit(bestflow, data=data_train)
	# Prepare challenge data for prediction
	challenge_preds = predict(final_fit, new_data = challengedata)$.pred + challengedata[[predvar]]
	names(challenge_preds)=challengedata$subject_id
	return(list(bestmodel=bestmodel,bestpred=bestpred,predictionsdf=predictionsdf,cv_fit=cv_fit,importance=importance,challenge_preds=challenge_preds))
}
findbest_lfc=function(data_train,predvar,challengedata){
	log_info(glue("Fitting model for {predvar}"))
	postpoint=glue("{predvar}_post")
	changevar=glue("{predvar}_lfc")
	data_train %>% print()
	colnames(data_train)
	data_train[[changevar]]=log2(data_train[[postpoint]]/data_train[[predvar]])
	data_train=data_train %>% drop_na(all_of(c(changevar)))
	#get rid of post columns and subject_id to prevent data leakage and overfitting
	outcomevals=data_train[[postpoint]]
	data_train= data_train %>% select(!any_of(c('subject_id')))
	data_train=data_train %>% select(!matches("post"))
	log_info(glue("data has {nrow(data_train)} rows"))
	folds <- vfold_cv(data_train, v = 5)
	#xgboost workflow
	form=as.formula(glue("{predvar}_lfc ~ ."))
	normalize_recipe=recipe(form, data = data_train) %>%
	step_normalize(all_numeric_predictors()) %>%
	step_dummy(all_nominal_predictors()) %>%
	step_novel(all_nominal_predictors())
	summary(normalize_recipe) %>% print(n=Inf)

	workflow_tb=workflow() %>% add_model(boost_tree(mode="regression",engine = "xgboost")) %>% add_recipe(normalize_recipe)
	tuneflow = workflow() %>%
		add_model(boost_tree(mode="regression",engine = "xgboost",learn_rate = tune(),trees = tune(),tree_depth = tune())) %>%
			add_recipe(normalize_recipe)
	bayes_controls=control_bayes(no_improve = 20,parallel_over = "resamples",verbose = TRUE,save_pred = TRUE,allow_par = F) #parallization fails
	res=tune_bayes(tuneflow,resamples=folds,metrics=metric_set(rmse,spearman_cor_abs),iter=bayes_iters,initial=10,control=bayes_controls)
	res_met=collect_metrics(res)
	res_met %>% filter(.metric=="spearman_cor_abs") %>% arrange(desc(mean)) %>% print(n=Inf)
	res_met %>% filter(.metric=="rmse") %>% arrange((mean)) %>% print(n=Inf)
	bestmodel=select_best(res,metric="rmse")
	bestmodel$predvar=glue("{predvar}_lfc")
	bestflow = workflow() %>%
		add_model(boost_tree(mode="regression",engine = "xgboost",learn_rate = bestmodel$learn_rate,trees = bestmodel$trees,tree_depth = bestmodel$tree_depth)) %>%
			add_recipe(normalize_recipe)
	cv_control=control_grid(save_pred = TRUE)
	cv_fit=fit(bestflow,data=data_train) %>% fit_resamples(resamples=folds,metrics=metric_set(spearman_cor_abs,rmse),control=cv_control)
	bestpred=bind_rows(cv_fit$.predictions) %>% arrange(.row)
	predictions_vec=bestpred[['.pred']]
	predictionsdf=data.frame(pred=predictions_vec,truth=outcomevals) %>% as_tibble()
	final_fit = fit(bestflow, data=data_train)
	importance = vi(extract_fit_parsnip(final_fit))
	#get top 20 features
	top20features=importance %>% arrange(desc(Importance)) %>% head(20) %>% pull(Variable)
	data_train=data_train %>% select(all_of(c(top20features,changevar)))
	#refit model with top 20 features
	normalize_recipe=recipe(form, data = data_train) %>%
	step_normalize(all_numeric_predictors()) %>%
	step_dummy(all_nominal_predictors()) %>%
	step_novel(all_nominal_predictors())
	summary(normalize_recipe) %>% print(n=Inf)
	workflow_tb=workflow() %>% add_model(boost_tree(mode="regression",engine = "xgboost")) %>% add_recipe(normalize_recipe)
	tuneflow = workflow() %>%
		add_model(boost_tree(mode="regression",engine = "xgboost",learn_rate = tune(),trees = tune(),tree_depth = tune())) %>%
			add_recipe(normalize_recipe)
	bayes_controls=control_bayes(no_improve = 20,parallel_over = "resamples",verbose = TRUE,save_pred = TRUE,allow_par = F) #parallization fails
	res=tune_bayes(tuneflow,resamples=folds,metrics=metric_set(rmse,spearman_cor_abs),iter=bayes_iters,initial=10,control=bayes_controls)
	res_met=collect_metrics(res)
	res_met %>% filter(.metric=="spearman_cor_abs") %>% arrange(desc(mean)) %>% print(n=Inf)
	res_met %>% filter(.metric=="rmse") %>% arrange((mean)) %>% print(n=Inf)
	bestmodel=select_best(res,metric="rmse")
	bestmodel$predvar=glue("{predvar}_lfc")
	bestflow = workflow() %>%
		add_model(boost_tree(mode="regression",engine = "xgboost",learn_rate = bestmodel$learn_rate,trees = bestmodel$trees,tree_depth = bestmodel$tree_depth)) %>%
			add_recipe(normalize_recipe)
	cv_control=control_grid(save_pred = TRUE)
	cv_fit=fit(bestflow,data=data_train) %>% fit_resamples(resamples=folds,metrics=metric_set(spearman_cor_abs,rmse),control=cv_control)
	bestpred=bind_rows(cv_fit$.predictions) %>% arrange(.row)
	predictions_vec=bestpred[['.pred']]
	predictionsdf=data.frame(pred=predictions_vec,truth=outcomevals) %>% as_tibble()
	final_fit = fit(bestflow, data=data_train)
	challenge_preds = predict(final_fit, new_data = challengedata)$.pred
	names(challenge_preds)=challengedata$subject_id
	return(list(bestmodel=bestmodel,bestpred=bestpred,predictionsdf=predictionsdf,cv_fit=cv_fit,importance=importance,challenge_preds=challenge_preds))
}

bestmodellist=list()
predlist_1=list()
predictiondflist=list()
cvfitlist=list()
importancelist=list()
predlist=list()
for(predvar in outcomecols){
	res=findbest(data,predvar,challengedata)
	bestmodellist[[predvar]]=res$bestmodel
	predlist_1[[predvar]]=res$bestpred
	predictiondflist[[predvar]]=res$predictionsdf
	cvfitlist[[predvar]]=res$cv_fit
	importancelist[[predvar]]=res$importance
	predlist[[predvar]]=res$challenge_preds
}
save.image(file = "changevar_environment.RData")

bestmodellist_lfc=list()
predlist_1_lfc=list()
predictiondflist_lfc=list()
cvfitlist_lfc=list()
importancelist_lfc=list()
predlist_lfc=list()
for(predvar in outcomecols){
	res=findbest_lfc(data,predvar,challengedata)
	bestmodellist_lfc[[predvar]]=res$bestmodel
	predlist_1_lfc[[predvar]]=res$bestpred
	predictiondflist_lfc[[predvar]]=res$predictionsdf
	cvfitlist_lfc[[predvar]]=res$cv_fit
	importancelist_lfc[[predvar]]=res$importance
	predlist_lfc[[predvar]]=res$challenge_preds
}


#prep submission
submission=read_tsv("/user/leuven/342/vsc34263/myvibdata/challenge/3rdChallengeSubmissionTemplate_10032024.tsv")
colnames(submission)
submission[['1.1) IgG-PT-D14-titer-Rank']] = order(predlist[['IgG_PT_anti_corrected']][as.character(submission$SubjectID)])
submission[['1.2) IgG-PT-D14-FC-Rank']] = order(predlist_lfc[['IgG_PT_anti_corrected']][as.character(submission$SubjectID)])
submission[['2.1) Monocytes-D1-Rank']] = order(predlist[['Monocytes_cellfreq_corrected']][as.character(submission$SubjectID)])
submission[['2.2) Monocytes-D1-FC-Rank']] = order(predlist_lfc[['Monocytes_cellfreq_corrected']][as.character(submission$SubjectID)])
submission[['3.1) CCL3-D3-Rank']] = order(predlist[['CCL3_gex_corrected']][as.character(submission$SubjectID)])
submission[['3.2) CCL3-D3-FC-Rank']] = order(predlist_lfc[['CCL3_gex_corrected']][as.character(submission$SubjectID)])

write_tsv(submission, "/user/leuven/342/vsc34263/myvibdata/challenge/submission.tsv")
