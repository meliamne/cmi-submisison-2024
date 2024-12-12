library(tidyverse)
library(rsample)
library(tidymodels)
library(ggplot2)

library(glue)
library(RColorBrewer)
library(logger)
library(pals)
library(shades)


polycol=pals::polychrome() %>% as.vector() %>% unname()
scale_colour_gradient_hot=scale_colour_gradientn(name = "", colors = rev(brewer.pal(11,"Spectral")))
scale_fill_gradient_hot=scale_fill_gradientn(name = "", colors = rev(brewer.pal(11,"Spectral")))

newdark_theme <- function() {
        # DarkTheme() +
        theme_dark() +
                theme(
                        axis.line = element_line(color = "white"),
                        axis.text = element_text(color = "white"),
                        axis.title = element_text(color = "white"),
                        legend.text = element_text(color = "white"),
                        legend.title = element_text(color = "white"),
                        legend.background = element_rect(fill = "#1f1f1f", color = "#1f1f1f"),
                        strip.text = element_text(color = "white"),
                        plot.background = element_rect(fill = "#1f1f1f",colour = "#1f1f1f"),
                        panel.background = element_rect(fill = "#1f1f1f", color = "#1f1f1f"),
                        plot.title = element_text(color = "white"),
                        panel.border = element_rect(fill = NA, color = "#1f1f1f"),
                        legend.key = element_rect(fill = "#1f1f1f",color = "#1f1f1f"),
                        rect=element_rect(fill = "#1f1f1f", color = "#1f1f1f"),
                )
}

data=read_tsv("/user/leuven/342/vsc34263/myvibdata/challenge/data_pp/merged/merged_data.tsv")
data$race=case_when(
	data$race=="White" ~ "White",
	data$race=="Asian" ~ "Asian",
	TRUE ~ NA
)
data$ethnicity=ifelse(data$ethnicity=="Unknown",NA,data$ethnicity)
data$timepoint_day_clean=ifelse(data$actual_day_relative_to_boost<=0,0,data$actual_day_relative_to_boost)
data$subject_id=factor(data$subject_id)
colnames(data)
data$PT_ratio=data$PT_P01579_tpol/data$PT_P05113_tpol

table(!(is.na(data$PT_ratio)),data$timepoint_day_clean)
table(!(is.na(data$PT_ratio)),data$subject_id)
data %>% filter(subject_id=="80") %>% select(subject_id,timepoint_day_clean,PT_ratio)
test=table(data$PT_ratio,data$timepoint_day_clean,data$subject_id)!=0
apply(test,3,sum)
data_agr %>% filter(subject_id=="80") %>% select(subject_id,timepoint_day_clean,timepoint_day_4step,PT_ratio)


split=initial_split(data,prop=0.8)
train=training(split)
test=testing(split)
folds <- vfold_cv(train, v = 5)



table(data$timepoint)
table(data$actual_day_relative_to_boost)
data$timepoint_day_clean=ifelse(data$actual_day_relative_to_boost<=0,0,data$actual_day_relative_to_boost)
data$actual_day_relative_to_boost_factor=factor(data$timepoint_day_clean)

data_agr

targetcols=c("IgG_PT_anti","Monocytes_cellfreq","PT_ratio")
targetcol=targetcols[1]

width=15

for(targetcol in targetcols){
	data_tmp = data %>% filter(!is.na(!!sym(targetcol)))
	basevalcol=glue("{targetcol}_base")
	changecol=glue("{targetcol}_change")
	changefrombaselinecol=glue("{targetcol}_changefrombaseline")
	fcchangefrombaselinecol=glue("{targetcol}_fcchangefrombaseline")
	lfcchangefrombaselinecol=glue("{targetcol}_lfcchangefrombaseline")
	outcomecol=glue("{targetcol}_at14")

	#boxplots over time
	plot = ggplot(data_tmp, aes(x = actual_day_relative_to_boost_factor, y = !!sym(targetcol))) +
		geom_boxplot() +
		newdark_theme() +
		labs(title = glue("{targetcol} over time"), x = "Days relative to boost", y = glue("{targetcol}")) +
		theme(legend.position = "none",
			panel.grid.major = element_line(colour = "grey80"),
			panel.grid.minor = element_blank()) +
		scale_x_discrete(breaks = c(0, 1, 3, 7, 14, 30))
	ggsave(glue("/user/leuven/342/vsc34263/myvibdata/challenge/plots/time/{targetcol}_over_time_boxplot.png"), plot, width = width, height = 10)
	plot = ggplot(data_tmp, aes(x = actual_day_relative_to_boost_factor, y = !!sym(targetcol))) +
		geom_boxplot() +
		newdark_theme() +
		labs(title = glue("{targetcol} over time"), x = "Days relative to boost", y = glue("{targetcol}")) +
		theme(legend.position = "none",
			panel.grid.major = element_line(colour = "grey80"),
			panel.grid.minor = element_blank()) +
		scale_x_discrete(breaks = c(0, 1, 3, 7, 14, 30))+
		facet_wrap(~dataset)
	ggsave(glue("/user/leuven/342/vsc34263/myvibdata/challenge/plots/time/{targetcol}_over_time_boxplot_wrap.png"), plot, width = width, height = 10)

	data_dist=data_tmp %>% 
		filter(timepoint_day_clean==0) %>%
		group_by(dataset) %>%
		summarise(mean=mean(!!sym(targetcol),na.rm=T),sd=sd(!!sym(targetcol),na.rm=T),n=n(),.groups="drop")
	data_dist



	#plot participants over time
	plot = ggplot(data_tmp, aes(x = timepoint_day_clean, y = !!sym(targetcol), color = subject_id)) +
		geom_point() +
		geom_line() +
		newdark_theme() +
		geom_smooth(colour = "black") +
		labs(title = glue("{targetcol} over time"), x = "Days relative to boost", y = glue("{targetcol}")) +
		theme(legend.position = "none",
			panel.grid.major = element_line(colour = "grey80"),
			panel.grid.minor = element_blank()) +
		scale_x_continuous(breaks = c(0, 1, 3, 7, 14, 30),limits=c(0,32))
	ggsave(glue("/user/leuven/342/vsc34263/myvibdata/challenge/plots/time/{targetcol}_over_time.png"), plot, width = width, height = 10)
	plot = ggplot(data_tmp, aes(x = timepoint_day_clean, y = !!sym(targetcol), color = subject_id)) +
		geom_point() +
		geom_line() +
		newdark_theme() +
		geom_smooth(colour = "black") +
		labs(title = glue("{targetcol} over time"), x = "Days relative to boost", y = glue("{targetcol}")) +
		theme(legend.position = "none",
			panel.grid.major = element_line(colour = "grey80"),
			panel.grid.minor = element_blank()) +
		scale_x_continuous(breaks = c(0, 1, 3, 7, 14, 30),limits=c(0,32))+
		facet_wrap(~dataset)
	ggsave(glue("/user/leuven/342/vsc34263/myvibdata/challenge/plots/time/{targetcol}_over_time_wrap.png"), plot, width = width, height = 10)


	data_baseline=data_tmp %>% filter(timepoint_day_clean==0)
	data_tmp[[basevalcol]]=data_baseline[[targetcol]][match(data_tmp$subject_id,data_baseline$subject_id)]
	data_outcome=data_tmp %>% filter(timepoint_day_clean==14)
	data_tmp[[outcomecol]]=data_outcome[[targetcol]][match(data_tmp$subject_id,data_outcome$subject_id)]
	data_tmp[[changefrombaselinecol]]=data_tmp[[targetcol]]-data_tmp[[basevalcol]]

	plotcol=glue("{targetcol}_changefrombaseline")
	plot = ggplot(data_tmp, aes(x = timepoint_day_clean, y = !!sym(changefrombaselinecol), color = subject_id)) +
		geom_point() +
		geom_line() +
		newdark_theme() +
		geom_smooth(colour = "black") +
		labs(title = glue("{plotcol} over time"), x = "Days relative to boost", y = glue("{plotcol}")) +
		theme(legend.position = "none",
			panel.grid.major = element_line(colour = "grey80"),
			panel.grid.minor = element_blank()) +
		scale_x_continuous(breaks = c(0, 1, 3, 7, 14, 30),limits=c(0,32))
	ggsave(glue("/user/leuven/342/vsc34263/myvibdata/challenge/plots/time/{targetcol}_changeover_time.png"), plot, width = width, height = 10)

	data_tmp[[fcchangefrombaselinecol]]=(data_tmp[[targetcol]]/data_tmp[[basevalcol]])

	plotcol=glue("{targetcol}_fcchangefrombaseline")
	plot = ggplot(data_tmp, aes(x = timepoint_day_clean, y = !!sym(fcchangefrombaselinecol), color = subject_id)) +
		geom_point() +
		geom_line() +
		newdark_theme() +
		geom_smooth(colour = "black") +
		labs(title = glue("{plotcol} over time"), x = "Days relative to boost", y = glue("{plotcol}")) +
		theme(legend.position = "none",
			panel.grid.major = element_line(colour = "grey80"),
			panel.grid.minor = element_blank()) +
		scale_x_continuous(breaks = c(0, 1, 3, 7, 14, 30),limits=c(0,32))
	ggsave(glue("/user/leuven/342/vsc34263/myvibdata/challenge/plots/time/{targetcol}_fcchangeover_time.png"), plot, width = width, height = 10)


	data_tmp[[lfcchangefrombaselinecol]]=log2(data_tmp[[targetcol]]/data_tmp[[basevalcol]])

	plotcol=glue("{targetcol}_lfcchangefrombaseline")
	plot = ggplot(data_tmp, aes(x = timepoint_day_clean, y = !!sym(lfcchangefrombaselinecol), color = subject_id)) +
		geom_point() +
		geom_line() +
		newdark_theme() +
		geom_smooth(colour = "black") +
		labs(title = glue("{plotcol} over time"), x = "Days relative to boost", y = glue("{plotcol}")) +
		theme(legend.position = "none",
			panel.grid.major = element_line(colour = "grey80"),
			panel.grid.minor = element_blank()) +
		scale_x_continuous(breaks = c(0, 1, 3, 7, 14, 30),limits=c(0,32))
	ggsave(glue("/user/leuven/342/vsc34263/myvibdata/challenge/plots/time/{targetcol}_lfcchangeover_time.png"), plot, width = width, height = 10)


	# data_sub=filter(data_tmp,timepoint_day_clean==14)
	# cor(data_sub[[basevalcol]],data_sub[[outcomecol]],use="pairwise.complete.obs",method="spearman")
	# cor(data_sub[[basevalcol]],data_sub[[changefrombaselinecol]],use="pairwise.complete.obs",method="spearman")
	# cor(data_sub[[basevalcol]],data_sub[[fcchangefrombaselinecol]],use="pairwise.complete.obs",method="spearman")
	# cor(data_sub[[basevalcol]],data_sub[[lfcchangefrombaselinecol]],use="pairwise.complete.obs",method="spearman")


	data_tmp$timepoint_day_4step=case_when(
		data_tmp$timepoint_day_clean<=0 ~ 0,
		data_tmp$timepoint_day_clean<=2 ~ 1,
		data_tmp$timepoint_day_clean<=6 ~ 3,
		TRUE ~ 7
	)
	data_agr=data_tmp %>% group_by(subject_id,timepoint_day_4step) %>% summarize_all(mean,na.rm=T) %>% ungroup()

	plot=ggplot(data_agr,aes(x=timepoint_day_4step,y=!!sym(targetcol),color=subject_id))+
		geom_point()+
		geom_line()+
		newdark_theme()+
		geom_smooth(colour="black")+
		labs(title=glue("{targetcol} over time"),x="Days relative to boost",y=glue("{targetcol}"))+
		theme(legend.position="none",
			panel.grid.major = element_line(colour = "grey80"),
			panel.grid.minor = element_blank())+
		scale_x_continuous(breaks=c(0,1,3,7))
	ggsave(glue("/user/leuven/342/vsc34263/myvibdata/challenge/plots/time/{targetcol}_over_time_agr.png"),plot,width=10,height=10)

	# predicting mean of post 7 days score is easier than predicting the actual score
	data_agr[[changecol]]=data_agr[[outcomecol]]-data_agr[[basevalcol]]
	data_agr[[fcchangefrombaselinecol]]=data_agr[[outcomecol]]/data_agr[[basevalcol]]
	data_agr[[lfcchangefrombaselinecol]]=log2(data_agr[[outcomecol]]/data_agr[[basevalcol]])
	data_agr_sub=filter(data_agr,timepoint_day_4step==7)

	cor(data_agr_sub[[basevalcol]],data_agr_sub[[outcomecol]],use="pairwise.complete.obs",method="spearman")
	cor(data_agr_sub[[basevalcol]],data_agr_sub[[changecol]],use="pairwise.complete.obs",method="spearman")
	cor(data_agr_sub[[basevalcol]],data_agr_sub[[fcchangefrombaselinecol]],use="pairwise.complete.obs",method="spearman")
	cor(data_agr_sub[[basevalcol]],data_agr_sub[[lfcchangefrombaselinecol]],use="pairwise.complete.obs",method="spearman")

	joindata = filter(data_tmp,timepoint_day_clean==14) %>% select(!!sym(targetcol),subject_id)
	colnames(joindata)=c(outcomecol,'subject_id')
	data_agr_sub=left_join(data_agr_sub,joindata,by="subject_id")
	joindata = filter(data,timepoint_day_clean==0) 

}