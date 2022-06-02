setwd("T:/vaccine/p704/analysis/manuscripts/NeutTiterBiomarker/code")

## ----- load package ---------
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)

# packageVersion("dplyr")
# [1] '1.0.4'

# packageVersion("ggplot2")
# [1] '3.3.3'

# packageVersion("tidyr")
# [1] '1.1.2'

# packageVersion("gridExtra")
# [1] '2.3'

# [1] "R version 4.0.2 (2020-06-22)"

##----- functions -------------
prettyNum0 <- function(x){sprintf("%.5g", x)}
source("triple-bnab-BH-titer.R")

## ----- read in data ---------
dat_ID_map_704 <- read.csv("../data/id_pubid_704.csv",stringsAsFactors = F)
dat_ID_map_703 <- read.csv("../data/id_pubid_703.csv",stringsAsFactors = F)
dat_ID_map <- bind_rows(dat_ID_map_704,dat_ID_map_703)

dat_704_pk <- read.csv("../data/nonmem_704.csv")
dat_703_pk <- read.csv("../data/nonmem_703.csv")
dat_monolix <- read.csv("../data/dat_amp_caseCtl_base.csv")

dat_703_placebo <- read.csv("../data/taenr_703.csv")
dat_704_placebo <- read.csv("../data/taenr_704.csv")
dat_placebo <- bind_rows(dat_703_placebo,dat_704_placebo) %>% 
  rename(PUB_ID=pub_id)

dat_ic50_ic80_703 <- read.csv("../data/v703_primary_cases_IC50_IC80.csv")
dat_ic50_ic80_704 <- read.csv("../data/v704_primary_cases_IC50_IC80.csv")
dat_ic50_ic80_raw <- bind_rows(dat_ic50_ic80_703,dat_ic50_ic80_704) 

dat_est_daily_all <- read.csv("../data/dat_est_daily_caseCtl.csv")
dat_est_daily_case <- readRDS("../data/dat_est_daily.rds")
conc_pred_medianWt <- readRDS("../data/conc_Ccpred_medianWt.rds")
dat_pred_bs <- readRDS("../data/conc_pred_daily_1000bs.rds")

dat_infTime_pfitter_days_703_sens1 <- read.csv("../data/v703_gpdx_tafd_20220202.csv")
dat_infTime_pfitter_days_704_sens1 <- read.csv("../data/v704_gpdx_tafd_20220202.csv")
dat_infTime_pfitter_days_sens1 <- bind_rows(dat_infTime_pfitter_days_703_sens1,dat_infTime_pfitter_days_704_sens1) %>% 
  select(pub_id,date,est_tafd,ci_low_95_tafd,ci_high_95_tafd,ci_low_99_tafd,ci_high_99_tafd,date_tafd,prob)

dat_eliCtl <- read.csv("../data/amp_eligible_ctrls.csv")

dat_ic_8mAb_704_raw <- read.csv("../data/VTN704_Non_Par_breakthrough_NAb_20210331.csv") 
dat_ic_8mAb_703_raw <- read.csv("../data/VTN703_Non_Par_breakthrough_NAb_20210628.csv") 

dat_sim_ss_16wkAR <- readRDS("../data/dat_sim_concSS_vrc07_pgt121_pgdm1400_forPE_16wkAR_1f2.5f5f.rds")
dat_sim_ss_24wkAR <- readRDS("../data/dat_sim_concSS_vrc07_pgt121_pgdm1400_forPE_24wkAR_1f2.5f5f.rds")
CcPred_vrc01_ss <- readRDS("../data/dat_sim_concSS_vrc01_d30.rds")
CcPred_vrc01_ss_d40 <- readRDS("../data/dat_sim_concSS_vrc01_d40.rds")

## --- Fig3: Distributions of predicted VRC01 serum ID80 titers (PT80s) 
##     against viruses acquired by placebo recipients, within each virus
##     neutralization IC80 sensitivity category ----------------

## get the trial information (dose, study)
dat_704_pk$study <- 704
dat_703_pk$study <- 703

dat_pk <- bind_rows(dat_704_pk,dat_703_pk) %>% 
  mutate(i=NULL) %>% 
  rename(dose=TRT)

dat_infu <- dat_pk %>% 
  filter(AMT!=0) %>% 
  mutate(CSCTRFLN=ifelse(CSCTRFLN==1,"Infected","Not Infected")
         ,CSCTRFLN=factor(CSCTRFLN,levels=c("Infected","Not Infected"))
         ,dose=paste0("VRC01 ",dose," mg/kg")
         ,study=ifelse(study==704,"HVTN704/HPTN085","HVTN703/HPTN081")
         ,study=factor(study,levels=c("HVTN704/HPTN085","HVTN703/HPTN081"))) 

dat_cov <- filter(dat_infu,TIME==0) %>% 
  select(ID,dose,CSCTRFLN,AGE:study) %>%
  mutate(sex=ifelse(sex==1,"male","female")
         ,race=ifelse(race==5,"White"
                      ,ifelse(race==3,"Black","Others"))
  )

## find partial control and exclude them
dat_ctl_wk80 <- filter(dat_pk,AMT==0) %>%
  filter(CSCTRFLN==0,AVISITN==2301) %>%
  select(ID,TIME)
dat_ctl_nowk80 <- filter(dat_pk,AMT==0) %>%
  filter(CSCTRFLN==0,AVISITN==2101,!ID%in%dat_ctl_wk80$ID) %>%
  select(ID,TIME) %>%
  mutate(TIME=TIME+56)

dat_ctl_wk80 <- bind_rows(dat_ctl_nowk80,dat_ctl_wk80)

## IC50/IC80 of AMP 
dat_ic50_ic80 <- dat_ic50_ic80_raw%>%
  select(protocol,rx_code,pub_id,gmt50ms,gmt80ms,gmt50ls,gmt80ls) %>%
  mutate(protocol=ifelse(protocol=="HVTN 703","HVTN703/HPTN081","HVTN704/HPTN085")
         ,rx_code=ifelse(rx_code=="T1","VRC01 10 mg/kg",ifelse(rx_code=="T2","VRC01 30 mg/kg","Control"))) %>%
  rename(PUB_ID=pub_id,dose=rx_code,study=protocol) %>%
  mutate(gmt50ms=as.numeric(ifelse(gmt50ms==">100","200",gmt50ms))
         ,gmt80ms=as.numeric(ifelse(gmt80ms==">100","200",gmt80ms))
         ,gmt50ls=as.numeric(ifelse(gmt50ls==">100","200",gmt50ls))
         ,gmt80ls=as.numeric(ifelse(gmt80ls==">100","200",gmt80ls)))

dat_ic50_ic80_placebo <- dat_ic50_ic80 %>%
  filter(dose=="Control"
         ,!is.na(gmt50ms)
  )

## estimated daily-grid concentrations for control only
dat_est_daily_ctl <- dat_est_daily_all %>%
  left_join(dat_ID_map) %>%
  left_join(select(dat_cov,ID,study,dose,CSCTRFLN)) %>%
  right_join(dat_ctl_wk80) %>%
  filter(time<=TIME) %>% 
  mutate(study="study-pooled")

dat_ic80_placebo_threeGrp <- dat_ic50_ic80_placebo %>%
  mutate(IC_grp=ifelse(gmt80ls<1,"IC80 < 1 \U00B5g/ml"
                       ,ifelse(gmt80ls<3&gmt80ls>=1,"IC80 1-3 \U00B5g/ml","IC80 \u2265 3 \U00B5g/ml"))
         ,IC_grp=factor(IC_grp,levels=c("IC80 \u2265 3 \U00B5g/ml","IC80 1-3 \U00B5g/ml","IC80 < 1 \U00B5g/ml"))
         ,study="study-pooled"
  ) %>%
  select(study,IC_grp,gmt80ls)

## calculate daily-grid ID80 by conc/IC80
dat_ID80_daily_ctl <- dat_est_daily_ctl %>%
  left_join(dat_ic80_placebo_threeGrp) %>%
  mutate(ID80=Cc/gmt80ls)

## add VE information
dat_ID80_daily_ctl_VE <- dat_ID80_daily_ctl %>%
  mutate(ve= ifelse(gmt80ls<1,75.4
                    ,ifelse(gmt80ls<3&gmt80ls>=1,4.2,3.3))
         ,ve_lc= ifelse(gmt80ls<1,45.5
                        ,ifelse(gmt80ls<3&gmt80ls>=1,-50,-48))
         ,ve_uc= ifelse(gmt80ls<1,88.9
                        ,ifelse(gmt80ls<3&gmt80ls>=1,56.0,36.8))
         ,IC_grp_short= ifelse(gmt80ls<1,"less 1"
                               ,ifelse(gmt80ls<3&gmt80ls>=1,"1 to 3","greater equal 3"))
         ,IC_grp_short=factor(IC_grp_short,levels=c("greater equal 3","1 to 3","less 1"))
  )

dat_ID80_daily_ctl_VE_sum <- dat_ID80_daily_ctl_VE %>%
  group_by(ve,ve_lc,ve_uc,IC_grp_short) %>%
  summarise(median_id80=median(ID80)
            ,lq_id80=quantile(ID80,0.25)
            ,uq_id80=quantile(ID80,0.75))

dat_ID80_daily_ctl_VE_sum_arrow <- filter(dat_ID80_daily_ctl_VE_sum,ve_lc==-50)


## start plot 
cairo_pdf("../figures/fig3_violin_boxplot_dailyID80_byVE.pdf",height=6,width = 8)

ggplot(dat_ID80_daily_ctl_VE,aes(x=ve,y=ID80,group=ve,color=IC_grp_short,fill=IC_grp_short))+
  geom_violin(width=30,size=0.8,alpha=0.2,fill="white")+
  geom_boxplot(outlier.color=NA,coef=0,width=5,alpha=0.3,data=filter(dat_ID80_daily_ctl_VE,gmt80ls>=3),show.legend = FALSE)+
  geom_boxplot(outlier.color=NA,coef=0,width=5,alpha=0.3,data=filter(dat_ID80_daily_ctl_VE,gmt80ls<3),show.legend = FALSE)+
  geom_point(aes(x=ve,y=median_id80),data=dat_ID80_daily_ctl_VE_sum,size=3,color="black")+
  geom_segment(aes(x=ve_lc,y=median_id80,xend=ve_uc,yend=median_id80)
               ,data=dat_ID80_daily_ctl_VE_sum
               ,color="black",size=0.6)+
  geom_segment(aes(x=ve_uc,y=median_id80,xend=ve_lc,yend=median_id80)
               ,data=dat_ID80_daily_ctl_VE_sum_arrow
               ,arrow = arrow(length = unit(0.2, "cm"),type="closed")
               ,color="black",size=0.6)+
  geom_text(data=filter(dat_ID80_daily_ctl_VE_sum,IC_grp_short=="greater equal 3")
            ,aes(y=median_id80,label=sprintf("%.1f", round(median_id80,1))),vjust=7,size=4.5,show.legend = FALSE, fontface = "bold",hjust=1)+
  geom_text(data=filter(dat_ID80_daily_ctl_VE_sum,IC_grp_short=="greater equal 3")
            ,aes(y=lq_id80,label=sprintf("%.1f", round(lq_id80,1))),vjust=7,size=4.5,show.legend = FALSE, fontface = "bold",hjust=1.1)+
  geom_text(data=filter(dat_ID80_daily_ctl_VE_sum,IC_grp_short=="greater equal 3")
            ,aes(y=uq_id80,label=sprintf("%.1f", round(uq_id80,1))),vjust=7,size=4.5,show.legend = FALSE, fontface = "bold",hjust=1.1)+
  geom_text(data=filter(dat_ID80_daily_ctl_VE_sum,IC_grp_short=="1 to 3")
            ,aes(y=median_id80,label=sprintf("%.1f", round(median_id80,1))),vjust=5,size=4.5,show.legend = FALSE, fontface = "bold",hjust=1)+
  geom_text(data=filter(dat_ID80_daily_ctl_VE_sum,IC_grp_short=="1 to 3")
            ,aes(y=lq_id80,label=sprintf("%.1f", round(lq_id80,1))),vjust=5,size=4.5,show.legend = FALSE, fontface = "bold",hjust=1.1)+
  geom_text(data=filter(dat_ID80_daily_ctl_VE_sum,IC_grp_short=="1 to 3")
            ,aes(y=uq_id80,label=sprintf("%.1f", round(uq_id80,1))),vjust=5,size=4.5,show.legend = FALSE, fontface = "bold",hjust=1)+
  geom_text(data=filter(dat_ID80_daily_ctl_VE_sum,IC_grp_short=="less 1")
            ,aes(y=median_id80,label=sprintf("%.1f", round(median_id80,1))),vjust=-4,size=4.5,show.legend = FALSE, fontface = "bold",hjust=0.5)+
  geom_text(data=filter(dat_ID80_daily_ctl_VE_sum,IC_grp_short=="less 1")
            ,aes(y=lq_id80,label=sprintf("%.1f", round(lq_id80,1))),vjust=-4,size=4.5,show.legend = FALSE, fontface = "bold",hjust=1.1)+
  geom_text(data=filter(dat_ID80_daily_ctl_VE_sum,IC_grp_short=="less 1")
            ,aes(y=uq_id80,label=sprintf("%.1f", round(uq_id80,1))),vjust=-4,size=4.5,show.legend = FALSE, fontface = "bold",hjust=0.5)+
  scale_y_log10(
    breaks = c(0.01,0.1,1,10,100,1000,10000)
    ,limits=c(0.01,10000)
    ,labels=prettyNum0
  )+
  scale_x_continuous(breaks=c(-50,-20,0,20,40,60,80,100),limits = c(-50,120))+
  scale_color_manual("",values=c("greater equal 3"="red"
                                 ,"1 to 3"="blue"
                                 ,"less 1"="forestgreen"),labels=c("IC80 \u2265 3 \U00B5g/ml"
                                                                   ,"IC80 1-3 \U00B5g/ml"
                                                                   ,"IC80 < 1 \U00B5g/ml"
                                 ))+
  scale_fill_manual("",values=c("greater equal 3"="red"
                                ,"1 to 3"="blue"
                                ,"less 1"="forestgreen")
                    ,labels=c("IC80 \u2265 3 \U00B5g/ml"
                              ,"IC80 1-3 \U00B5g/ml"
                              ,"IC80 < 1 \U00B5g/ml"
                    ))+
  xlab("Prevention Efficacy (%)")+
  ylab("Predicted VRC01 Serum ID80 Titers Against Placebo Virus Isolates")+
  coord_flip()+
  theme_bw(base_size = 15)+
  theme(legend.position = "bottom")

dev.off()

## ------ Extended Fig6. Distributions of predicted serum ID50 titers (PT50s) 
##        against viruses acquired by placebo recipients, within each virus neutralization
##        IC50 sensitivity category. ------------------------------------------------

## calculate daily-grid ID50 by conc/IC50
dat_ic50_placebo_threeGrp <- dat_ic50_ic80_placebo %>% 
  mutate(IC_grp=ifelse(gmt50ls<1,"IC50 < 1 \U00B5g/ml"
                       ,ifelse(gmt50ls<3&gmt50ls>=1,"IC50 1-3 \U00B5g/ml","IC50 \u2265 3 \U00B5g/ml"))
         ,IC_grp=factor(IC_grp,levels=c("IC50 \u2265 3 \U00B5g/ml","IC50 1-3 \U00B5g/ml","IC50 < 1 \U00B5g/ml"))
         ,study="study-pooled"
  ) %>% 
  select(study,IC_grp,gmt50ls)

dat_ID50_daily_ctl <- dat_est_daily_ctl %>%
  left_join(dat_ic50_placebo_threeGrp) %>% 
  mutate(ID50=Cc/gmt50ls) 

## add VE information
dat_ID50_daily_ctl_VE <- dat_ID50_daily_ctl %>%
  mutate(ve= ifelse(gmt50ls<1,39.0
                    ,ifelse(gmt50ls<3&gmt50ls>=1,23.9,0.8))
         ,ve_lc= ifelse(gmt50ls<1,-1.2
                        ,ifelse(gmt50ls<3&gmt50ls>=1,-47.3,-50))
         ,ve_uc= ifelse(gmt50ls<1,63.2
                        ,ifelse(gmt50ls<3&gmt50ls>=1,60.7,42.7))
         ,IC_grp_short= ifelse(gmt50ls<1,"less 1"
                               ,ifelse(gmt50ls<3&gmt50ls>=1,"1 to 3","greater equal 3"))
         ,IC_grp_short=factor(IC_grp_short,levels=c("greater equal 3","1 to 3","less 1"))
  )

dat_ID50_daily_ctl_VE_sum <- dat_ID50_daily_ctl_VE %>%
  group_by(ve,ve_lc,ve_uc,IC_grp_short) %>%
  summarise(median_id50=median(ID50)
            ,lq_id50=quantile(ID50,0.25)
            ,uq_id50=quantile(ID50,0.75))

dat_ID50_daily_ctl_VE_sum_arrow <- filter(dat_ID50_daily_ctl_VE_sum,ve_lc==-50)

## start plot
cairo_pdf("../figures/extFig6_violin_boxplot_dailyID50_byVE.pdf",height=6,width = 8)

ggplot(dat_ID50_daily_ctl_VE,aes(x=ve,y=ID50,group=ve,color=IC_grp_short,fill=IC_grp_short))+
  geom_violin(width=15,size=0.8,alpha=0.2,fill="white")+
  geom_boxplot(outlier.color=NA,coef=0,width=4,alpha=0.3,data=filter(dat_ID50_daily_ctl_VE,gmt50ls>=3),show.legend = FALSE)+
  geom_boxplot(outlier.color=NA,coef=0,width=4,alpha=0.3,data=filter(dat_ID50_daily_ctl_VE,gmt50ls<3),show.legend = FALSE)+
  geom_point(aes(x=ve,y=median_id50),data=dat_ID50_daily_ctl_VE_sum,size=3,color="black")+
  geom_segment(aes(x=ve_lc,y=median_id50,xend=ve_uc,yend=median_id50)
               ,data=dat_ID50_daily_ctl_VE_sum
               ,color="black")+
  geom_segment(aes(x=ve_uc,y=median_id50,xend=ve_lc,yend=median_id50)
               ,data=dat_ID50_daily_ctl_VE_sum_arrow
               ,arrow = arrow(length = unit(0.2, "cm"),type="closed")
               ,color="black")+
  geom_text(data=filter(dat_ID50_daily_ctl_VE_sum,IC_grp_short=="greater equal 3")
            ,aes(y=median_id50,label=sprintf("%.1f", round(median_id50,1))),vjust=3,size=4.5,show.legend = FALSE, fontface = "bold",hjust=1)+
  geom_text(data=filter(dat_ID50_daily_ctl_VE_sum,IC_grp_short=="greater equal 3")
            ,aes(y=lq_id50,label=sprintf("%.1f", round(lq_id50,1))),vjust=3,size=4.5,show.legend = FALSE, fontface = "bold",hjust=1)+
  geom_text(data=filter(dat_ID50_daily_ctl_VE_sum,IC_grp_short=="greater equal 3")
            ,aes(y=uq_id50,label=sprintf("%.1f", round(uq_id50,1))),vjust=3,size=4.5,show.legend = FALSE, fontface = "bold",hjust=1)+
  geom_text(data=filter(dat_ID50_daily_ctl_VE_sum,IC_grp_short=="1 to 3")
            ,aes(y=median_id50,label=sprintf("%.1f", round(median_id50,1))),vjust=3.5,size=4.5,show.legend = FALSE, fontface = "bold",hjust=1)+
  geom_text(data=filter(dat_ID50_daily_ctl_VE_sum,IC_grp_short=="1 to 3")
            ,aes(y=lq_id50,label=sprintf("%.1f", round(lq_id50,1))),vjust=3.5,size=4.5,show.legend = FALSE, fontface = "bold",hjust=1)+
  geom_text(data=filter(dat_ID50_daily_ctl_VE_sum,IC_grp_short=="1 to 3")
            ,aes(y=uq_id50,label=sprintf("%.1f", round(uq_id50,1))),vjust=3.5,size=4.5,show.legend = FALSE, fontface = "bold",hjust=0.5)+
  geom_text(data=filter(dat_ID50_daily_ctl_VE_sum,IC_grp_short=="less 1")
            ,aes(y=median_id50,label=sprintf("%.1f", round(median_id50,1))),vjust=-3,size=4.5,show.legend = FALSE, fontface = "bold",hjust=1)+
  geom_text(data=filter(dat_ID50_daily_ctl_VE_sum,IC_grp_short=="less 1")
            ,aes(y=lq_id50,label=sprintf("%.1f", round(lq_id50,1))),vjust=-3,size=4.5,show.legend = FALSE, fontface = "bold",hjust=1)+
  geom_text(data=filter(dat_ID50_daily_ctl_VE_sum,IC_grp_short=="less 1")
            ,aes(y=uq_id50,label=sprintf("%.1f", round(uq_id50,1))),vjust=-3,size=4.5,show.legend = FALSE, fontface = "bold",hjust=0.5)+
  scale_y_log10(
    breaks = c(0.01,0.1,1,10,100,1000,10000)
    ,limits=c(0.01,10000)
    ,labels=prettyNum0
  )+
  scale_x_continuous(breaks=c(-50,-20,0,20,40,60,80),limits = c(-50,80))+
  scale_color_manual("",values=c("greater equal 3"="red"
                                 ,"1 to 3"="blue"
                                 ,"less 1"="forestgreen"),labels=c("IC50 \u2265 3 \U00B5g/ml"
                                                                   ,"IC50 1-3 \U00B5g/ml"
                                                                   ,"IC50 < 1 \U00B5g/ml"
                                 ))+
  scale_fill_manual("",values=c("greater equal 3"="red"
                                ,"1 to 3"="blue"
                                ,"less 1"="forestgreen")
                    ,labels=c("IC50 \u2265 3 \U00B5g/ml"
                              ,"IC50 1-3 \U00B5g/ml"
                              ,"IC50 < 1 \U00B5g/ml"
                    ))+
  xlab("Prevention Efficacy (%)")+
  ylab("Predicted VRC01 Serum ID50 Titers Against Placebo Virus Isolates")+
  coord_flip()+
  theme_bw(base_size = 15)+
  theme(legend.position = "bottom")

dev.off()


## --------- Fig4 panel A: Predicted serum ID80 titers (PT80s) to autologous acquired viruses 
##           at HIV-1 acquistion among cases, and to placebo-recipient acquired viruses 
##           among non-cases -------------------------------------------------

# exclude 6 partical controls
# exclude 2 ppts who don't have RNA positive test: 703-1530 (dose 30), 703-2886 (dose 30)

# PDI excluded: 704-1654 (dose30), 704-1350 (dose 30), 704-0714 (dose 30)
# excluded 703-1357, not cases any more for primary analysis
# estimated infection time == 0, 704-2189 (dose 10)

## estimated infection time 
dat_infTime_sens1 <- dat_infTime_pfitter_days_sens1 %>%
  rename(TAFD_MEDIAN=est_tafd,TAFD_LB=ci_low_95_tafd,TAFD_UB=ci_high_95_tafd,TAFD_day=date_tafd,mass=prob,PUB_ID=pub_id) %>%
  mutate(PUB_ID=gsub("_","-",PUB_ID))

dat_infTime <- dat_infTime_sens1

## Geometric mean of IC50/IC80
dat_ic50_ic80_placebo_pool <- dat_ic50_ic80_placebo %>%
  mutate(study="Pooled AMP Trials") %>%
  bind_rows(dat_ic50_ic80_placebo)

dat_ctl_IC50_IC80_gm <- dat_ic50_ic80_placebo_pool%>%
  group_by(study) %>%
  summarise(gmt50ms=exp(mean(log(gmt50ms)))
            ,gmt80ms=exp(mean(log(gmt80ms)))
            ,gmt50ls=exp(mean(log(gmt50ls)))
            ,gmt80ls=exp(mean(log(gmt80ls))))

## predicted ID80 for Cases at estimated infection time 
dat_case_EstInfTime_cc_nab <- dat_est_daily_case %>%
  # bind_rows(data.frame(ID=unique(dat_est_daily_case$ID),time=0,Cc=0.01)) %>%
  left_join(dat_ID_map) %>%
  left_join(distinct(dat_infTime,PUB_ID,TAFD_MEDIAN)) %>%
  mutate(TAFD_MEDIAN=ifelse(TAFD_MEDIAN<0,0,TAFD_MEDIAN)) %>%
  filter(time==TAFD_MEDIAN) %>%
  left_join(select(dat_cov,ID,study,dose)) %>%
  mutate(case_ctl="Infected") %>%
  left_join(select(dat_ic50_ic80,PUB_ID,gmt50ms,gmt80ms,gmt50ls,gmt80ls)) %>%
  mutate(Cc=ifelse(Cc<=0.01,0.01,Cc)
         ,ID50=Cc/gmt50ms
         ,ID50_ls=Cc/gmt50ls
         ,ID80=Cc/gmt80ms
         ,ID80_ls=Cc/gmt80ls
  ) %>%
  filter(!PUB_ID%in%c("703-1357","704-1654","704-0714","704-1350"))

## predicted ID50/ID80 for control = median CC/GM of IC50/IC80 
dat_est_median_ctl <- dat_est_daily_all %>%
  left_join(dat_ID_map) %>%
  left_join(select(dat_cov,ID,study,dose,CSCTRFLN)) %>%
  right_join(dat_ctl_wk80) %>%
  filter(time<=TIME) %>%
  group_by(study,dose,PUB_ID) %>%
  summarise(Cc_median=median(Cc))

dat_est_median_ctl_cc_nab <- dat_est_median_ctl %>%
  ungroup() %>%
  left_join(dat_ctl_IC50_IC80_gm) %>%
  mutate(ID50=Cc_median/gmt50ms
         ,ID80=Cc_median/gmt80ms
         ,ID50_ls=Cc_median/gmt50ls
         ,ID80_ls=Cc_median/gmt80ls
         ,case_ctl="Not Infected") %>%
  rename(Cc=Cc_median) %>%
  left_join(dat_ID_map)

dat_est_nab_cc <- bind_rows(dat_case_EstInfTime_cc_nab
                            ,dat_est_median_ctl_cc_nab)
dat_est_nab_cc_pool <- dat_est_nab_cc %>%
  mutate(study="Pooled AMP Trials") %>%
  bind_rows(dat_est_nab_cc)

## prepare data for ID80 plot
dat <- dat_est_nab_cc_pool %>%
  filter(study=="Pooled AMP Trials",case_ctl=="Infected") %>%
  select(PUB_ID,ID80_ls,case_ctl,dose) %>%
  rename(ID80=ID80_ls,CSCTRFLN=case_ctl) %>%
  bind_rows(select(dat_ID80_daily_ctl,PUB_ID,ID80,CSCTRFLN,dose)) %>%
  mutate(study="Pooled AMP Trials")

dat_sum <- dat %>%
  filter(!is.na(ID80)) %>% 
  mutate(resp=ifelse(ID80==1,0,1)) %>% 
  group_by(study,dose,CSCTRFLN) %>%
  summarise(median_id80=median(ID80,na.rm=T)
            ,lq_id80=quantile(ID80,0.25,na.rm=T)
            ,uq_id80=quantile(ID80,0.75,na.rm=T)
            ,resp_rate=paste0("Resp.% = ",round(sum(resp)/n()*100,1),"%")) %>%
  mutate(x.pos=ifelse(CSCTRFLN=="Infected",0.6,1.5))

## start plot
cairo_pdf("../figures/fig4_panelA_violinPlot_estID80Ls_ctlDaily_percentile.pdf",height=5,width = 6)

ggplot(dat,aes(x=CSCTRFLN,y=ID80,color=CSCTRFLN))+
  geom_violin(width=0.5,show.legend = FALSE)+
  geom_boxplot(outlier.shape = NA, coef = 0,width=0.15
               ,aes(fill=CSCTRFLN),alpha=0.2
  )+
  geom_jitter(position=position_jitter(width=0.1,height=0),aes(color=CSCTRFLN)
              ,data=filter(dat,CSCTRFLN=="Infected")
              ,size=1)+
  geom_text(aes(x=x.pos,y=median_id80,label=sprintf("%.1f",round(median_id80,1))),data=dat_sum,fontface="bold",size=3.5,show.legend = FALSE)+
  geom_text(aes(x=x.pos,y=lq_id80,label=sprintf("%.1f",round(lq_id80,1))),data=dat_sum,fontface="bold",size=3.5,show.legend = FALSE)+
  geom_text(aes(x=x.pos,y=uq_id80,label=sprintf("%.1f",round(uq_id80,1))),data=dat_sum,fontface="bold",size=3.5,show.legend = FALSE)+
  scale_y_log10(breaks=c(0.001,0.01,0.1,1,10,100,1000,10000),labels=prettyNum0)+
  scale_color_manual("",values=c("Infected"="red","Not Infected"="blue"))+
  scale_fill_manual("",values=c("Infected"="red","Not Infected"="blue"))+
  facet_grid(~dose)+
  xlab("")+
  ylab("Predicted VRC01 Serum ID80 Titers Against Autologous\n(Infected) or Placebo (Not Infected) Virus Isolates")+
  ggtitle("A")+
  theme_bw()+
  theme(legend.position="bottom"
        ,axis.title = element_text(size=11)
        ,axis.text=element_text(size=12)
        ,legend.text = element_text(size=10)
        ,strip.text = element_text(size=13)
        ,plot.title = element_text(size = 15, face = "bold")
        ,plot.title.position = "plot"
        ,plot.margin = unit(c(0.2,1,0.1,1),"cm")
  )

dev.off()


## --------- ExtFig7 panel A: Predicted serum ID50 titers (PT50s) to autologous acquired viruses 
##           at HIV-1 acquistion among cases, and to placebo-recipient acquired viruses 
##           among non-cases -------------------------------------------------

dat <- dat_est_nab_cc_pool %>%
  filter(study=="Pooled AMP Trials",case_ctl=="Infected") %>%
  select(PUB_ID,ID50_ls,case_ctl,dose) %>%
  rename(ID50=ID50_ls,CSCTRFLN=case_ctl) %>%
  bind_rows(select(dat_ID50_daily_ctl,PUB_ID,ID50,CSCTRFLN,dose)) %>%
  mutate(study="Pooled AMP Trials"
  )

dat_sum <- dat %>%
  filter(!is.na(ID50)) %>% 
  mutate(resp=ifelse(ID50==1,0,1)) %>% 
  group_by(study,dose,CSCTRFLN) %>%
  summarise(median_id50=median(ID50,na.rm=T)
            ,lq_id50=quantile(ID50,0.25,na.rm=T)
            ,uq_id50=quantile(ID50,0.75,na.rm=T)
            ,resp_rate=paste0("Resp.% = ",sprintf("%.1f",round(sum(resp)/n()*100,1)),"%"))%>%
  mutate(x.pos=ifelse(CSCTRFLN=="Infected",0.6,1.5))

cairo_pdf("../figures/extFig7_panelA_violinPlot_estID50Ls_ctlDaily_percentile.pdf",height=5,width = 6)

P_A_id50 <- ggplot(dat,aes(x=CSCTRFLN,y=ID50,color=CSCTRFLN))+
  geom_violin(width=0.5,show.legend = FALSE)+
  geom_boxplot(outlier.shape = NA, coef = 0,width=0.1
               ,aes(fill=CSCTRFLN),alpha=0.2
  )+
  geom_jitter(position=position_jitter(width=0.1,height=0),aes(color=CSCTRFLN)
              ,data=filter(dat,CSCTRFLN=="Infected")
              ,size=1)+
  geom_text(aes(x=x.pos,y=median_id50,label=sprintf("%.1f",round(median_id50,1))),data=dat_sum,fontface="bold",size=3.5,show.legend = FALSE)+
  geom_text(aes(x=x.pos,y=lq_id50,label=sprintf("%.1f",round(lq_id50,1))),data=dat_sum,fontface="bold",size=3.5,show.legend = FALSE)+
  geom_text(aes(x=x.pos,y=uq_id50,label=sprintf("%.1f",round(uq_id50,1))),data=dat_sum,fontface="bold",size=3.5,show.legend = FALSE)+
  scale_y_log10(breaks=c(0.001,0.01,0.1,1,10,100,1000,10000),labels=prettyNum0)+
  scale_color_manual("",values=c("Infected"="red","Not Infected"="blue"))+
  scale_fill_manual("",values=c("Infected"="red","Not Infected"="blue"))+
  facet_grid(~dose)+
  xlab("")+
  ylab("Predicted VRC01 Serum ID50 Titers Against Autologous\n(Infected) or Placebo (Not Infected) Virus Isolates")+
  ggtitle("A")+
  theme_bw()+
  theme(legend.position="bottom"
        ,axis.title = element_text(size=11)
        ,axis.text=element_text(size=12)
        ,legend.text = element_text(size=10)
        ,strip.text = element_text(size=13)
        ,plot.title = element_text(size = 15, face = "bold")
        ,plot.title.position = "plot"
        ,plot.margin = unit(c(0.2,1,0.1,1),"cm")
  )

print(P_A_id50)

dev.off()


## --------- Fig4 panel B: Predicted serum ID80 titers (PT80s) to autologous acquired viruses 
##           at HIV-1 acquistion among cases, and to placebo-recipient acquired viruses 
##           among non-cases -------------------------------------------------

## weights for controls which comes from all eligible control data: 
## number in each stratum of eligible ctl/total eligible ctl 

dat_eliCtl_noPDI <- filter(dat_eliCtl,is.na(pdistatus),delta2==1)

n_total <- nrow(dat_eliCtl_noPDI)

dat_wt_dPool_sPool <- dat_eliCtl_noPDI %>% 
  group_by(Protocol,rx_code) %>% 
  summarise(n=n()) %>% 
  mutate(f=n/n_total
         ,rx_code=ifelse(rx_code=="T1","VRC01 10 mg/kg","VRC01 30 mg/kg")
         ,Protocol=ifelse(Protocol=="HVTN 703","HVTN703/HPTN081","HVTN704/HPTN085")
         ,n=NULL) %>% 
  rename(study=Protocol,dose=rx_code)

dat_wt_dPool <- dat_eliCtl_noPDI %>% 
  group_by(Protocol) %>% 
  mutate(n_total_dPool=n()) %>% 
  group_by(Protocol,rx_code,n_total_dPool) %>% 
  summarise(n=n()) %>% 
  mutate(f=n/n_total_dPool
         ,rx_code=ifelse(rx_code=="T1","VRC01 10 mg/kg","VRC01 30 mg/kg")
         ,Protocol=ifelse(Protocol=="HVTN 703","HVTN703/HPTN081","HVTN704/HPTN085")
         ,n=NULL
         ,n_total_dPool=NULL) %>% 
  rename(study=Protocol,dose=rx_code)

dat_wt_sPool <- dat_eliCtl_noPDI %>% 
  group_by(rx_code) %>% 
  mutate(n_total_sPool=n()) %>% 
  group_by(Protocol,rx_code,n_total_sPool) %>% 
  summarise(n=n()) %>% 
  mutate(f=n/n_total_sPool
         ,rx_code=ifelse(rx_code=="T1","VRC01 10 mg/kg","VRC01 30 mg/kg")
         ,Protocol=ifelse(Protocol=="HVTN 703","HVTN703/HPTN081","HVTN704/HPTN085")
         ,n=NULL
         ,n_total_sPool=NULL
  ) %>% 
  rename(study=Protocol,dose=rx_code)

## function for marker method weight adjustment 
f_markerMethod_pool <- function(dat_sum0,dat_w_case_dpool,dat_w_case_spool,dat_w_case_dspool){
  #---- case -------
  dat_sum_case0 <- dat_sum0 %>%
    filter(case_ctl=="Infected")
  
  dat_sum_dpool_case <- dat_sum_case0 %>% 
    left_join(dat_w_case_dpool) %>% 
    mutate(mu_wt=mu*f
           ,var_wt=var*f^2) %>% 
    group_by(study,case_ctl) %>% 
    summarise(mu=sum(mu_wt)
              ,var=sum(var_wt)
              ,dose="Pooled VRC01")
  
  dat_sum_spool_case <- dat_sum_case0 %>% 
    left_join(dat_w_case_spool) %>% 
    mutate(mu_wt=mu*f
           ,var_wt=var*f^2) %>% 
    group_by(dose,case_ctl) %>% 
    summarise(mu=sum(mu_wt)
              ,var=sum(var_wt)
              ,study="Pooled AMP Trials")
  
  dat_sum_dspool_case <- dat_sum_case0 %>% 
    left_join(dat_w_case_dspool) %>% 
    mutate(mu_wt=mu*f
           ,var_wt=var*f^2) %>% 
    group_by(case_ctl) %>% 
    summarise(mu=sum(mu_wt)
              ,var=sum(var_wt)
              ,study="Pooled AMP Trials"
              ,dose="Pooled VRC01")
  
  dat_sum_case <- bind_rows(dat_sum_case0,dat_sum_dpool_case,dat_sum_spool_case,dat_sum_dspool_case)
  
  #--- ctl -------
  dat_sum_ctl0 <- dat_sum0 %>%
    filter(case_ctl=="Not Infected")
  
  dat_sum_dpool_ctl <- dat_sum_ctl0 %>% 
    left_join(dat_wt_dPool) %>% 
    mutate(mu_wt=mu*f
           ,var_wt=var*f^2) %>% 
    group_by(study,case_ctl) %>% 
    summarise(mu=sum(mu_wt)
              ,var=sum(var_wt)) %>% 
    mutate(dose="Pooled VRC01")
  
  dat_sum_stypool_ctl <- dat_sum_ctl0 %>% 
    left_join(dat_wt_sPool) %>% 
    mutate(mu_wt=mu*f
           ,var_wt=var*f^2) %>% 
    group_by(dose,case_ctl) %>% 
    summarise(mu=sum(mu_wt)
              ,var=sum(var_wt)) %>% 
    mutate(study="Pooled AMP Trials")
  
  dat_sum_allpool_ctl <- dat_sum_ctl0 %>% 
    left_join(dat_wt_dPool_sPool)%>% 
    mutate(mu_wt=mu*f
           ,var_wt=var*f^2) %>% 
    group_by(case_ctl) %>% 
    summarise(mu=sum(mu_wt)
              ,var=sum(var_wt)) %>% 
    mutate(study="Pooled AMP Trials",dose="Pooled VRC01")
  
  
  dat_sum_ctl <- bind_rows(dat_sum_dpool_ctl,dat_sum_stypool_ctl,dat_sum_allpool_ctl,dat_sum_ctl0) 
  
  #----- combine case and control -----------
  dat_sum <- bind_rows(dat_sum_case,dat_sum_ctl)
  
  dat_sum_diff <- dat_sum %>%
    group_by(study,dose) %>%
    summarise(mu=mu[case_ctl=="Not Infected"]-mu[case_ctl=="Infected"]
              ,var=sum(var)
              ,case_ctl="ctl-case")
  
  dat_sum_ci <- bind_rows(dat_sum,dat_sum_diff) %>%
    ungroup() %>%
    mutate(se=sqrt(var)
           ,lc=mu-1.96*se
           ,uc=mu+1.96*se
           ,z = mu/se
           ,p = 2*(1- pnorm(abs(z)))
           ,type=paste0(dose,".",case_ctl)
           ,type=factor(type,levels=c("Pooled VRC01.ctl-case"
                                      ,"VRC01 30 mg/kg.ctl-case"
                                      ,"VRC01 10 mg/kg.ctl-case"
                                      ,"Pooled VRC01.Infected"
                                      ,"Pooled VRC01.Not Infected"
                                      ,"VRC01 30 mg/kg.Infected"
                                      ,"VRC01 30 mg/kg.Not Infected"
                                      ,"VRC01 10 mg/kg.Infected"
                                      ,"VRC01 10 mg/kg.Not Infected"
           ))
           ,case_ctl=ifelse(case_ctl=="ctl-case","(Not Infected)/Infected",case_ctl)
           ,case_ctl = factor(case_ctl,levels=c("Infected","Not Infected","(Not Infected)/Infected"))
           ,text = paste0(round(exp(mu),1)," (",round(exp(lc),1),", ",round(exp(uc),1),")")
    )
  return(dat_sum_ci)
}

## FP ID80 function
f_fp_id80 <- function(df){
  
  ggplot(df,aes(y=type,x=exp(mu),color=case_ctl))+
    geom_vline(aes(xintercept =1),linetype="dashed",color="forestgreen",size=0.8)+
    geom_point()+
    geom_errorbarh(aes(xmin=exp(lc),xmax=exp(uc)),height=0.3)+
    geom_text(aes(y=type,x=exp(mu),label=text),vjust=-0.8,size=3,color="black")+
    xlab("Geometric Mean of ID80 and Ratios (95% CI)")+
    ylab("")+
    ggtitle(unique(df$study))+
    scale_y_discrete(labels = c('VRC01 30 mg/kg.ctl-case' = "VRC01 30 mg/kg, (Not Infected)/Infected"
                                ,'VRC01 10 mg/kg.ctl-case' = "VRC01 10 mg/kg, (Not Infected)/Infected"
                                ,'Pooled VRC01.ctl-case' = "VRC01 Pooled, (Not Infected)/Infected"
                                ,'Pooled VRC01.Not Infected' = "VRC01 Pooled, Not Infected"
                                ,'Pooled VRC01.Infected' = "VRC01 Pooled, Infected"
                                ,'VRC01 30 mg/kg.Not Infected' = "VRC01 30 mg/kg, Not Infected"
                                ,'VRC01 30 mg/kg.Infected' = "VRC01 30 mg/kg, Infected"
                                ,'VRC01 10 mg/kg.Not Infected' = "VRC01 10 mg/kg, Not Infected"
                                ,'VRC01 10 mg/kg.Infected' = "VRC01 10 mg/kg, Infected"
    )
    ,expand = expand_scale(add = 0.8)
    
    )+
    scale_color_manual("",values=c("Infected"="red","Not Infected"="blue","(Not Infected)/Infected"="black"))+
    scale_x_log10(limits=c(exp(min(dat_sum_ci_id80$lc)),exp(max(dat_sum_ci_id80$uc)))
                  ,expand = expand_scale(add = 0.2))+
    theme_bw()+
    theme(legend.position="bottom"
          ,axis.title = element_text(size=11)
          ,axis.text.y=element_text(size=12,hjust=0)
          ,axis.text.x=element_text(size=12)
          ,legend.text = element_text(size=10)
          ,strip.text = element_text(size=13)
          ,plot.title = element_text(size = 15, face = "bold")
          ,plot.title.position = "plot"
          ,plot.margin = unit(c(0.2,1,0.2,0.2),"cm")
    )
}

## function of marker method from yanqing's method for Nab data 

f_markerMethod_yanqing_nab <- function(var_IC,var_nab){
  
  var_IC <- rlang::sym(var_IC)
  var_nab <- rlang::sym(var_nab)
  
  # mu from bootstrap 
  
  dat_mu_bs_sampleData_list <- lapply(1:1000, function(i){
    # browser()
    file <- paste0(i,".rds")
    dat <- readRDS(file.path("../data/Conc_pred_dailyGrid_sampleData/",file))
    
    #case
    dat_case <- filter(dat,CSCTRFLN==1) %>%
      inner_join(dat_infTime_g0_newID) %>%
      inner_join(dat_ic50_ic80_case_newID) %>%
      mutate(nab=Cc/!!var_IC
             ,log_nab=log(nab)) %>%
      group_by(ID,study,dose,CSCTRFLN,old_ID,b) %>%
      summarise(integ_logNab=sum(mass*log_nab)/sum(mass))
    
    dat_mu_case <- dat_case %>%
      group_by(study,dose,CSCTRFLN,b) %>%
      summarise(mu=mean(integ_logNab))
    
    #control
    dat_ctl <- filter(dat,CSCTRFLN==0) %>%
      inner_join(dat_ctl_wk80_newID) %>%
      filter(time<=TIME) %>%
      group_by(ID,study,dose,CSCTRFLN,old_ID,b) %>%
      summarise(median_Cc=median(Cc))
    
    dat_mu_ctl <- dat_ctl %>%
      left_join(dat_ctl_IC50_IC80_gm_studyN) %>%
      mutate(nab=median_Cc/!!var_IC
             ,log_nab=log(nab)
      ) %>%
      group_by(study,dose,CSCTRFLN,b) %>%
      summarise(mu=mean(log_nab) )
    
    dat_mu <- bind_rows(dat_mu_case,dat_mu_ctl)
    return(dat_mu)
    
  })
  
  
  dat_mu_bs_sampleData <- bind_rows(dat_mu_bs_sampleData_list) %>%
    ungroup() %>%
    mutate(study=ifelse(study==1,"HVTN704/HPTN085","HVTN703/HPTN081")
           ,dose=ifelse(dose==1,"VRC01 30 mg/kg","VRC01 10 mg/kg")
           ,CSCTRFLN=ifelse(CSCTRFLN==1,"Infected","Not Infected")) %>%
    rename(mu_bs=mu,case_ctl=CSCTRFLN)
  
  
  # raw data mu case and ctl
  dat_integ_case <- dat_est_daily_case %>%
    inner_join(dat_infTime_g0) %>%
    inner_join(dat_ic50_ic80_case) %>%
    left_join(select(dat_cov,ID,study,dose,CSCTRFLN)) %>%
    mutate(nab=Cc/!!var_IC
           ,log_nab=log(nab)) %>%
    group_by(ID,study,dose,CSCTRFLN) %>%
    summarise(integ_logNab=sum(mass*log_nab)/sum(mass))
  
  dat_mu_case <- dat_integ_case %>%
    group_by(study,dose,CSCTRFLN) %>%
    summarise(mu=mean(integ_logNab)) %>% 
    rename(case_ctl=CSCTRFLN)
  
  dat_mu_ctl <- dat_est_median_ctl_cc_nab %>%
    group_by(study,dose,case_ctl) %>%
    summarise(mu=mean(log(!!var_nab)) ) 
  
  dat_mu_raw <- bind_rows(dat_mu_case,dat_mu_ctl)
  
  dat_sum_yanqing <- left_join(dat_mu_raw,dat_mu_bs_sampleData) %>%
    group_by(study,dose,case_ctl,mu) %>%
    summarise(var=mean((mu_bs-mu)^2) )
  
  dat_w_case_dpool_yanqing <- dat_integ_case %>% 
    rename(case_ctl=CSCTRFLN) %>% 
    filter(!is.na(integ_logNab),case_ctl=="Infected") %>% 
    group_by(study,dose,case_ctl) %>% 
    summarise(n=n()) %>% 
    group_by(study,case_ctl) %>% 
    mutate(n_dpool=sum(n)
           ,f=n/n_dpool
           ,n=NULL
           ,n_dpool=NULL)
  
  dat_w_case_spool_yanqing <- dat_integ_case %>% 
    rename(case_ctl=CSCTRFLN) %>% 
    filter(!is.na(integ_logNab),case_ctl=="Infected") %>% 
    group_by(study,dose,case_ctl) %>% 
    summarise(n=n()) %>% 
    group_by(dose,case_ctl) %>% 
    mutate(n_spool=sum(n)
           ,f=n/n_spool
           ,n=NULL
           ,n_spool=NULL)
  
  dat_w_case_dspool_yanqing <- dat_integ_case %>% 
    rename(case_ctl=CSCTRFLN) %>% 
    filter(!is.na(integ_logNab),case_ctl=="Infected") %>% 
    group_by(study,dose,case_ctl) %>% 
    summarise(n=n()) %>% 
    group_by(case_ctl) %>% 
    mutate(n_dspool=sum(n)
           ,f=n/n_dspool
           ,n=NULL
           ,n_dspool=NULL)
  
  dat_sum_ci_bs <- f_markerMethod_pool(dat_sum_yanqing
                                       ,dat_w_case_dpool_yanqing
                                       ,dat_w_case_spool_yanqing
                                       ,dat_w_case_dspool_yanqing) %>% 
    mutate(outcome="mu")
  
  
  return(dat_sum_ci_bs)
}

## ID50/ID80 data prepare
dat_infTime_g0 <- dat_infTime %>%
  inner_join(dat_ID_map) %>%
  filter(mass>0,TAFD_day>0,!PUB_ID%in%c("703-1357","704-1654","704-1350","704-0714")) %>%
  select(ID,TAFD_day,mass) %>%
  rename(time=TAFD_day)

dat_infTime_g0_newID <- dat_infTime_g0 %>%
  rename(old_ID=ID)

dat_ctl_wk80_newID <- dat_ctl_wk80 %>%
  rename(old_ID=ID)

dat_ic50_ic80_case <- dat_ic50_ic80 %>%
  filter(dose!="Control",!is.na(gmt50ms)) %>%
  left_join(dat_ID_map) %>%
  select(ID,gmt50ms,gmt80ms,gmt50ls,gmt80ls)

dat_ic50_ic80_case_newID <- dat_ic50_ic80_case %>%
  rename(old_ID=ID)

dat_ctl_IC50_IC80_gm_studyN <- dat_ctl_IC50_IC80_gm %>%
  mutate(study=ifelse(study=="HVTN704/HPTN085",1,ifelse(study=="HVTN703/HPTN081",0,2))) 

## marker method calculation
dat_sum_ci_id80 <- f_markerMethod_yanqing_nab(var_IC = "gmt80ls",var_nab = "ID80_ls")

## start plot
cairo_pdf("../figures/fig4_panelB_fp_ID80Ls_yanqing_Pfitter.pdf",height=6,width = 8,onefile = T)
P_list <- plyr::dlply(dat_sum_ci_id80,plyr::.(study),f_fp_id80)
P_list[[3]]
dev.off()

## --------- ExtFig7 panel A: Predicted serum ID50 titers (PT50s) to autologous acquired viruses 
##           at HIV-1 acquistion among cases, and to placebo-recipient acquired viruses 
##           among non-cases -------------------------------------------------


## FP ID50 function
f_fp_id50 <- function(df){
  
  ggplot(df,aes(y=type,x=exp(mu),color=case_ctl))+
    geom_vline(aes(xintercept =1),linetype="dashed",color="forestgreen",size=0.8)+
    geom_point()+
    geom_errorbarh(aes(xmin=exp(lc),xmax=exp(uc)),height=0.3)+
    geom_text(aes(y=type,x=exp(mu),label=text),vjust=-0.8,size=3,color="black")+
    xlab("Geometric Mean of ID50 and Ratios (95% CI)")+
    ylab("")+
    ggtitle(unique(df$study))+
    scale_y_discrete(labels = c('VRC01 30 mg/kg.ctl-case' = "VRC01 30 mg/kg, (Not Infected)/Infected"
                                ,'VRC01 10 mg/kg.ctl-case' = "VRC01 10 mg/kg, (Not Infected)/Infected"
                                ,'Pooled VRC01.ctl-case' = "VRC01 Pooled, (Not Infected)/Infected"
                                ,'Pooled VRC01.Not Infected' = "VRC01 Pooled, Not Infected"
                                ,'Pooled VRC01.Infected' = "VRC01 Pooled, Infected"
                                ,'VRC01 30 mg/kg.Not Infected' = "VRC01 30 mg/kg, Not Infected"
                                ,'VRC01 30 mg/kg.Infected' = "VRC01 30 mg/kg, Infected"
                                ,'VRC01 10 mg/kg.Not Infected' = "VRC01 10 mg/kg, Not Infected"
                                ,'VRC01 10 mg/kg.Infected' = "VRC01 10 mg/kg, Infected"
    )
    ,expand = expand_scale(add = 0.8)
    
    )+
    scale_color_manual("",values=c("Infected"="red","Not Infected"="blue","(Not Infected)/Infected"="black"))+
    scale_x_log10(limits=c(exp(min(dat_sum_ci_id50$lc)),exp(max(dat_sum_ci_id50$uc)))
                  ,expand = expand_scale(add = 0.2))+
    theme_bw()+
    theme(legend.position="bottom"
          ,axis.title = element_text(size=11)
          ,axis.text.y=element_text(size=12,hjust=0)
          ,axis.text.x=element_text(size=12)
          ,legend.text = element_text(size=10)
          ,strip.text = element_text(size=13)
          ,plot.title = element_text(size = 15, face = "bold")
          ,plot.title.position = "plot"
          ,plot.margin = unit(c(0.2,1,0.2,0.2),"cm")
    )
}

## marker method calculation
dat_sum_ci_id50 <- f_markerMethod_yanqing_nab(var_IC = "gmt50ls",var_nab = "ID50_ls")

## start plot
cairo_pdf("../figures/extFig7_panelB_fp_ID50Ls_yanqing_Pfitter.pdf",height=6,width = 8,onefile = T)
P_list <- plyr::dlply(dat_sum_ci_id50,plyr::.(study),f_fp_id50)
P_list[[3]]
dev.off()

## -------- Extended fig5: Estimated serum VRC01 concentration at the estimated time 
##          of infection (median of the Bayesian posterior distribution of infection time) 
##          since last infusion in all primary endpoint HIV-1 cases, by randomization arm 
##          in (A) HVTN 704/HPTN 085 and (B) HVTN 703/HPTN 081. --------------------

dat_infTime_great0 <- dat_infTime %>%
  filter(mass>0,TAFD_day>0) %>%
  select(PUB_ID,TAFD_day,mass) %>%
  filter(PUB_ID%in%dat_case_EstInfTime_cc_nab$PUB_ID) %>%
  left_join(dat_ID_map)

set.seed(123)
dat_pred_bsPost_90CI <- plyr::ddply(dat_infTime_great0,plyr::.(ID),function(df){
  cur_id <- unique(df$ID)
  time_sample <- sample(df$TAFD_day,size=1000,prob = df$mass,replace = T)
  cur_dat_pred_bs <- filter(dat_pred_bs,ID==cur_id)
  dat_pred_1000000 <- data.frame(time=rep(time_sample,1000),b=rep(1:1000,each=length(time_sample))) %>%
    left_join(cur_dat_pred_bs)
  dat_pred_90CI <- dat_pred_1000000 %>%
    group_by(ID) %>%
    summarise(Cc_ub=quantile(Cc,0.95)
              ,Cc_lb=quantile(Cc,0.05))
  return(dat_pred_90CI)
}
)

dat_monolix_noConc <- filter(dat_monolix,TIME==-1)
dat_cases_noPk <- dat_monolix_noConc %>%
  left_join(dat_ID_map)

dat_estConcAtInfTime_sinceLastInfu <- dat_case_EstInfTime_cc_nab %>%
  left_join(select(dat_infu,ID,TIME,AMT)) %>%
  mutate(infTimeSinceLastInfu=time-TIME) %>%
  filter(infTimeSinceLastInfu>0) %>%
  group_by(PUB_ID,ID,Cc) %>%
  summarise(infTimeSinceLastInfu=min(infTimeSinceLastInfu)) %>%
  left_join(select(dat_cov,ID,study,dose))

dat_estConcAtInfTime_sinceLastInfu_placebo <- dat_infTime %>%
  distinct(PUB_ID,TAFD_MEDIAN) %>%
  right_join(dat_placebo) %>%
  mutate(infTimeSinceLastInfu=TAFD_MEDIAN-TAENR) %>%
  filter(infTimeSinceLastInfu>0) %>%
  group_by(PUB_ID) %>%
  summarise(infTimeSinceLastInfu=min(infTimeSinceLastInfu)) %>%
  mutate(Cc=0.01
         ,study=ifelse(grepl("703-",PUB_ID),"HVTN703/HPTN081","HVTN704/HPTN085")
         ,dose="Control")

dat_estConcAtInfTime_sinceLastInfu <- bind_rows(dat_estConcAtInfTime_sinceLastInfu,dat_estConcAtInfTime_sinceLastInfu_placebo)

dat_pred_bsPost_90CI_time <- dat_pred_bsPost_90CI %>%
  left_join(dat_estConcAtInfTime_sinceLastInfu) %>%
  filter(!PUB_ID%in%dat_cases_noPk$PUB_ID) %>%
  filter(PUB_ID!="704-2277")

## prediction interval plot function
f_predInterval <- function(df){
  title <- unique(df$study)
  cur_dat_pred_bsPost_90CI_time <- filter(dat_pred_bsPost_90CI_time,study==title)
  cur_dat_estConcAtInfTime_sinceLastInfu <- filter(dat_estConcAtInfTime_sinceLastInfu,study==title)
  
  ggplot()+
    geom_line(data=df,aes(y=conc_median,x=time),size=0.3)+
    geom_ribbon(aes(ymax=conc_ub,ymin=conc_lb,x=time),color="lightgrey",alpha=0.1,data=df)+
    geom_errorbar(aes(ymin=Cc_lb, ymax=Cc_ub,x=infTimeSinceLastInfu)
    , width=3
    ,size=0.5
    ,data=cur_dat_pred_bsPost_90CI_time
    )+
    geom_point(data=cur_dat_estConcAtInfTime_sinceLastInfu,aes(y=Cc,x=infTimeSinceLastInfu),size=2)+
    scale_x_continuous(limits=c(0,70),breaks=(0:10)*7,labels=0:10)+
    scale_y_log10(limits=c(0.01,1000),breaks=c(0.01,0.1,1,10,100,500),labels=c("\u2264 0.01","0.1","1","10","100","500"))+
    scale_color_manual("",values=c("red","blue"))+
    ggtitle(title)+
    xlab("Weeks since last infusion")+
    ylab("VRC01 serum concentration (mcg/ml)")+
    facet_grid(~dose)+
    theme_bw()+
    theme(legend.position="bottom"
          ,axis.title = element_text(size=18)
          ,axis.text=element_text(size=14)
          ,legend.text = element_text(size=13)
          ,legend.title=element_text(size=15)
          ,strip.text = element_text(size=19)
          ,plot.title = element_text(size = 22, face = "bold",hjust = 0.5)
    )
}

## start plot
P_list <- plyr::dlply(conc_pred_medianWt %>%
                        group_by(ID,time) %>%
                        summarise(conc_lb=quantile(Cc,0.025)
                                  ,conc_ub=quantile(Cc,0.975)
                                  ,conc_median=median(Cc)) %>%
                        left_join(select(dat_cov,ID,study,dose))
                      ,plyr::.(study)
                      ,f_predInterval)

cairo_pdf("../figures/extFig5_predIntervalConcAndTime_allCase_Pfitter_byIC80Grp.pdf",height=10,width = 10,onefile = T)
grid.arrange(P_list[[1]],P_list[[2]],ncol=1)
dev.off()


##------------ fig5: Predicted neutralization coverage and (C, D) 
##            geometric mean predicted serum ID80 titer (PT80) against 
##            viruses circulating in each of the AMP trials for the bnAb
##            regimen PGT121LS + PGDM1400LS + VRC07-523LS 20+20+20 mg/kg 
##            delivered intravenously every 16 weeks and evaluated in study 
##            cohorts of the same sizes as the AMP trials. ---------------

## simulated PK for vrc07ls, pgdm1400 and pgt121, dose 20+20+20 and 40+40+40
dat_sim_16wkAR_2.5fd <- filter(dat_sim_ss_16wkAR,drug=="VRC07-523LS"|(drug%in%c("PGDM1400","PGT121")&fold==2.5),time<=112)


## data preparation for IC50 and IC80 (against 8 mAbs) 
dat_ic_8mAb_703 <- dat_ic_8mAb_703_raw %>% 
  select(protnum,isolate,poscrit,titer,mab_name)

dat_ic_8mAb_704 <- dat_ic_8mAb_704_raw %>% 
  select(protnum,isolate,poscrit,titer,mab_name)

dat_ic_8mAb <- bind_rows(dat_ic_8mAb_703,dat_ic_8mAb_704) %>% 
  filter(mab_name%in%c("PGDM1400","PGT121","VRC07-523LS","PGT121.414LS","VRC01")) %>%
  mutate(pub_id=gsub("HIV ","",isolate)
         ,pub_id=substr(pub_id,2,9)
         ,pub_id=gsub("_","-",pub_id)) %>% 
  left_join(select(dat_eliCtl,rx_code,pub_id)) %>% 
  filter(rx_code=="C3") %>% 
  mutate(mab_name=ifelse(mab_name=="PGT121.414LS","PGT121",mab_name)
         ,titer_num=as.numeric(ifelse(titer=="<<1.14e-02","0.0057"
                                      ,ifelse(titer==">>25.0","50"
                                              ,ifelse(titer=="<<4.57e-03","0.002285"
                                                      ,ifelse(titer=="<<0.011431","0.0057155",titer)))))) %>% 
  rename(study=protnum,drug=mab_name) %>% 
  select(study,isolate,poscrit,titer_num,drug)


## calculate coverage and make plots for coverage and titers 

gen_fig_coverage_titer <- function(dat_sim_pk,dat_ic80,fold, hl_fold,poscrit){
  
  if(poscrit=="50"){
    scale_secY <- 7
  }else{
    scale_secY <- 6
  }
  
  # calculate the mean concentration and then the coverage
  dat_IPRE_gm <-dat_sim_pk %>%
    mutate(IPRE=ifelse(Cc_ss<=0,0.0001,Cc_ss)) %>%
    group_by(drug,study,time,dose) %>%
    summarise(IPRE=exp(mean(log(IPRE)))) %>%
    ungroup() %>%
    mutate(dose=paste0(dose," mg/kg"))
  
  dat_cov_conc_all <- dat_IPRE_gm %>%
    left_join(dat_ic80) %>%
    mutate(cov_ind=ifelse(IPRE>titer_num*fold,1,0)
    )
  
  dat_cov_comb <- filter(dat_cov_conc_all,drug!="VRC01") %>%
    select(-IPRE,-titer_num) %>%
    mutate(drug=ifelse(drug=="VRC07-523LS","VRC07.523LS",drug)) %>%
    spread(key="drug",value="cov_ind") %>%
    mutate(ind_overlap=ifelse((!is.na(PGDM1400))&(!is.na(PGT121))&(!is.na(VRC07.523LS)),1,0)
           ,PGDM1400=ifelse(is.na(PGDM1400),0,PGDM1400)
           ,PGT121=ifelse(is.na(PGT121),0,PGT121)
           ,VRC07.523LS=ifelse(is.na(VRC07.523LS),0,VRC07.523LS)
           ,cov_sum=PGDM1400+PGT121+VRC07.523LS
           ,dose_n = gsub(" mg/kg","",dose)
           ,dose=paste0(dose_n,"+",dose_n,"+",dose_n," mg/kg")
    )
  
  dat_cov_comb_1Act <- dat_cov_comb %>%
    mutate(cov_1act = ifelse(cov_sum>=1,1,0)) %>%
    group_by(study,time,dose) %>%
    summarise(n_all=n()
              ,n_overlap=sum(ind_overlap)
              ,cov=mean(cov_1act)) %>%
    select(-n_all) %>%
    mutate(drug="PGT121 + PGDM1400 + VRC07-523LS"
           ,value_type="coverage"
           ,active="(1-Active)") %>%
    rename(value=cov)
  
  
  dat_cov_comb_2Act <- dat_cov_comb %>%
    mutate(cov_2act = ifelse(cov_sum>=2,1,0)) %>%
    group_by(study,time,dose) %>%
    summarise(n_all=n()
              ,n_overlap=sum(ind_overlap)
              ,cov=mean(cov_2act)) %>%
    select(-n_all) %>%
    mutate(drug="PGT121 + PGDM1400 + VRC07-523LS"
           ,value_type="coverage"
           ,active="(2-Active)") %>%
    rename(value=cov)
  
  dat_cov_comb_3Act <- dat_cov_comb %>%
    mutate(cov_3act = ifelse(cov_sum>=3,1,0)) %>%
    group_by(study,time,dose) %>%
    summarise(n_all=n()
              ,n_overlap=sum(ind_overlap)
              ,cov=mean(cov_3act)) %>%
    select(-n_all) %>%
    mutate(drug="PGT121 + PGDM1400 + VRC07-523LS"
           ,value_type="coverage"
           ,active="(3-Active)") %>%
    rename(value=cov)
  
  dat_cov_vrc01 <- filter(dat_cov_conc_all,drug=="VRC01") %>%
    group_by(study,drug,dose,time) %>%
    summarise(cov=mean(cov_ind)) %>%
    mutate(value_type="coverage") %>%
    rename(value=cov)
  
  poscrit.num <- as.numeric(poscrit)
  
  dat_nab_gm <- dat_IPRE_gm %>%
    left_join(filter(dat_ic_8mAb_gm,poscrit==poscrit.num)) %>% 
    mutate(nab=IPRE/titer_num
           ,nab=log10(nab)/scale_secY+0.5
           ,value_type=paste0("ID",poscrit)
    ) %>%
    rename(value=nab) %>% 
    select(-poscrit,-titer_num,-IPRE)
  
  
  if(poscrit=="50"){
    dat_plot <- bind_rows(dat_cov_vrc01,dat_cov_comb_1Act,dat_cov_comb_2Act,dat_cov_comb_3Act,dat_nab_gm) %>%
      mutate(fold=fold
             ,scale_secY=scale_secY
             ,route="IV"
             ,drug_dose=ifelse(is.na(active)
                               ,paste0(route," ",drug," ",dose)
                               ,paste0(route," ",drug," ",dose," ",active))
             ,type=ifelse(value_type=="coverage","Coverage","ID50")
             ,drug_dose_type=paste0(type,": ",drug_dose)
             ,drug_dose_type=factor(drug_dose_type,levels=c("Coverage: IV PGT121 + PGDM1400 + VRC07-523LS 20+20+20 mg/kg (1-Active)"
                                                            ,"Coverage: IV PGT121 + PGDM1400 + VRC07-523LS 20+20+20 mg/kg (2-Active)"
                                                            ,"Coverage: IV PGT121 + PGDM1400 + VRC07-523LS 20+20+20 mg/kg (3-Active)"
                                                            ,"Coverage: IV VRC01 30 mg/kg"
                                                            ,"ID50: IV PGT121 20 mg/kg"
                                                            ,"ID50: IV PGDM1400 20 mg/kg"
                                                            ,"ID50: IV VRC07-523LS 20 mg/kg"
                                                            ,"ID50: IV VRC01 30 mg/kg")))
  }else{
    dat_plot <- bind_rows(dat_cov_vrc01,dat_cov_comb_1Act,dat_cov_comb_2Act,dat_cov_comb_3Act,dat_nab_gm) %>%
      mutate(fold=fold
             ,scale_secY=scale_secY
             ,route="IV"
             ,drug_dose=ifelse(is.na(active)
                               ,paste0(route," ",drug," ",dose)
                               ,paste0(route," ",drug," ",dose," ",active))
             ,type=ifelse(value_type=="coverage","Coverage","ID80")
             ,drug_dose_type=paste0(type,": ",drug_dose)
             ,drug_dose_type=factor(drug_dose_type,levels=c("Coverage: IV PGT121 + PGDM1400 + VRC07-523LS 20+20+20 mg/kg (1-Active)"
                                                            ,"Coverage: IV PGT121 + PGDM1400 + VRC07-523LS 20+20+20 mg/kg (2-Active)"
                                                            ,"Coverage: IV PGT121 + PGDM1400 + VRC07-523LS 20+20+20 mg/kg (3-Active)"
                                                            ,"Coverage: IV VRC01 30 mg/kg"
                                                            ,"ID80: IV PGT121 20 mg/kg"
                                                            ,"ID80: IV PGDM1400 20 mg/kg"
                                                            ,"ID80: IV VRC07-523LS 20 mg/kg"
                                                            ,"ID80: IV VRC01 30 mg/kg")))
  }
  
  
  #---------- start coverage plot for comb and VRC01 -------------------
  p_list <- plyr::dlply(dat_plot,plyr::.(study),function(df){

    #second y axis scaler
    scale_secY <- unique(df$scale_secY)
    
    #average coverage
    dat_avg_8wk <- df %>%
      filter(time<=8*7,value_type=="coverage") %>%
      group_by(drug,active) %>%
      summarise(mean_cov_8wk=paste0(round(mean(value)*100,0),"%"))
    
    dat_avg_12wk <- df %>%
      filter(time<=12*7,value_type=="coverage") %>%
      group_by(drug,active) %>%
      summarise(mean_cov_12wk=paste0(round(mean(value)*100,0),"%"))
    
    dat_avg_16wk <- df %>%
      filter(time<=16*7,value_type=="coverage") %>%
      group_by(drug,active) %>%
      summarise(mean_cov_16wk=paste0(round(mean(value)*100,0),"%"))
    
    #footnote
    cap_avg_cov_vrc01 <- paste0("VRC01: Average coverage in the first 8, 12 and 16 weeks are "
                                ,dat_avg_8wk$mean_cov_8wk[dat_avg_8wk$drug=="VRC01"]
                                , ", ",dat_avg_12wk$mean_cov_12wk[dat_avg_12wk$drug=="VRC01"]
                                ,", and ", dat_avg_16wk$mean_cov_16wk[dat_avg_16wk$drug=="VRC01"], ", respectively.")
    
    cap_avg_cov_comb_1Act <- paste0("PGT121 + PGDM1400 + VRC07-523LS (1-Active): Average coverage in the first 8, 12 and 16 weeks\nare "
                                    ,dat_avg_8wk$mean_cov_8wk[dat_avg_8wk$drug=="PGT121 + PGDM1400 + VRC07-523LS"&dat_avg_8wk$active=="(1-Active)"]
                                    , ", ",dat_avg_12wk$mean_cov_12wk[dat_avg_12wk$drug=="PGT121 + PGDM1400 + VRC07-523LS"&dat_avg_12wk$active=="(1-Active)"]
                                    ,", and ", dat_avg_16wk$mean_cov_16wk[dat_avg_16wk$drug=="PGT121 + PGDM1400 + VRC07-523LS"&dat_avg_16wk$active=="(1-Active)"], ", respectively.")
    
    cap_avg_cov_comb_2Act <- paste0("PGT121 + PGDM1400 + VRC07-523LS (2-Active): Average coverage in the first 8, 12 and 16 weeks\nare "
                                    ,dat_avg_8wk$mean_cov_8wk[dat_avg_8wk$drug=="PGT121 + PGDM1400 + VRC07-523LS"&dat_avg_8wk$active=="(2-Active)"]
                                    , ", ",dat_avg_12wk$mean_cov_12wk[dat_avg_12wk$drug=="PGT121 + PGDM1400 + VRC07-523LS"&dat_avg_12wk$active=="(2-Active)"]
                                    ,", and ", dat_avg_16wk$mean_cov_16wk[dat_avg_16wk$drug=="PGT121 + PGDM1400 + VRC07-523LS"&dat_avg_16wk$active=="(2-Active)"], ", respectively.")
    
    cap_avg_cov_comb_3Act <- paste0("PGT121 + PGDM1400 + VRC07-523LS (3-Active): Average coverage in the first 8, 12 and 16 weeks\nare "
                                    ,dat_avg_8wk$mean_cov_8wk[dat_avg_8wk$drug=="PGT121 + PGDM1400 + VRC07-523LS"&dat_avg_8wk$active=="(3-Active)"]
                                    , ", ",dat_avg_12wk$mean_cov_12wk[dat_avg_12wk$drug=="PGT121 + PGDM1400 + VRC07-523LS"&dat_avg_12wk$active=="(3-Active)"]
                                    ,", and ", dat_avg_16wk$mean_cov_16wk[dat_avg_16wk$drug=="PGT121 + PGDM1400 + VRC07-523LS"&dat_avg_16wk$active=="(3-Active)"], ", respectively.")
    
    
    cap_avg_cov <- paste0(cap_avg_cov_comb_1Act,"\n",cap_avg_cov_comb_2Act,"\n",cap_avg_cov_comb_3Act,"\n",cap_avg_cov_vrc01)
    
    # y label and plot title
    y.lab <- paste0("Predicted coverage (ID",poscrit," > ",fold,")")
    
    if(hl_fold==1){
      title <- paste0("IV PGT121 + PGDM1400 + VRC07-523LS 20+20+20 mg/kg\nHVTN",unique(df$study)," placebo viruses (m=",unique(df$n_overlap)[2],")")
    }else{
      title <- paste0("IV PGT121LS + PGDM1400LS + VRC07-523LS 20+20+20 mg/kg\nHVTN",unique(df$study)," placebo viruses (m=",unique(df$n_overlap)[2],")\nLS half-life ",hl_fold,"-fold over WT for PGT121 and PGDM1400")
    }
    
    
    #plot 
    ggplot(df,aes(x=time,y=value))+
      geom_line(size=1,aes(color=drug_dose_type,linetype=drug_dose_type))+
      scale_x_continuous(breaks=0:12*2*7,labels=0:12*2)+
      scale_y_continuous(breaks=seq(0.1,1,by=0.1),limits = c(0,1.01)
                         ,labels = function(b) { paste0(round(b * 100, 0), "%")}
                         ,expand = c(0, 0)
                         ,sec.axis = sec_axis(~.
                                              ,name=paste0("Geometric mean of predicted ID",poscrit," titers")
                                              ,breaks=log10(c(1,5,10,50,100,200,500))/scale_secY+0.5
                                              ,labels=c("1","5","10","50","100","200","500"))
                         
      )+
      xlab("Weeks post product administration at steady state")+
      ylab(y.lab)+
      ggtitle(title)+
      labs(caption=cap_avg_cov)+
      scale_color_manual("",values=c("darkblue","blue","lightblue","red","forestgreen","brown","purple","orange"))+
      scale_linetype_manual("",values=c(1,1,1,1,3,3,3,3)
      )+
      guides(color=guide_legend(ncol=2,nrow=4))+
      theme_bw()+
      theme(axis.title=element_text(size=11)
            ,axis.text = element_text(size=10)
            ,plot.margin=unit(c(1.5,1.5,0.5,1.5), "lines")
            ,plot.title = element_text(hjust=0,size=12,vjust=6,face="bold")
            ,plot.caption = element_text(face="bold",hjust=0,size=8)
            ,plot.caption.position =  "plot"
            # ,plot.title.position = "plot"
            ,axis.title.y.right = element_text(vjust=2)
            ,legend.text = element_text(size=6.4)
            ,legend.position="bottom"
            ,legend.box="vertical"
            ,legend.spacing.y = unit(0, 'cm'))
    
    
  })
  return(p_list)
  
}

## plot data
dat_ic_8mAb_gm <- dat_ic_8mAb %>% 
  group_by(study,poscrit,drug) %>% 
  summarise(titer_num=exp(mean(log(titer_num))))

dat_sim_16wkAR_2.5fd_d20 <- dat_sim_16wkAR_2.5fd %>% 
  bind_rows(filter(CcPred_vrc01_ss,time<=112)) %>% 
  filter(dose%in%c(20,30))


## start plot
p_list <- gen_fig_coverage_titer (dat_sim_pk=dat_sim_16wkAR_2.5fd_d20
                                    ,dat_ic80=filter(dat_ic_8mAb,poscrit==80)
                                    ,fold=200
                                    ,hl_fold=2.5
                                    ,poscrit="80"
                                    
                                    
)

cairo_pdf("../figures/fig5_panelA_cov_ID80_iv_placeboVirus_2.5fd_16wkAR_cut200.pdf",height = 7,width=6,onefile = T)
p_list
dev.off()


## ID80/ID50 over time curve for individual bNAb and combined.
dat_ic_8mAb_gm <- dat_ic_8mAb %>% 
  group_by(study,poscrit,drug) %>% 
  summarise(titer_num=exp(mean(log(titer_num))))

## ID80 plot function
gen_fig_nabTiterCurve <- function(dat_sim_pk,hl_fold,poscrit){
  dat_IPRE_gm <-dat_sim_pk %>%
    mutate(IPRE=ifelse(Cc_ss<=0,0.0001,Cc_ss)) %>%
    group_by(drug,study,time,dose) %>%
    summarise(IPRE=exp(mean(log(IPRE)))) %>%
    ungroup() %>%
    mutate(dose=paste0(dose," mg/kg"))
  
  poscrit.num <- as.numeric(poscrit)
  dat_nab_gm <- dat_IPRE_gm %>%
    left_join(filter(dat_ic_8mAb_gm,poscrit==poscrit.num)) %>% 
    mutate(nab=IPRE/titer_num
    ) %>% 
    select(drug,study,time,dose,poscrit,nab)
  
  dat_nab_comb <- dat_nab_gm %>% 
    filter(drug!="VRC01") %>% 
    spread(key="drug",value="nab") %>% 
    mutate(nab_bliss=ifelse(poscrit==50
                            ,calc_3bnab_BHtiter(PGDM1400, PGT121, `VRC07-523LS`)
                            ,calc_3bnab_BHtiter(PGDM1400*4, PGT121*4, `VRC07-523LS`*4,titer_target=0.8))) %>% 
    rename(nab=nab_bliss) %>% 
    mutate(drug="IV PGT121LS + PGDM1400LS + VRC07-523LS 20+20+20 mg/kg"
           ,type="combine") %>% 
    select(drug,study,time,dose,poscrit,nab,type) 
  
  dat_nab_plot <- dat_nab_gm %>% 
    mutate(drug=ifelse(drug=="PGT121","PGT121LS"
                       ,ifelse(drug=="PGDM1400","PGDM1400LS",drug))
           ,drug=paste0("IV ",drug," ",dose)
           ,type="individual") %>% 
    bind_rows(dat_nab_comb) %>% 
    mutate(drug=factor(drug,levels=c("IV PGT121LS 20 mg/kg"
                                     ,"IV PGDM1400LS 20 mg/kg"
                                     ,"IV VRC07-523LS 20 mg/kg"
                                     ,"IV VRC01 30 mg/kg"
                                     ,"IV PGT121LS + PGDM1400LS + VRC07-523LS 20+20+20 mg/kg")))
  
  p_list <- plyr::dlply(dat_nab_plot,plyr::.(study),function(df){
    title <- paste0("HVTN ",unique(df$study))
    ggplot(df,aes(x=time,y=nab))+
      geom_line(size=1,aes(color=drug,linetype=drug))+
      scale_x_continuous(breaks=0:12*2*7,labels=0:12*2)+
      scale_y_log10(breaks=c(1,10,100,300,1000,2000,4000))+
      xlab("Weeks post product administration at steady state")+
      ylab(paste0("Geometric mean of PT",poscrit))+
      scale_color_manual("",values=c("forestgreen","blue","purple","orange","red"))+
      scale_linetype_manual("",values=c(2,2,2,2,1)
      )+
      guides(color=guide_legend(ncol=2,nrow=4))+
      ggtitle(title)+
      theme_bw()+
      theme(axis.title=element_text(size=14)
            ,axis.text = element_text(size=12)
            ,plot.margin=unit(c(1.5,1.5,0.5,1.5), "lines")
            ,plot.title = element_text(size=15)
            ,legend.text = element_text(size=7)
            ,legend.position="bottom"
            ,legend.box="vertical"
            ,legend.spacing.y = unit(0, 'cm'))
  }) 
  
  p_list
  
}

# start plot
cairo_pdf("../figures/fig5_panelB_ID80Curve_2.5fd_16wkAR.pdf",height = 6,width=6,onefile = T)
gen_fig_nabTiterCurve (dat_sim_pk=dat_sim_16wkAR_2.5fd_d20
                       ,hl_fold=2.5
                       ,poscrit="80"
)
dev.off()


##------------ Extended Fig8: Predicted neutralization coverage and (C, D) 
##            geometric mean predicted serum ID50 titer (PT50) against 
##            viruses circulating in each of the AMP trials for the bnAb
##            regimen PGT121LS + PGDM1400LS + VRC07-523LS 20+20+20 mg/kg 
##            delivered intravenously every 16 weeks and evaluated in study 
##            cohorts of the same sizes as the AMP trials. ---------------

## coverage plot
p_list <- gen_fig_coverage_titer (dat_sim_pk=dat_sim_16wkAR_2.5fd_d20
                                  ,dat_ic80=filter(dat_ic_8mAb,poscrit==50)
                                  ,fold=600
                                  ,hl_fold=2.5
                                  ,poscrit="50"
                                  
                                  
)

cairo_pdf("../figures/extFig8_panelA_cov_ID50_iv_placeboVirus_2.5fd_16wkAR_cut600.pdf",height = 7,width=6,onefile = T)
p_list
dev.off()

# ID50 plot
cairo_pdf("../figures/extFig8_panelB_ID50Curve_2.5fd_16wkAR.pdf",height = 6,width=6,onefile = T)
gen_fig_nabTiterCurve (dat_sim_pk=dat_sim_16wkAR_2.5fd_d20
                       ,hl_fold=2.5
                       ,poscrit="50"
                       
                       
)
dev.off()


##---------- fig6/Extended fig9: Predicted serum ID80/ID50 titer (PT80/PT50)-predicted prevention efficacy 
##           over time in the context of viruses circulating in each of the AMP trials
##          for the bnAb regimen PGDM1400LS + PGT121LS + VRC07-523LS at 20+20+20 mg/kg 
##          or 40+40+40 mg/kg, delivered intravenously every 16 weeks and evaluated in 
##          study cohorts of the same sizes as the AMP trials.------------------------

## function to estimate PE 

#original pegu_fn from Allan's code
log_logistic = function(b,c,d,e,f) {
  function(dose) { c + (d - c)/((1 + exp(b * (log(dose/e))))^f) }
}

# 2 times of titers compared to Allan's original numbers
# Code from Bryan for fitting Pegu efficacy curve (needs library drc which doesn't load in this version of R)
# ID80
pegu_dat = data.frame(
  titer = c(22*2, 48*2, 103*2),
  neut = c(0.5, 0.75, 0.9)
)
empirical_pd_fit = drc::drm(neut ~ titer, data = pegu_dat,
                            fct = drc::LL.5(names = c("slope", "lower", "upper", "inflection", "asymmetry"),
                                            fixed = c(NA, 0, 1, NA, NA)))
# Coefficients:
# slope:(Intercept)  inflection:(Intercept)
# -1.4649                 52.0219
# asymmetry:(Intercept)
# 0.8418

# ID50
pegu_dat = data.frame(
  titer = c(87*2, 208*2, 437*2),
  neut = c(0.5, 0.75, 0.9)
)
empirical_pd_fit = drc::drm(neut ~ titer, data = pegu_dat,
                            fct = drc::LL.5(names = c("slope", "lower", "upper", "inflection", "asymmetry"),
                                            fixed = c(NA, 0, 1, NA, NA)))

# Coefficients:
#      slope:(Intercept)  inflection:(Intercept)
#                -1.6683                412.1532
#  asymmetry:(Intercept)
#                 0.4197

pegu_fn_id80 = log_logistic(-1.4649, 0, 1, 52.0219, 0.8418)
pegu_fn_id50 = log_logistic(-1.6683, 0, 1, 412.1532, 0.4197)


## PE for VRC01 dose 30
dat_pe_vrc01 <- plyr::ddply(CcPred_vrc01_ss,plyr::.(id),function(df){
  dat_pe_ind <- df %>%
    left_join(dat_ic_8mAb) %>%
    mutate(nab=Cc_ss/titer_num
           ,pe=ifelse(poscrit==50,pegu_fn_id50(nab),pegu_fn_id80(nab))) %>%
    select(id,time,dose,drug,isolate,pe,poscrit,study) %>%
    group_by(dose,id,time,poscrit,study) %>%
    summarise(pe=mean(pe)
              ,n_virus=n())
  return(dat_pe_ind)
})

## PE for VRC01 dose 40
dat_pe_vrc01_d40 <- plyr::ddply(CcPred_vrc01_ss_d40,plyr::.(id),function(df){
  dat_pe_ind<- df %>%
    left_join(dat_ic_8mAb) %>%
    mutate(nab=Cc_ss/titer_num
           ,pe=ifelse(poscrit==50,pegu_fn_id50(nab),pegu_fn_id80(nab))) %>%
    select(id,time,dose,drug,isolate,pe,poscrit,study) %>%
    group_by(dose,id,time,poscrit,study) %>%
    summarise(pe=mean(pe)
              ,n_virus=n())
  return(dat_pe_ind)
  
})

## PE for VRC07ls, pgdm1400 and pgt121 and PE figures
gen_fig_pe <- function(dat_sim,hl_fold,dat_pe_vrc01_time,AR_week,dat_pe_vrc01_d40_time){
  dat_pe_individal_list <- plyr::ddply(dat_sim,plyr::.(dose,id,time),function(df){
    print(unique(df$id))
    print(unique(df$time))
    # browser()
    pe <- df %>%
      left_join(dat_ic_8mAb) %>%
      mutate(nab=Cc_ss/titer_num) %>%
      select(id,time,dose,drug,isolate,nab,poscrit,study) %>%
      spread(key="drug",value="nab") %>%
      na.omit() %>%
      group_by(poscrit,isolate) %>%  # for bh can't calculate some matrix, so calculate by each row
      mutate(nab_add=PGDM1400+PGT121+`VRC07-523LS`
             ,nab_bliss=ifelse(poscrit==50
                               ,calc_3bnab_BHtiter(PGDM1400, PGT121, `VRC07-523LS`)
                               ,calc_3bnab_BHtiter(PGDM1400*4, PGT121*4, `VRC07-523LS`*4,titer_target=0.8))
             ,pe_add=ifelse(poscrit==50,pegu_fn_id50(nab_add),pegu_fn_id80(nab_add))
             ,pe_bh=ifelse(poscrit==50,pegu_fn_id50(nab_bliss),pegu_fn_id80(nab_bliss))) %>%
      group_by(dose,id,time,poscrit,study) %>%
      summarise(pe_add=mean(pe_add)
                ,pe_bh=mean(pe_bh)
                ,n_virus=n())
    # if(class(pe_try)=="try-error")
    #  {pe <- data.frame(dose=unique(df$dose),id=unique(df$id))}
    return(pe)
  })
  
  dat_pe_individal <- bind_rows(dat_pe_individal_list)

  dat_pe_vrc01_sub <- bind_rows(dat_pe_vrc01_time,dat_pe_vrc01_d40_time) %>%
    mutate(dose=ifelse(dose==30,"IV VRC01 30 mg/kg","IV VRC01 40 mg/kg"))
  
  dat_pe <- dat_pe_individal %>%
    mutate(dose=ifelse(dose==20,"IV 20+20+20 mg/kg","IV 40+40+40 mg/kg")) %>%
    rename(pe=pe_bh) %>%
    select(-pe_add) %>%
    bind_rows(dat_pe_vrc01_sub) %>%
    group_by(time,poscrit,dose,study,n_virus) %>%
    summarise(median_pe=median(pe,na.rm=T)
              ,lb_pe=quantile(pe, 0.025,na.rm=T)
              ,ub_pe=quantile(pe,0.975,na.rm=T)
    )
  
  
  p_list <- plyr::dlply(dat_pe,plyr::.(poscrit,study),function(df){
    # browser()
    
    if(hl_fold==1){
      title <- paste0("Predicted PE at steady state of IV PGT121 + PGDM1400 + VRC07-523LS\nagainst HVTN",unique(df$study)," placebo viruses (m = ",unique(df$n_virus),"), based on PE vs. ID",unique(df$poscrit)," curve in\nAMP and Pegu et al..")
    }else{
      title <- paste0("Predicted PE at steady state of IV PGT121LS + PGDM1400LS + VRC07-523LS\nagainst HVTN",unique(df$study)," placebo viruses (m = ",unique(df$n_virus),"), based on PE vs. ID",unique(df$poscrit)," curve in AMP\nand Pegu et al., and assuming LS half-life ",hl_fold,"-fold over WT for PGT121 and PGDM1400")
    }
    
    ggplot(df,aes(x=time,y=median_pe))+
      geom_line(aes(color=dose),size=0.8)+
      geom_ribbon(aes(ymin=lb_pe,ymax=ub_pe,fill=dose),alpha=0.1)+
      scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1),labels =scales:: percent)+
      scale_x_continuous(breaks=0:12*2*7,labels=0:12*2)+
      scale_color_manual("",values=c("blue","dark blue","orange","red"))+
      scale_fill_manual("",values=c("blue","dark blue","orange","red"))+
      ylab("Predicted Prevention Efficacy")+
      xlab("Weeks post product administration at steady state")+
      ggtitle(title)+
      theme_bw(base_size = 14)+
      theme(title=element_text(size=11)
            # ,legend.position = c(0.2,0.4)
            ,axis.title=element_text(size=15)
      )
    
  })
  
  p_list
  
}


## start plot 
p_list <- gen_fig_pe(dat_sim=dat_sim_16wkAR_2.5fd
                       ,hl_fold=2.5
                       ,dat_pe_vrc01_time=filter(dat_pe_vrc01,time<=112)
                       ,dat_pe_vrc01_d40_time=filter(dat_pe_vrc01_d40,time<=112)
                       ,AR_week=16)

cairo_pdf("../figures/fig6_extFig9_pe_combNabBH_2.5fd_16wkAR.pdf",height = 5.5,width=8,onefile = T)
p_list
dev.off()


##------------ Extended fig10: Predicted serum ID80 titer (PT80)-predicted prevention 
##            efficacy over time in the context of viruses circulating in each of the AMP 
##            trials for the bnAb regimen PGDM1400LS + PGT121LS + VRC07-523LS at 20+20+20 mg/kg 
##            or 40+40+40 mg/kg, delivered intravenously every 24 weeks and evaluated in study cohorts
##            of the same sizes as the AMP trials.  ------------------------------------

## pk data for 24 week infusion interval
dat_sim_24wkAR_2.5fd <- filter(dat_sim_ss_24wkAR,drug=="VRC07-523LS"|(drug%in%c("PGDM1400","PGT121")&fold==2.5))

## start plot 
p_list <- gen_fig_pe(dat_sim=dat_sim_24wkAR_2.5fd
                       ,hl_fold=2.5
                       ,dat_pe_vrc01_time=filter(dat_pe_vrc01,time<=168)
                       ,dat_pe_vrc01_d40_time=filter(dat_pe_vrc01_d40,time<=168)
                       ,AR_week=24)

cairo_pdf("../figures/extFig10_pe_combNabBH_2.5fd_24wkAR.pdf",height = 5.5,width=8,onefile = T)
p_list
dev.off()

















