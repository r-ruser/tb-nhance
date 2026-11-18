options(survey.lonely.psu='adjust')
b<-'.';o<-file.path(b,'IJHEH_Exposure_Nutrition_Phenotype_20261118');dir.create(o,showWarnings=FALSE);for(n in c('Source_data','code','logs','Figures'))dir.create(file.path(o,n),showWarnings=FALSE)
.libPaths(c(file.path(b,'TB_metabolome_exposome_20261115_2130/analysis_MWAS_ExWAS/Rlib'),.libPaths()));library(survey)
a<-readRDS(file.path(b,'NHANES_TB_Stage3_Inference_Audit_20261116/数据/rebuilt_private.rds'))$data
a$BMI<-as.numeric(a$BMXBMI);a$albumin<-as.numeric(a$LBXSAL);a$vitD<-as.numeric(a$LBXVIDMS)
now<-ifelse(a$WHD020%in%c(7777,9999,777,999),NA,a$WHD020);past<-ifelse(a$WHD050%in%c(7777,9999,777,999),NA,a$WHD050);loss<-past-now
unintent<-ifelse(is.na(loss),NA,ifelse(loss<10,0,ifelse(a$WHQ060==2,1,ifelse(a$WHQ060==1,0,NA))))
a$risk<-ifelse(!is.na(a$BMI)&!is.na(unintent),as.integer(a$BMI<18.5|unintent==1),NA)
vv<-c('BMI','albumin','vitD','risk','LBXSCU','LBXSSE','LBXSZN','LBXSIR','URXUIO')
labs<-c('BMI (kg/m²)','Serum albumin (g/dL)','25-hydroxyvitamin D (nmol/L)','Low-BMI/unintentional-weight-loss proxy (%)','Serum copper (ug/dL)','Serum selenium (ug/L)','Serum zinc (ug/dL)','Serum iron (ug/dL)','Urinary iodine (ug/L)')
rr<-list()
for(j in seq_along(vv)){v<-vv[j];wt<-if(v%in%c('LBXSCU','LBXSSE','LBXSZN','URXUIO'))'WTSA2YR'else'WTMEC2YR';ii<-which(is.finite(a[[wt]])&a[[wt]]>0);a$xx<-a[[v]];de<-svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=as.formula(paste0('~',wt)),nest=TRUE,data=a[ii,])
 for(g in c('Q0T0','Q1T0','Q0T1','Q1T1')){dd<-subset(de,age>=20&group==g&is.finite(xx));est<-svymean(~xx,dd);bb<-unname(coef(est));se<-unname(SE(est));w<-as.numeric(weights(dd));rr[[length(rr)+1]]<-data.frame(variable=v,label=labs[j],weight=wt,group=g,n=nrow(dd$variables),risk_n=if(v=='risk')sum(dd$variables$risk)else NA,value=bb,low=bb-qt(.975,17)*se,high=bb+qt(.975,17)*se,kish_ess=sum(w)^2/sum(w^2),inference_df=17)}
}
r<-do.call(rbind,rr);stopifnot(nrow(r)==36);stopifnot(r$risk_n[r$variable=='risk'&r$group=='Q1T0']==15,r$risk_n[r$variable=='risk'&r$group=='Q0T1']==19)
write.csv(r,file.path(o,'Source_data/18_nutrition_four_state_descriptive.csv'),row.names=FALSE);capture.output(sessionInfo(),file=file.path(o,'logs/nutrition_descriptive_sessionInfo.txt'));print(r[r$variable%in%c('BMI','albumin','vitD','risk')&r$group%in%c('Q1T0','Q0T1'),]);cat('Completed descriptive profiles; no exposure-effect models changed.\n')
