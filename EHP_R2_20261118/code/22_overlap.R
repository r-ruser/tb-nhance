options(survey.lonely.psu="adjust")
.libPaths(c("TB_metabolome_exposome_20261115_2130/analysis_MWAS_ExWAS/next_stage_validation_20261116/Rlib","TB_metabolome_exposome_20261115_2130/analysis_MWAS_ExWAS/Rlib",.libPaths()))
library(survey);library(svyVGAM)
out<-"EHP_R2_20261118";o<-file.path(out,"B_results")
saveout<-function(z,n)write.csv(z,file.path(o,paste0(n,".csv")),row.names=FALSE,na="")
a<-readRDS("NHANES_TB_Stage3_Inference_Audit_20261116/数据/rebuilt_private.rds")$data
reg<-read.csv(file.path(out,"Source_data/01_expanded_registry_188.csv"))
a$race6<-factor(a$RIDRETH3,levels=c(3,1,2,4,6,7),labels=c("NH White","Mexican American","Other Hispanic","NH Black","NH Asian","Other or multiracial"))
a$duration<-factor(ifelse(a$birth=="NonUS"&a$DMDYRSUS%in%1:9,ifelse(a$DMDYRSUS<=3,"<10 y",ifelse(a$DMDYRSUS<=5,"10-<20 y","20+ y")),NA))
a$g<-factor(a$group,levels=c("Q0T0","Q1T0","Q0T1","Q1T1"))
cv<-c("age","sex_male","race","birth","smoking","educ","PIR")
flags<-list();counts<-list();crude<-list()
for(i in seq_len(nrow(reg))){
 r<-reg[i,];v<-r$variable; C<-cv
 if(r$smoking_biomarker)C<-setdiff(C,"smoking");if(r$urinary)C<-c(C,"logucr")
 x<-a[[v]]
 if(r$detection_category %in% "<50% detected"){lc<-a[[r$lc_variable]];x<-ifelse(lc %in% 0:1,as.integer(lc==0),NA)}else if(r$transform=="log2")x<-log2(ifelse(x>0,x,NA))
 b<-a[,unique(c("SEQN","SDMVSTRA","SDMVPSU",C,"race6","duration","y","g","smoking","birth","age")),drop=F]
 b$x<-x;b$wt<-a[[r$weight]]
 b<-subset(b,is.finite(wt)&wt>0)
 d<-svydesign(~SDMVPSU,strata=~SDMVSTRA,weights=~wt,data=b,nest=TRUE)
 cc<-subset(d,age>=20 & !is.na(g) & is.finite(x) & complete.cases(b[,C,drop=F]))
 dd<-subset(cc,!is.na(y))

 z<-dd$variables
 binary<-r$detection_category %in% "<50% detected"
 # For a continuous exposure assess support across its median as a diagnostic only.
 med<-if(nrow(z))median(z$x) else NA_real_
 xx<-if(binary)z$x else as.integer(z$x>med)
 cats<-intersect(c("sex_male","race","birth","smoking","educ"),C)
 parts<-lapply(cats,function(cov){
  g<-droplevels(factor(z[[cov]]))
  tt<-table(g,factor(xx,levels=0:1))
  if(!nrow(tt))return(NULL)
  data.frame(variable=v,covariate=cov,level=rownames(tt),lower=tt[,1],upper=tt[,2])
 })
 cell<-do.call(rbind,parts)
 counts[[v]]<-cell
 status<-if(length(unique(z$x))<2)"no contrast" else if(any(cell$lower<5|cell$upper<5))"positivity_limited" else "no flag"
 flags[[v]]<-data.frame(variable=v,n=nrow(z),comparison=if(binary)"detection" else "above vs at/below median",cutoff_model_scale=if(binary)NA else med,minimum=if(nrow(cell))min(cell$lower,cell$upper)else NA,status=status)
 if(v %in% c("LBXBCD","LBX2DF","LBXVFN")){
  fit<-svyglm(y~x,dd,family=quasibinomial());b0<-coef(fit)[2];se<-sqrt(vcov(fit)[2,2])
  crude[[v]]<-data.frame(variable=v,n=nrow(z),events=sum(z$y),OR=exp(b0),low=exp(b0-qt(.975,17)*se),high=exp(b0+qt(.975,17)*se),se=se)
 }
}
for(k in c("flags","counts","crude"))saveout(do.call(rbind,get(k)),paste0("r2_",k))
