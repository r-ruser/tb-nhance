# Run from the NHANES project root, passing the public CDC ZIP as argument 1.
.libPaths(c("TB_metabolome_exposome_20261115_2130/analysis_MWAS_ExWAS/Rlib",.libPaths()))
library(survey)
out<-"IDP_TBESC_20261122/B_results"
save_table<-function(x,n)write.csv(x,file.path(out,paste0("external_",n,".csv")),row.names=FALSE,na="")
archive<-commandArgs(TRUE)[1]
read_part<-function(name)read.delim(unz(archive,paste0("PartA03062025/",name,".csv")),sep="|",na.strings="",stringsAsFactors=FALSE)
m<-read_part("PartAMain"); labs<-subset(read_part("AllCombinedLabs"),d_main_set==1)
stopifnot(nrow(m)==22020,!anyDuplicated(m$studyid),!anyDuplicated(labs$studyid),all(labs$studyid %in% m$studyid))
a<-merge(m,labs,by="studyid",all.x=TRUE)
a$qvalid<-a$qft_pn %in% 0:1
a$tvalid<-a$tst_in_window %in% 1 & !is.na(a$tst_result) & a$tst_result>=0
a$paired<-as.integer(a$qvalid & a$tvalid)
a$state<-ifelse(a$paired==0,NA,paste0("Q",a$qft_pn,"T",as.integer(a$tst_result>=10)))
a$foreign<-ifelse(is.na(a$birth_country),NA,as.integer(a$birth_country!="USA"))
a$agegroup<-factor(a$d_age_cat)
a$sex<-factor(ifelse(a$d_sex %in% 1:2,a$d_sex,NA))
a$site<-factor(a$d_site)
a$BCG<-ifelse(a$bcg %in% 0:1,a$bcg,NA)
a$contact<-ifelse(a$close_contact %in% 0:1,a$close_contact,NA)
a$smoke<-ifelse(a$smoked_100==0,0,ifelse(a$current_smoker %in% 0:1,a$current_smoker,NA))
a$bmi_raw<-703.0696*a$weight_lb/a$height_in^2
a$bmi5<-ifelse(a$height_in>=48 & a$height_in<=84 & a$weight_lb>=60 & a$weight_lb<=600 & a$bmi_raw>=12 & a$bmi_raw<=80,a$bmi_raw/5,NA)
a$income3<-ifelse(a$income %in% 1:3,"USD 1000 or less",ifelse(a$income %in% 4:6,"USD 1001-4000",ifelse(a$income %in% 7:8,"Above USD 4000","Unknown")))
a$direction<-ifelse(a$state %in% c("Q0T1","Q1T0"),as.integer(a$state=="Q1T0"),NA)
a$in25<-a$d_age_cat %in% 6:8
b<-subset(a,in25)
stopifnot(nrow(b)==14233,sum(a$paired)==21002)
save_table(as.data.frame(with(b,table(qvalid,tvalid))),"availability")
save_table(as.data.frame(with(b,table(state,useNA="always"))),"states")
save_table(data.frame(phase=c("All enrolled","Age suppressed","Age 15-24","Age 25+","Age 25+ main set","Age 25+ paired","Age 25+ discordant","Age 25+ missing birthplace","Age 25+ both anthropometric fields present but excluded by range checks"),n=c(nrow(a),sum(is.na(a$d_age_cat)),sum(a$d_age_cat %in% 5),nrow(b),sum(b$d_main_set %in% 1),sum(b$paired),sum(!is.na(b$direction)),sum(is.na(b$foreign)),sum(!is.na(b$weight_lb)&!is.na(b$height_in)&is.na(b$bmi5)))),"flow")
tb<-svydesign(~site,weights=~1,data=b)
nh<-readRDS(Sys.glob("NHANES_TB_Stage3*/数据/rebuilt_private.rds")[1])$data
nh<-subset(nh,WTMEC2YR>0);nh$paired<-as.integer(!is.na(nh$group));nh$foreign<-as.integer(nh$birth=="NonUS")
nd<-svydesign(~SDMVPSU,strata=~SDMVSTRA,weights=~WTMEC2YR,nest=TRUE,data=nh)
rates<-function(des,name){
 res<-list()
 one<-function(dd,label){
  z<-svymean(~paired,dd);p<-unname(coef(z));s<-as.numeric(SE(z));df<-degf(des)
  ci<-plogis(qlogis(p)+c(-1,1)*qt(.975,df)*s/(p*(1-p)))
  data.frame(cohort=name,contrast=label,n=nrow(dd$variables),paired=sum(dd$variables$paired),estimate=100*p,low=100*ci[1],high=100*ci[2],df=df)
 }
 res[[1]]<-one(des,"All");res[[2]]<-one(subset(des,foreign==0),"US-born");res[[3]]<-one(subset(des,foreign==1),"Foreign-born")
 dd<-subset(des,!is.na(foreign));fit<-svyglm(paired~foreign,dd)
 k<-unname(coef(fit)[2]);s<-sqrt(vcov(fit)[2,2]);df<-degf(des)
 res[[4]]<-data.frame(cohort=name,contrast="Foreign-born minus US-born (pp)",n=nrow(dd$variables),paired=sum(dd$variables$paired),estimate=100*k,low=100*(k-qt(.975,df)*s),high=100*(k+qt(.975,df)*s),df=df)
 do.call(rbind,res)
}
save_table(rbind(rates(tb,"TBESC enrolled"),rates(subset(nd,age>=25 & !(TBQ070 %in% 1)),"NHANES no documented reaction"),rates(subset(nd,age>=25 & TBQ070 %in% 2),"NHANES documented no reaction")),"completion")
fitrow<-function(fit,term,design,name,odds=TRUE){
 k<-unname(coef(fit)[term]);se<-sqrt(vcov(fit)[term,term]);q<-qt(.975,degf(design));z<-c(k,k-q*se,k+q*se)
 if(odds)z<-exp(z) else z<-100*z
 data.frame(model=name,term=term,n=nrow(design$variables),events=sum(model.response(model.frame(fit))),sites=degf(design)+1,estimate=z[1],low=z[2],high=z[3])
}
d<-subset(tb,complete.cases(foreign,sex,agegroup))
sfit<-svyglm(paired~foreign+sex+agegroup+site,d)
res<-list(fitrow(sfit,"foreign",d,"TBESC within-site adjusted difference (pp)",FALSE))
d<-subset(tb,complete.cases(direction,bmi5,agegroup,sex,foreign,BCG,contact,smoke))
for(full in c(FALSE,TRUE)){
 f<-if(full)direction~bmi5+foreign+sex+agegroup+site+smoke+BCG+contact else direction~bmi5+foreign+sex+agegroup+site+smoke
 fit<-svyglm(f,d,family=quasibinomial())
 for(v in intersect(c("bmi5","foreign","BCG","contact"),names(coef(fit))))res[[length(res)+1]]<-fitrow(fit,v,d,if(full)"Discordance plus BCG and TB contact" else "Discordance matched base")
}
save_table(do.call(rbind,res),"models")
summ<-function(dd,groups){
 do.call(rbind,lapply(split(dd,groups),function(g)data.frame(n=nrow(g),paired=sum(g$paired),pct=100*mean(g$paired))))
}
inc<-summ(b,b$income3);save_table(cbind(income=rownames(inc),inc),"income")
capture.output(sessionInfo(),file=file.path(out,"external_session.txt"))
# Sensitivity to site mix and to the inclusion of US territories in the US-born category.
leave<-sapply(levels(b$site),function(s){
 z<-subset(b,site!=s & !is.na(foreign));100*diff(tapply(z$paired,z$foreign,mean))
})
b$territory<-b$birth_country %in% c("GUM","PRI","VIR","ASM","MNP")
territory<-update(tb,foreign=ifelse(b$territory,0,foreign))
save_table(rbind(data.frame(check="Leave-one-site-out difference range (pp)",estimate=NA,low=min(leave),high=max(leave)),data.frame(check="US territories grouped with USA (pp)",rates(territory,"Alternate coding")[4,c("estimate","low","high")])),"sensitivity")
# Mean BMI difference between the paired and incomplete groups.
profile<-function(des,v,name){
 des$variables$value<-des$variables[[v]]
 dd<-subset(des,is.finite(value));f<-svyglm(value~paired,dd);k<-coef(f)[2];s<-sqrt(vcov(f)[2,2]);q<-qt(.975,degf(des))
 data.frame(cohort=name,n=nrow(dd$variables),paired=sum(dd$variables$paired),incomplete=nrow(dd$variables)-sum(dd$variables$paired),
 mean_incomplete=unname(coef(f)[1]),mean_paired=unname(sum(coef(f))),difference=unname(k),low=unname(k-q*s),high=unname(k+q*s))
}
save_table(rbind(profile(update(tb,bmi=bmi5*5),"bmi","TBESC enrolled"),profile(subset(nd,age>=25 & !(TBQ070 %in% 1)),"BMXBMI","NHANES no documented reaction")),"bmi")
# Case-count audit only; these counts are not time-to-event risks.
save_table(as.data.frame(with(subset(b,!is.na(state)),table(state,d_case_status))),"case_counts")
nd<-update(nd,direction=y,bmi5=BMXBMI/5,agegroup=cut(age,c(25,45,65,Inf),right=FALSE),sex=sex_male,smoke=as.integer(smoking=="Current"))
common<-function(des,name){
 dd<-subset(des,complete.cases(direction,bmi5,agegroup,sex,foreign,smoke))
 fit<-svyglm(direction~bmi5+agegroup+sex+foreign+smoke,dd,family=quasibinomial())
 do.call(rbind,lapply(c("bmi5","foreign"),function(term){
  k<-coef(fit)[term];se<-sqrt(vcov(fit)[term,term]);q<-qt(.975,degf(des))
  data.frame(cohort=name,n=nrow(dd$variables),QFT_only=sum(dd$variables$direction),term=term,OR=exp(k),low=exp(k-q*se),high=exp(k+q*se))
 }))
}
save_table(rbind(common(tb,"TBESC"),common(subset(nd,age>=25),"NHANES")),"common_models")
