# Run from the NHANES project root. All exports are aggregate.
options(survey.lonely.psu="adjust")
.libPaths(c("TB_metabolome_exposome_20261115_2130/analysis_MWAS_ExWAS/Rlib",.libPaths()))
library(survey)
library(haven)
out <- "IDP_Paired_Testing_20261121"
write_result <- function(x,n) write.csv(x,file.path(out,"B_results",paste0(n,".csv")),row.names=FALSE,na="")
a <- readRDS(Sys.glob("NHANES_TB_Stage3*/数据/rebuilt_private.rds")[1])$data
raw <- Sys.glob("NHANES离线数据库/*/2011-2012")[1]
for (file in c("Examination/tbx_g.xpt","Laboratory/tb_g.xpt","Demographics/demo_g.xpt")) {
 z <- read_xpt(file.path(raw,file))
 stopifnot(!anyDuplicated(z$SEQN),!anyDuplicated(a$SEQN))
 for (v in intersect(c("TBQ070","LBXTBIN","TBDRUIND","WTMEC2YR","INDFMPIR","DMDBORN4"),names(z))) {
  stopifnot(isTRUE(all.equal(as.numeric(a[[v]]),as.numeric(z[[v]][match(a$SEQN,z$SEQN)]),check.attributes=FALSE)))
 }
}
a <- subset(a,WTMEC2YR>0)
a$complete <- as.integer(a$LBXTBIN %in% 1:2 & is.finite(a$TBDRUIND))
stopifnot(all(a$complete==!is.na(a$group)))
a$eligible <- a$age>=20 & !(a$TBQ070 %in% 1)
a$negative_history <- a$TBQ070 %in% 2
categorize <- function(x) factor(ifelse(is.na(x),"Missing",as.character(x)),levels=c(levels(x),"Missing"))
a$income <- categorize(cut(a$PIR,c(-Inf,1,2,Inf),right=FALSE,labels=c("<1","1 to <2","2 or above")))
a$nativity <- categorize(factor(a$birth,levels=c("US","NonUS"),labels=c("US-born","Foreign-born")))
a$bmi <- categorize(cut(a$BMXBMI,c(-Inf,18.5,25,30,Inf),right=FALSE,labels=c("Underweight","Normal range","Overweight","Obesity")))
a$albumin <- categorize(cut(a$LBXSAL,c(-Inf,3.5,Inf),right=FALSE,labels=c("<3.5","3.5 or above")))
a$vitaminD <- categorize(cut(a$LBXVIDMS,c(-Inf,50,Inf),right=FALSE,labels=c("<50","50 or above")))
validwt <- function(x) ifelse(x %in% c(7777,9999,777,999),NA,x)
loss <- validwt(a$WHD050)-validwt(a$WHD020)
unintent <- ifelse(loss<10,0,ifelse(a$WHQ060 %in% 1:2,as.integer(a$WHQ060==2),NA))
a$risk <- ifelse(is.na(a$BMXBMI)|is.na(unintent),NA,as.integer(a$BMXBMI<18.5|unintent==1))
a$nutritionRisk <- categorize(factor(a$risk,levels=0:1,labels=c("Absent","Present")))
d <- svydesign(~SDMVPSU,strata=~SDMVSTRA,weights=~WTMEC2YR,nest=TRUE,data=a)
ds <- subset(d,eligible)
crit <- qt(.975,17)
stopifnot(degf(ds)==17,nrow(ds$variables)==5235)
rates <- function(des,v=NULL,ref=NULL) {
 des$variables$stratum <- if(is.null(v))factor("All") else droplevels(des$variables[[v]])
 lev <- levels(des$variables$stratum)
 if(is.null(ref))ref <- lev[1]
 rates <- if(length(lev)==1) svymean(~complete,des) else svyby(~complete,~stratum,des,svymean,covmat=TRUE)
 b <- coef(rates); V <- vcov(rates); i <- match(ref,lev)
 do.call(rbind,lapply(seq_along(lev),function(j){
  x <- des$variables$stratum==lev[j]; w<-weights(des)[x]
  p <- b[j]; se<-sqrt(V[j,j]); delta<-b[j]-b[i]; sd<-sqrt(max(0,V[j,j]+V[i,i]-2*V[j,i]))
  n<-sum(x); hits<-sum(des$variables$complete[x])
  interval<-if(hits>0 && hits<n) plogis(qlogis(p)+c(-1,1)*crit*se/(p*(1-p))) else c(NA,NA)
  data.frame(variable=if(is.null(v))"Overall" else v,level=lev[j],reference=ref,n=n,paired=hits,incomplete=n-hits,kish=sum(w)^2/sum(w^2),
   pct=100*p,low=100*interval[1],high=100*interval[2],difference=100*delta,
   diff_low=if(j!=i && n>=5 && hits>0 && hits<n)100*(delta-crit*sd)else NA,
   diff_high=if(j!=i && n>=5 && hits>0 && hits<n)100*(delta+crit*sd)else NA)
 }))
}
domains<-list(all_adult_MEC=subset(d,age>=20),no_known_contraindication=ds,documented_no_prior_reaction=subset(ds,negative_history))
flow<-do.call(rbind,lapply(names(domains),function(n)cbind(domain=n,rates(domains[[n]]))))
write_result(flow,"completion_flow")
refs<-c(income="2 or above",nativity="US-born",bmi="Normal range",albumin="3.5 or above",vitaminD="50 or above",nutritionRisk="Absent")
groups<-do.call(rbind,lapply(names(refs),function(v)rates(ds,v,refs[v])))
write_result(groups,"completion_groups")
strict <- subset(ds,negative_history)
write_result(do.call(rbind,lapply(names(refs),function(v)rates(strict,v,refs[v]))),"completion_groups_documented")

prof <- lapply(c("PIR","BMXBMI","LBXSAL","LBXVIDMS","risk","foreign"),function(v){
 des<-ds; des$variables$foreign<-ifelse(is.na(des$variables$birth),NA,as.integer(des$variables$birth=="NonUS"))
 des$variables$value<-des$variables[[v]]*ifelse(v %in% c("risk","foreign"),100,1)
 des<-subset(des,is.finite(value))
 s<-svyby(~value,~complete,des,svymean,covmat=TRUE)
 b<-coef(s);se<-SE(s);z<-svycontrast(s,c(-1,1))
 data.frame(variable=v,n_incomplete=sum(des$variables$complete==0),n_complete=sum(des$variables$complete==1),incomplete=b[1],inc_low=b[1]-crit*se[1],inc_high=b[1]+crit*se[1],complete=b[2],comp_low=b[2]-crit*se[2],comp_high=b[2]+crit*se[2],difference=coef(z),low=as.numeric(coef(z)-crit*SE(z)),high=as.numeric(coef(z)+crit*SE(z)))
})
write_result(do.call(rbind,prof),"completion_profiles")
dflow<-subset(a,age>=20)
write_result(as.data.frame(with(dflow,table(reaction=ifelse(is.na(TBQ070),"Unknown",ifelse(TBQ070==1,"Yes","No")),qft=ifelse(is.na(LBXTBIN),"Missing",ifelse(LBXTBIN==3,"Indeterminate","Interpretable")),tst=ifelse(is.na(TBDRUIND),"No reading","Reading")))),"completion_disposition")
# Conditional response weighting in the observable eligibility domain. All eligible adults
# remain in the response model. Missing covariates have separate indicator levels.
b<-ds$variables
rf<-glm(complete~age+sex_male+race+nativity+income+bmi,
 data=b,weights=WTMEC2YR/mean(WTMEC2YR),family=quasibinomial())
pr<-predict(rf,type="response");stopifnot(length(pr)==nrow(b),rf$converged)
a$prob<-pr[match(a$SEQN,b$SEQN)]
write_result(data.frame(n=nrow(b),min=min(pr),median=median(pr),max=max(pr),floor=.05,n_below_floor=sum(pr<.05)),"completion_propensity")
reg<-read.csv(file.path(out,"Source_data/01_expanded_registry_188.csv"))
primary<-read.csv(file.path(out,"B_results/primary188.csv"))
stopifnot(identical(primary$variable,reg$variable))
fit_marker<-function(i) {
 r<-reg[i,]; v<-r$variable
 covs<-c("age","sex_male","race","birth","smoking","educ","PIR")
 if(r$smoking_biomarker)covs<-setdiff(covs,"smoking")
 if(r$urinary)covs<-c(covs,"logucr")
 b<-a[,unique(c(covs,"SDMVPSU","SDMVSTRA","age","y","eligible","prob"))]
 x<-a[[v]]
 if(r$detection_category %in% "<50% detected"){
  lc<-a[[r$lc_variable]];x<-ifelse(lc %in% 0:1,as.integer(lc==0),NA)
 }else if(r$transform=="log2") x<-log2(ifelse(x>0,x,NA))
 b$x<-x;b$w<-a[[r$weight]];b<-subset(b,w>0)
 b$domain<-b$eligible & complete.cases(b[,c("y","x",covs)])
 b$w<-b$w/ifelse(b$eligible,pmax(.05,b$prob),1)
 des<-svydesign(~SDMVPSU,strata=~SDMVSTRA,weights=~w,nest=TRUE,data=b)
 design<-subset(des,domain); z<-design$variables
 row<-data.frame(variable=v,primary_status=primary$status[i],primary_OR=primary$OR[i],n=nrow(z),OR=NA,low=NA,high=NA,p=NA,status="invariant")
 if(length(unique(z$x))<2)return(row)
 fit<-try(svyglm(reformulate(c("x",covs),"y"),design,family=quasibinomial()),silent=TRUE)
 if(inherits(fit,"try-error")){row$status<-"fit_error";return(row)}
 k<-coef(fit)["x"];s<-sqrt(vcov(fit)["x","x"])
 row$status<-if(is.finite(k+s)&&fit$converged&&abs(k)<5&&s>0&&s<2)"ok" else "unstable"
 if(row$status=="ok")row[1,c("OR","low","high","p")]<-list(exp(k),exp(k-crit*s),exp(k+crit*s),2*pt(-abs(k/s),17))
 row
}
ipw<-do.call(rbind,lapply(seq_len(nrow(reg)),fit_marker))
ipw$q188<-p.adjust(ifelse(is.finite(ipw$p),ipw$p,1),"BH")
stopifnot(all(ipw$n[is.finite(primary$n)]==primary$n[is.finite(primary$n)]))
write_result(ipw,"completion_ipw188")
capture.output(sessionInfo(),file=file.path(out,"B_results/completion_session.txt"))
