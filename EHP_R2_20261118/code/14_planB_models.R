options(survey.lonely.psu="adjust")
.libPaths(c("TB_metabolome_exposome_20261115_2130/analysis_MWAS_ExWAS/next_stage_validation_20261116/Rlib","TB_metabolome_exposome_20261115_2130/analysis_MWAS_ExWAS/Rlib",.libPaths()))
library(survey);library(svyVGAM)
out<-"EHP_PlanB_20261118";o<-file.path(out,"B_results")
saveout<-function(z,n)write.csv(z,file.path(o,paste0(n,".csv")),row.names=FALSE,na="")
a<-readRDS("NHANES_TB_Stage3_Inference_Audit_20261116/数据/rebuilt_private.rds")$data
reg<-read.csv(file.path(out,"Source_data/01_expanded_registry_188.csv"))
s<-read.csv(file.path(out,"Supplementary_Data_S3_outcome_screen_188.csv"));stopifnot(nrow(s)==188,!anyDuplicated(s$variable),identical(s$variable,reg$variable))
for(n in c("inference_df","residual_df")){
 df<-s[[n]];ok<-which(s$status=="ok"&is.finite(df)&df>0)
 p<-rep(NA_real_,188);p[ok]<-2*pt(-abs(s$b[ok]/s$se[ok]),df[ok])
 s[[paste0("p_",n)]]<-p;s[[paste0("q_",n)]]<-p.adjust(ifelse(is.finite(p),p,1),"BH")
}
stopifnot(max(abs(s$q_188-s$q_inference_df),na.rm=TRUE)<1e-10)
s$unit_comparison<-ifelse(reg$detection_category %in% "<50% detected","Detected vs below LOD",ifelse(reg$transform=="log2","Per doubling","Per original unit"))
s$p_normal<-ifelse(s$status=="ok",2*pnorm(-abs(s$b/s$se)),NA_real_);s$q_normal<-p.adjust(ifelse(is.finite(s$p_normal),s$p_normal,1),"BH");s$q_188<-s$q_inference_df;s$q_183<-NULL;saveout(s,"primary188")
a$race6<-factor(a$RIDRETH3,levels=c(3,1,2,4,6,7),labels=c("NH White","Mexican American","Other Hispanic","NH Black","NH Asian","Other or multiracial"))
a$duration<-factor(ifelse(a$birth=="NonUS"&a$DMDYRSUS%in%1:9,ifelse(a$DMDYRSUS<=3,"<10 y",ifelse(a$DMDYRSUS<=5,"10-<20 y","20+ y")),NA))
a$g<-factor(a$group,levels=c("Q0T0","Q1T0","Q0T1","Q1T1"))
cv<-c("age","sex_male","race","birth","smoking","educ","PIR")
tabout<-list();smoke<-list();cells<-list();strat<-list();race<-list();multi<-list();diagn<-list()
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
 fitone<-function(ds,co,term="x"){
  n<-nrow(ds$variables);e<-sum(ds$variables$y);w<-weights(ds)
  ans<-data.frame(variable=v,n=n,events=e,kish=if(n)sum(w)^2/sum(w^2)else NA,OR=NA,low=NA,high=NA,p=NA,status="insufficient")
  if(min(e,n-e)<5||length(unique(ds$variables$x))<2)return(ans)
  f<-try(suppressWarnings(svyglm(reformulate(c("x",co),"y"),ds,family=quasibinomial())),silent=T)
  if(inherits(f,"try-error"))return(ans)
  if(!term%in%names(coef(f)))return(ans)
  B<-coef(f)[term];S<-sqrt(vcov(f)[term,term])
  if(!is.finite(B+S)||abs(B)>5||S>2||S<=0||!f$converged||f$rank<length(coef(f))){ans$status<-"unstable";return(ans)}
  ans[1,c("OR","low","high","p","status")]<-list(exp(B),exp(B-qt(.975,17)*S),exp(B+qt(.975,17)*S),2*pt(-abs(B/S),17),"ok");ans
 }
 race[[v]]<-fitone(subset(dd,!is.na(race6)),c(setdiff(C,"race"),"race6"))
 if(v%in%c("LBX2DF","LBXVFN")){
  z<-dd$variables
  for(g in c("All","US","NonUS")){
   sub<-if(g=="All") z else z[z$birth==g,,drop=FALSE]
   tb<-as.data.frame(table(detected=factor(sub$x,0:1),QFT_only=factor(sub$y,0:1)));tabout[[paste(v,g)]]<-cbind(variable=v,birth=g,tb)
  }
  smoke[[v]]<-cbind(variable=v,as.data.frame(table(detected=factor(z$x,0:1),smoking=z$smoking)))
  cells[[v]]<-do.call(rbind,lapply(split(z,z$x),function(h)data.frame(variable=v,detected=h$x[1],n=nrow(h),events=sum(h$y),weighted_event_pct=100*weighted.mean(h$y,h$wt),kish=sum(h$wt)^2/sum(h$wt^2))))
  for(adj in c(TRUE,FALSE))diagn[[paste(v,adj)]]<-cbind(model=if(adj)"smoking adjusted" else "smoking omitted",fitone(dd,if(adj)C else setdiff(C,"smoking")))
  # Unweighted Firth diagnostic via the Jeffreys penalized binomial likelihood.
  X<-model.matrix(reformulate(c("x",C)),z);y<-z$y
  loss<-function(beta){eta<-drop(X%*%beta);p<-plogis(eta);I<-crossprod(X,X*(p*(1-p)));detI<-determinant(I,logarithm=TRUE);if(detI$sign<=0)return(1e30)
   -sum(y*eta-pmax(eta,0)-log1p(exp(-abs(eta))))-.5*as.numeric(detI$modulus)}
  f<-optim(rep(0,ncol(X)),loss,method="BFGS",hessian=T,control=list(maxit=1200,reltol=1e-12));j<-match("x",colnames(X));se<-sqrt(solve(f$hessian)[j,j]);bb<-f$par[j]
  diagn[[paste(v,"firth")]]<-data.frame(model="unweighted Firth, Wald CI",variable=v,n=nrow(z),events=sum(y),kish=NA,OR=exp(bb),low=exp(bb-1.96*se),high=exp(bb+1.96*se),p=NA,status=if(f$convergence==0)"ok" else "failed")
 }
 if(v=="LBXBCD"){
  co<-c(setdiff(C,c("birth","race")),"race6")
  for(country in c("US","NonUS"))strat[[country]]<-cbind(model=country,fitone(subset(dd,birth==country),co))
  strat[["duration"]]<-cbind(model="foreign born, duration adjusted",fitone(subset(dd,birth=="NonUS"&!is.na(duration)),c(co,"duration")))
  strat[["interaction"]]<-cbind(model="birth interaction ratio of ORs",fitone(dd,c(co,"birth","x:birth"),"x:birthNonUS"))
 }
 if(r$weight=="WTMEC2YR"){
  f<-try(suppressWarnings(svy_vglm(reformulate(c("x",C),"g"),family=VGAM::multinomial(refLevel="Q0T0"),design=cc)),silent=T)
  Ls<-list(Q1T0_vs_Q0T0=c(1,0,0),Q0T1_vs_Q0T0=c(0,1,0),Q1T1_vs_Q0T0=c(0,0,1),Q1T0_vs_Q0T1=c(1,-1,0))
  for(k in names(Ls)){
   row<-data.frame(variable=v,contrast=k,n=nrow(cc$variables),beta=NA,se=NA,OR=NA,low=NA,high=NA,status="failed")
   if(!inherits(f,"try-error")){
    idx<-match(paste0("x:",1:3),names(coef(f)));L<-Ls[[k]];bb<-sum(coef(f)[idx]*L);ss<-sqrt(drop(L%*%vcov(f)[idx,idx]%*%L))
    if(is.finite(bb+ss)&&ss>0) row[1,4:9]<-list(bb,ss,exp(bb),exp(bb-qt(.975,17)*ss),exp(bb+qt(.975,17)*ss),if(abs(bb)>5||ss>2)"unstable"else "ok")
   };multi[[paste(v,k)]]<-row
  }
 }
}
zz<-do.call(rbind,race);zz$q188<-p.adjust(ifelse(zz$status=="ok",zz$p,1),"BH");saveout(zz,"race6_188")
for(n in c("tabout","smoke","cells","strat","multi","diagn"))saveout(do.call(rbind,get(n)),n)
old<-read.csv("NHANES_TB_Stage3_Inference_Audit_20261116/数据/10_多项模型.csv");saveout(merge(subset(do.call(rbind,multi),variable=="LBXBCD"),old,by="contrast",suffixes=c("_package","_custom")),"multinomial_check")
capture.output(sessionInfo(),file=file.path(o,"sessionInfo.txt"))
