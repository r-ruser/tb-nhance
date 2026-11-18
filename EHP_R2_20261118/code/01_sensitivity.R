options(survey.lonely.psu='adjust')
root <- 'IJHEH_Validation_20261118/Validation'
prior <- 'IJHEH_Unified_Circular_Figures_20261118/Source_data'
.libPaths(c('TB_metabolome_exposome_20261115_2130/analysis_MWAS_ExWAS/Rlib',.libPaths()))
library(survey)
library(haven)
a <- readRDS('NHANES_TB_Stage3_Inference_Audit_20261116/数据/rebuilt_private.rds')$data
reg <- read.csv(file.path(prior,'01_expanded_registry_188.csv'))
diet <- read_xpt('NHANES离线数据库/nhanes1999-2020数据合集/2011-2012/Dietary/dr1tot_g.xpt')
stopifnot(!anyDuplicated(diet$SEQN))
ix <- match(a$SEQN,diet$SEQN)
for(v in c('DRD340','DRD360','DR1DRSTZ','DR1TKCAL','DR1TPROT','DR1TPFAT','DR1TFIBE','WTDRD1')) a[[v]] <- as.numeric(diet[[v]][ix])
yes <- function(z) ifelse(z %in% c(1,2),as.integer(z==1),NA)
a$log_cot <- ifelse(a$LBXCOT>0,log2(a$LBXCOT),NA)
a$fish <- yes(a$DRD360); a$shellfish <- yes(a$DRD340)
a$prior_tb <- yes(a$TBQ040); a$contact <- yes(a$TBQ060)
a$kcal <- ifelse(a$DR1DRSTZ==1 & a$DR1TKCAL>0,a$DR1TKCAL/1000,NA)
a$protein <- a$DR1TPROT/a$kcal; a$pufa <- a$DR1TPFAT/a$kcal; a$fiber <- a$DR1TFIBE/a$kcal
adult <- with(a,age>=20 & WTMEC2YR>0 & is.finite(WTMEC2YR))
paired <- adult & !is.na(a$group)
stopifnot(sum(paired)==4215,sum(paired & !is.na(a$y))==398)
for(v in c('log_cot','fish','shellfish','prior_tb','contact','kcal','protein','pufa','fiber','LBXTBA','LBXTBN','LBXTBM')){
 for(g in c('Q0T0','Q1T0','Q0T1','Q1T1')){
  j<-which(paired & a$group==g); z<-a[[v]][j]
  cat(v,g,length(j),sum(is.finite(z)),sum(z==1,na.rm=TRUE),'\n')
 }
}
delta <- a$LBXTBA-a$LBXTBN
mit <- a$LBXTBM-a$LBXTBN
discord <- function(q,t) ifelse(q!=t,as.numeric(q),NA_real_)
rules <- list()
for(t in c(5,10,15)) for(q in c(.20,.35,.70)){
 qp <- delta>=q-1e-10 & delta>=.25*a$LBXTBN-1e-10
 qp[!paired | is.na(mit) | is.na(a$LBXTBN) | a$LBXTBN>8 | (!(qp %in% TRUE) & mit<.5)] <- NA
 rules[[sprintf('Q%.2f_T%d',q,t)]] <- discord(qp,a$TBDRUIND>=t)
}
rules[['coded_T5']]<-discord(a$LBXTBIN==1,a$TBDRUIND>=5)
rules[['coded_T15']]<-discord(a$LBXTBIN==1,a$TBDRUIND>=15)
rules[['away_from_cutpoints']]<-ifelse(abs(delta-.35)>.15 & abs(a$TBDRUIND-10)>2,a$y,NA)
cts <- do.call(rbind,lapply(names(rules),function(k) data.frame(scenario=k,qft_only=sum(paired & rules[[k]]==1,na.rm=TRUE),tst_only=sum(paired & rules[[k]]==0,na.rm=TRUE))))
write.csv(cts,file.path(root,'01_threshold_counts.csv'),row.names=FALSE)
qc <- paired & is.finite(delta) & is.finite(mit) & a$LBXTBN<=8
pred <- ifelse(delta>=.35-1e-10 & delta>=.25*a$LBXTBN-1e-10,1,ifelse(mit>=.5,2,3))
write.csv(as.data.frame(table(coded=a$LBXTBIN[qc],reconstructed=pred[qc])),file.path(root,'02_classifier_audit.csv'),row.names=FALSE)
extra <- c('log_cot','fish','shellfish','prior_tb','contact')
base <- c('age','sex_male','race','educ','PIR','birth','smoking')
nutrition <- list(BMI=a$BMXBMI,albumin=a$LBXSAL,vitamin_D=a$LBXVIDMS)
targets <- c(rules,list(QFT_antigen=asinh(delta/.35),TST_induration=log1p(a$TBDRUIND),QFT_mitogen=asinh(mit/.5)))
targets[['TST_any']] <- as.numeric(a$TBDRUIND>0)
targets[['coded_T10']] <- a$y
fit <- function(x,y,cv,wt,scenario,variable,scale,linear=FALSE,domain=rep(TRUE,nrow(a))){
 ok <- which(is.finite(wt)&wt>0&!is.na(a$SDMVSTRA)&!is.na(a$SDMVPSU))
 d <- a[ok,c('SDMVSTRA','SDMVPSU',cv),drop=FALSE];d$x<-x[ok];d$y<-y[ok];d$wt<-wt[ok]
 d$keep <- paired[ok]&domain[ok]&complete.cases(d)&is.finite(d$x)&is.finite(d$y)
 j<-which(d$keep); m<-length(j); ne<-if(linear) NA_integer_ else sum(d$y[j]==1)
 w<-d$wt[j]
 cellmin<-if(scale=='Detection' && !linear && m) min(table(factor(d$x[j],0:1),factor(d$y[j],0:1))) else NA_integer_
 r<-data.frame(scenario,variable,scale,linear,n=m,events=ne,minor_cell=cellmin,kish=if(m)sum(w)^2/sum(w^2)else NA,design_df=NA,beta=NA,se=NA,low=NA,high=NA,p=NA,status='nonestimable')
 if(m<30 || length(unique(d$x[j]))<2 || (!linear && min(ne,m-ne)<10) || (!is.na(cellmin)&&cellmin<5) || (scale=='Detection' && min(table(d$x[j]))<5)) return(r)
 design<-svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~wt,nest=TRUE,data=d)
 sub<-subset(design,keep)
 r$design_df<-degf(sub)
 f<-tryCatch(suppressWarnings(svyglm(reformulate(c('x',cv),'y'),sub,family=if(linear) gaussian() else quasibinomial())),error=function(e)NULL)
 if(is.null(f)||!f$converged||!('x'%in%names(coef(f)))||f$rank<length(coef(f))||r$design_df<10)return(r)
 bb<-unname(coef(f)['x']);ss<-sqrt(vcov(f)['x','x'])
 if(!is.finite(bb)||!is.finite(ss)||ss<=0||(!linear&&(abs(bb)>5||ss>2))){r$status<-'unstable';return(r)}
 r[c('beta','se','low','high','p','status')]<-list(bb,ss,bb-qt(.975,17)*ss,bb+qt(.975,17)*ss,2*pt(-abs(bb/ss),17),'ok')
 r
}
results<-list()
add<-function(...)results[[length(results)+1L]]<<-fit(...)
for(i in seq_len(nrow(reg))){
 rr<-reg[i,]; v<-rr$variable;x<-as.numeric(a[[v]]);scale<-'Original'
 if(!is.na(rr$detection_category)&&rr$detection_category=='<50% detected'){
 lc<-a[[rr$lc_variable]];x<-ifelse(lc%in%c(0,1),as.integer(lc==0),NA);scale<-'Detection'
 } else if(rr$transform=='log2'){x<-ifelse(x>0,log2(x),NA);scale<-'Doubling'}
 cv<-base
 if(isTRUE(rr$smoking_biomarker))cv<-setdiff(cv,'smoking')
 if(isTRUE(rr$urinary))cv<-c(cv,'logucr')
 w<-a[[rr$weight]]
 for(k in names(targets)) add(x,targets[[k]],cv,w,k,v,scale,k%in%c('QFT_antigen','QFT_mitogen','TST_induration'))
 expanded<-c(cv,setdiff(extra,if(isTRUE(rr$smoking_biomarker))'log_cot'else ''))
 keep<-complete.cases(a[,expanded])
 for(yname in c('coded_T10','QFT_antigen','TST_induration')){
  add(x,targets[[yname]],cv,w,paste0(yname,'_matched_base'),v,scale,yname!='coded_T10',keep)
  add(x,targets[[yname]],expanded,w,paste0(yname,'_expanded'),v,scale,yname!='coded_T10',keep)
 }
 for(block in c('cotinine','seafood','TB_history')){
  cvb<-switch(block,cotinine=if(isTRUE(rr$smoking_biomarker))character()else 'log_cot',seafood=c('fish','shellfish'),TB_history=c('prior_tb','contact'))
  if(length(cvb)) {
   block_keep<-complete.cases(a[,cvb,drop=FALSE])
   add(x,a$y,cv,w,paste0(block,'_matched_base'),v,scale,FALSE,block_keep)
   add(x,a$y,c(cv,cvb),w,paste0(block,'_adjusted'),v,scale,FALSE,block_keep)
  }
 }
 if(rr$weight=='WTMEC2YR'){
  cd<-c(expanded,'kcal','protein','pufa','fiber');keepd<-complete.cases(a[,cd])
  add(x,a$y,expanded,a$WTDRD1,'recall_matched_base',v,scale,FALSE,keepd)
  add(x,a$y,cd,a$WTDRD1,'recall_adjusted',v,scale,FALSE,keepd)
 }
 if(i%%25==0) cat('Markers processed',i,'\n')
}
for(v in names(nutrition)){
 x<-nutrition[[v]];x<-(x-mean(x[paired],na.rm=TRUE))/sd(x[paired],na.rm=TRUE)
 cv<-c(base,extra); keep<-complete.cases(a[,cv])
 for(k in c('coded_T10','QFT_antigen','TST_induration')){
  add(x,targets[[k]],base,a$WTMEC2YR,paste0('nutrition_',k,'_base'),v,'Per SD',k!='coded_T10',keep)
  add(x,targets[[k]],cv,a$WTMEC2YR,paste0('nutrition_',k,'_expanded'),v,'Per SD',k!='coded_T10',keep)
 }
}
out<-do.call(rbind,results)
out$q<-ave(ifelse(out$status=='ok',out$p,1),out$scenario,FUN=function(p) p.adjust(p,n=if(length(p)>20)188 else length(p),method='BH'))
write.csv(out,file.path(root,'03_outcome_sensitivity.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(root,'R_sessionInfo_outcome.txt'))
