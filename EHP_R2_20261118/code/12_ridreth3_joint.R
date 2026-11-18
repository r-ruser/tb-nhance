root <- 'EHP_PlanB_20261118/B_results'
prior <- 'IJHEH_Unified_Circular_Figures_20261118/Source_data'
a <- readRDS('NHANES_TB_Stage3_Inference_Audit_20261116/数据/rebuilt_private.rds')$data
spec <- readRDS('TB_metabolome_exposome_20261115_2130/临床升级_双组学_20261116/output/analysis_data_private.rds')
for(nm in names(spec$mods)){
 z<-a[,spec$mods[[nm]],drop=FALSE]; a[[paste0('fa_',nm)]]<-log2(rowSums(z))
}
reg<-read.csv(file.path(prior,'01_expanded_registry_188.csv'))
result<-read.csv(file.path(prior,'05_full_spectrum_primary_6392.csv'))
idx<-which(a$age>=20 & a$WTFAS2YR>0)
f<-a[idx,]; f$race<-factor(f$RIDRETH3);stopifnot(all(f$RIDRETH3 %in% c(1,2,3,4,6,7)))
G<-list(cut(f$age,c(20,40,60,Inf),right=FALSE),factor(f$sex_male),factor(f$race))
stopifnot(all(complete.cases(as.data.frame(G))))
strata<-split(seq_len(nrow(f)),f$SDMVSTRA)
psu<-lapply(strata,function(ii)split(ii,f$SDMVPSU[ii]))
stopifnot(all(lengths(psu)>=2))
k<-sum(lengths(psu)); w<-matrix(f$WTFAS2YR,nrow(f),k+1); scale<-numeric(k); col<-1
for(h in psu) for(g in h){
 col<-col+1; num<-length(h); w[unlist(h),col]<-w[unlist(h),col]*num/(num-1);w[g,col]<-0;scale[col-1]<-(num-1)/num
}
# Independent iterative proportional fitting; reference totals are recalculated inside each PSU replicate.
rake <- function(weight,eligible){
 target<-lapply(G,function(g)as.numeric(rowsum(weight,g,reorder=TRUE)))
 ww<-weight*eligible
 for(j in seq_len(300)){
  for(l in seq_along(G)){
   g<-G[[l]]; observed<-as.numeric(rowsum(ww,g,reorder=TRUE)); desired<-target[[l]]
   if(any(observed==0 & desired>0))stop('Empty calibration level')
   fac<-ifelse(observed>0,desired/observed,1);ww<-ww*fac[as.integer(g)]
  }
  err<-max(sapply(seq_along(G),function(l)max(abs(as.numeric(rowsum(ww,G[[l]],reorder=TRUE))-target[[l]])/pmax(1,target[[l]]))))
  if(err<1e-10)return(list(w=ww,error=err))
 }
 stop('IPF failed')
}
weights<-list(WTMEC2YR=w);audit<-list()
for(wt in setdiff(unique(reg$weight),'WTMEC2YR')){
 eligible<-is.finite(f[[wt]])&f[[wt]]>0
 ipf<-lapply(seq_len(ncol(w)),function(j)rake(w[,j],eligible))
 weights[[wt]]<-do.call(cbind,lapply(ipf,'[[','w'))
 audit[[wt]]<-data.frame(component=wt,eligible=sum(eligible),min_positive=min(weights[[wt]][weights[[wt]]>0]),max_weight=max(weights[[wt]]),PSU_replicates=k,max_relative_margin_error=max(sapply(ipf,'[[','error')))
}
write.csv(do.call(rbind,audit),file.path(root,'04_joint_IPF_audit.csv'),row.names=FALSE)
result$independent_beta<-result$jk_se<-result$jk_low<-result$jk_high<-result$jk_p<-NA_real_
result$failed_reps<-NA_integer_
base<-c('age','sex_male','race','educ','PIR','birth','smoking')
for(i in seq_len(nrow(reg))){
 rr<-reg[i,];v<-rr$variable;x<-f[[v]]
 if(!is.na(rr$detection_category)&&rr$detection_category=='<50% detected') {
  lc<-f[[rr$lc_variable]];x<-ifelse(lc%in%c(0,1),as.integer(lc==0),NA)
 } else if(rr$transform=='log2')x<-ifelse(x>0,log2(x),NA)
 cv<-base
 if(isTRUE(rr$smoking_biomarker))cv<-setdiff(cv,'smoking')
 if(isTRUE(rr$urinary))cv<-c(cv,'logucr')
 ww<-weights[[rr$weight]]
 for(j in which(result$variable==v & result$status=='ok')){
  z<-f[[result$metabolite[j]]]
  if(result$layer[j]=='individual')z<-log2(ifelse(z>0,z,NA))
  keep<-which(ww[,1]>0 & !is.na(f$group)&complete.cases(f[,cv])&is.finite(x)&is.finite(z))
  stopifnot(length(keep)==result$n[j])
  dat<-f[keep,cv,drop=FALSE];dat$x<-x[keep]
  X<-model.matrix(reformulate(c('x',cv)),dat)
  Y<-(z[keep]-result$log2mean[j])/result$log2sd[j]
  slope<-function(wt) {
   # Drop absent covariate levels in individual replicates; an aliased exposure remains unavailable.
   fit<-lm.wfit(X,Y,wt,tol=1e-9)
   unname(fit$coefficients[2])
  }
  bet<-apply(ww[keep,,drop=FALSE],2,slope)
  result$independent_beta[j]<-bet[1]
  result$failed_reps[j]<-sum(!is.finite(bet[-1]))
  if(any(!is.finite(bet)))next
  ss<-sqrt(sum(scale*(bet[-1]-bet[1])^2))
  if(ss<=0)next
  result[j,c('jk_se','jk_low','jk_high','jk_p')]<-list(ss,bet[1]-qt(.975,17)*ss,bet[1]+qt(.975,17)*ss,2*pt(-abs(bet[1]/ss),17))
 }
 if(i%%20==0)cat('Markers checked',i,'\n')
}
result$jk_q<-ave(ifelse(is.finite(result$jk_p),result$jk_p,1),result$layer,FUN=function(z)p.adjust(z,'BH'))
write.csv(result,file.path(root,'05_joint_jackknife_6392.csv'),row.names=FALSE)
stopifnot(nrow(result)==6392,all(vapply(audit,function(z)z$max_relative_margin_error<1e-9,TRUE)))
capture.output(sessionInfo(),file=file.path(root,'R_sessionInfo_joint.txt'))
