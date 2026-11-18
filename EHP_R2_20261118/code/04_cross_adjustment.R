s<-readLines('EHP_R2_20261118/code/01_sensitivity.R')
eval(parse(text=s[seq_len(grep('^adult <-',s)-1)]))
diet_cov<-a[,c('log_cot','fish','shellfish','prior_tb','contact')]
s<-readLines('EHP_R2_20261118/code/02_joint_jackknife.R')
eval(parse(text=s[seq_len(grep('^result\\$independent_beta',s)-1)]))
f<-cbind(f,diet_cov[idx,])
base<-c('age','sex_male','race','educ','PIR','birth','smoking')
extras<-names(diet_cov)
rows<-list()
for(i in seq_len(nrow(reg))){
 rr<-reg[i,]; v<-rr$variable; x<-f[[v]]
 binary<-!is.na(rr$detection_category)&&rr$detection_category=='<50% detected'
 if(binary){lc<-f[[rr$lc_variable]];x<-ifelse(lc%in%c(0,1),as.integer(lc==0),NA)}
 else if(rr$transform=='log2')x<-ifelse(x>0,log2(x),NA)
 cv<-base; add<-extras
 if(isTRUE(rr$smoking_biomarker)){cv<-setdiff(cv,'smoking');add<-setdiff(add,'log_cot')}
 if(isTRUE(rr$urinary))cv<-c(cv,'logucr')
 ww<-weights[[rr$weight]]
 eligible<-ww[,1]>0 & !is.na(f$group) & complete.cases(f[,c(cv,add)]) & is.finite(x)
 for(j in which(result$variable==v)){
  z<-f[[result$metabolite[j]]];if(result$layer[j]=='individual')z<-log2(ifelse(z>0,z,NA))
  keep<-which(eligible&is.finite(z)); wj<-ww[keep,,drop=FALSE]; w0<-wj[,1]
  r<-result[j,c('variable','family','layer','metabolite','exposure_weight','scale')]
  r$n<-length(keep);r$kish<-if(length(keep))sum(w0)^2/sum(w0^2)else NA_real_;r$baseline_status<-result$status[j]
  for(k in c('base','expanded')){
   r$scenario<-k;r$beta<-r$se<-r$low<-r$high<-r$p<-NA_real_;r$status<-'nonestimable'
   if(result$status[j]=='ok' && length(keep)>=50 && length(unique(x[keep]))>=2 && (!binary||min(table(x[keep]))>=5)){
    dat<-f[keep,c(cv,add),drop=FALSE];dat$x<-x[keep]
    X<-model.matrix(reformulate(c('x',cv,if(k=='expanded')add)),dat)
    Y<-(z[keep]-result$log2mean[j])/result$log2sd[j]
    if(qr(X)$rank==ncol(X)){
     coef<-apply(wj,2,function(w)unname(lm.wfit(X,Y,w)$coefficients[2]))
     if(all(is.finite(coef))){
      se<-sqrt(sum(scale*(coef[-1]-coef[1])^2))
      if(se>0){r$status<-'ok';r$beta<-coef[1];r$se<-se;r$low<-coef[1]-qt(.975,17)*se;r$high<-coef[1]+qt(.975,17)*se;r$p<-2*pt(-abs(coef[1]/se),17)}
     }
    }
   }
   rows[[length(rows)+1]]<-r
  }
 }
 if(i%%20==0)cat('Adjusted lipid markers',i,'\n')
}
r<-do.call(rbind,rows)
r$q<-ave(ifelse(r$status=='ok',r$p,1),r$scenario,r$layer,FUN=function(p)p.adjust(p,'BH'))
write.csv(r,file.path(root,'08_joint_adjusted_12784.csv'),row.names=FALSE)
