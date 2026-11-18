script<-readLines('EHP_R2_20261118/code/01_sensitivity.R')
pos<-grep('^results<-list',script)
# Reuse source linkage, outcome definitions and analysis domains.
eval(parse(text=script[seq_len(pos-1)]))
library(splines)
results<-curves<-list()
display<-c('LBXBCD','LBXBPB','LBXTHG','LBXBSE','LBXSIR','LBXSSE','LBXSZN','LBXSCU')
for(i in seq_len(nrow(reg))){
 rr<-reg[i,]; v<-rr$variable
 cv<-base
 if(isTRUE(rr$smoking_biomarker))cv<-setdiff(cv,'smoking')
 if(isTRUE(rr$urinary))cv<-c(cv,'logucr')
 cv<-c(cv,setdiff(extra,if(isTRUE(rr$smoking_biomarker))'log_cot'else ''))
 x<-as.numeric(a[[v]])
 if(rr$transform=='log2')x<-ifelse(x>0,log2(x),NA)
 binary <- !is.na(rr$detection_category)&&rr$detection_category=='<50% detected'
 ok<-which(is.finite(a[[rr$weight]])&a[[rr$weight]]>0)
 for(k in c('coded_T10','QFT_antigen','TST_induration')){
  y<-targets[[k]]
  d<-a[ok,c('SDMVSTRA','SDMVPSU',cv),drop=FALSE]
  d$x<-x[ok];d$y<-y[ok];d$wt<-a[[rr$weight]][ok]
  d$keep<-paired[ok]&complete.cases(d)&is.finite(d$x)&is.finite(d$y)
  j<-which(d$keep); w<-d$wt[j]; n<-length(j)
  z<-data.frame(variable=v,outcome=k,n=n,events=if(k=='coded_T10')sum(d$y[j])else NA,kish=if(n)sum(w)^2/sum(w^2)else NA,status='nonestimable',knots=NA_character_,p_overall=NA,p_nonlinear=NA)
  if(binary){z$status<-'binary_detection';results[[length(results)+1]]<-z;next}
  knots<-if(n)as.numeric(quantile(d$x[j],c(.05,.35,.65,.95)))else numeric()
  if(n>=50 && length(unique(knots))==4 && (k!='coded_T10'||min(z$events,n-z$events)>=20)){
   basis<-function(x)ns(x,knots=knots[2:3],Boundary.knots=knots[c(1,4)])
   # Reparameterize the natural cubic basis so the first term is explicitly linear.
   X<-basis(d$x[j]); coeff<-qr.solve(X-X[rep(1,n),],d$x[j]-d$x[j][1])
   drop<-which.max(abs(coeff)); cols<-setdiff(1:3,drop)
   d$s1<-d$s2<-NA_real_; d$s1[j]<-X[,cols[1]];d$s2[j]<-X[,cols[2]]
   des<-subset(svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~wt,nest=TRUE,data=d),keep)
   f<-tryCatch(suppressWarnings(svyglm(reformulate(c('x','s1','s2',cv),'y'),des,family=if(k=='coded_T10')quasibinomial()else gaussian())),error=function(e)NULL)
   if(!is.null(f)&&f$converged&&degf(des)>=10&&f$rank==length(coef(f))&&all(c('x','s1','s2')%in%names(coef(f)))){
    tt<-c('x','s1','s2');bb<-coef(f)[tt]; V<-vcov(f)[tt,tt]
    wald<-function(ii)tryCatch(pf(drop(crossprod(bb[ii],solve(V[ii,ii],bb[ii])))/length(ii),length(ii),17,lower.tail=FALSE),error=function(e)NA_real_)
    z$p_overall<-wald(1:3);z$p_nonlinear<-wald(2:3);z$knots<-paste(signif(knots,8),collapse='|')
    if(all(is.finite(c(z$p_overall,z$p_nonlinear)))&&all(diag(V)>0)) {
     z$status<-'ok'
     if(v%in%display){
      ref<-median(d$x[j]); g<-seq(knots[1],knots[4],length.out=61)
      M<-cbind(g-ref,sweep(basis(g)[,cols,drop=FALSE],2,basis(ref)[1,cols],'-'))
      est<-drop(M%*%bb);se<-sqrt(pmax(0,rowSums((M%*%V)*M)))
      curves[[length(curves)+1]]<-data.frame(variable=v,outcome=k,x=g,reference=ref,beta=est,se,low=est-qt(.975,17)*se,high=est+qt(.975,17)*se,n)
     }
    }
   }
  }
  results[[length(results)+1]]<-z
 }
}
r<-do.call(rbind,results)
for(v in c('overall','nonlinear'))r[[paste0('q_',v)]]<-ave(ifelse(r$status=='ok',r[[paste0('p_',v)]],1),r$outcome,FUN=function(p)p.adjust(p,'BH'))
write.csv(r,file.path(root,'06_RCS_tests_564.csv'),row.names=FALSE)
write.csv(do.call(rbind,curves),file.path(root,'07_RCS_curve_coordinates.csv'),row.names=FALSE)
