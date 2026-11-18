.libPaths(c('TB_metabolome_exposome_20261115_2130/analysis_MWAS_ExWAS/Rlib',.libPaths()))
library(survey);library(ggplot2)
options(survey.lonely.psu='adjust')
out<-'EHP_R2_20261118';o<-file.path(out,'B_results');fig<-file.path(out,'Figures')
saveout<-function(x,n)write.csv(x,file.path(o,paste0(n,'.csv')),row.names=F,na='')
saveplot<-function(g,n,w,h){for(ext in c('svg','pdf','png'))ggsave(file.path(fig,paste0(n,'.',ext)),g,width=w,height=h,dpi=350,bg='white',device=if(ext=='pdf')function(filename,width,height,...)quartz(type='pdf',file=filename,width=width,height=height)else ext)}
s<-read.csv(file.path(o,'primary188.csv'))
lip<-read.csv(file.path(out,'Source_data/Figure3_all_34_discordance.csv'))
z<-read.csv(file.path(out,'Supplementary_Data_S1_full_cross_omic_6392.csv'));stopifnot(nrow(z)==6392,!anyDuplicated(z[c('variable','metabolite')]))
i<-match(z$variable,s$variable);j<-match(z$metabolite,lip$metabolite)
z$direct<-s$b[i];z$lipid_outcome<-lip$beta[j];z$direct_status<-s$status[i];z$lipid_status<-lip$status[j]
z$direction<-ifelse(z$status!='ok'|z$direct_status!='ok'|z$lipid_status!='ok'|!is.finite(z$beta*z$direct*z$lipid_outcome),'unavailable',ifelse(z$beta*z$direct*z$lipid_outcome==0,'zero',ifelse(z$beta*z$direct*z$lipid_outcome>0,'aligned','opposed')))
saveout(z,'bridge6392');saveout(subset(z,variable%in%c('LBXBCD','LBX2DF','LBXVFN')),'profiles102')
a1<-as.data.frame.matrix(table(interaction(z$variable,z$layer,drop=T),factor(z$direction,levels=c('aligned','opposed','zero','unavailable'))))
a1$variable<-sub('[.].*$','',rownames(a1));a1$layer<-sub('^[^.]*[.]','',rownames(a1));a1$family<-s$family[match(a1$variable,s$variable)];a1$fraction<-with(a1,aligned/(aligned+opposed));saveout(a1,'bridge_by_marker')
g<-ggplot(subset(a1,is.finite(fraction)),aes(fraction,family))+geom_vline(xintercept=.5,lty=2,col='grey65')+geom_point(position=position_jitter(width=0,height=.12,seed=409),alpha=.55,size=1.1,col='#21667A')+facet_wrap(~layer,nrow=1)+labs(x='Fraction of evaluable nonzero sign comparisons aligned',y=NULL)+theme_bw(base_size=8)
saveplot(g,'FigureS3_direction',7,4.5)
f<-subset(lip,layer=='module');crit<-qt(.975,17);delta<-uniroot(function(x)pt(-crit,17,ncp=x)+1-pt(crit,17,ncp=x)-.8,c(0,10))$root
f$detect_upper<-exp(delta*f$se);f$detect_lower<-1/f$detect_upper;saveout(f,'module_precision')
a<-readRDS('NHANES_TB_Stage3_Inference_Audit_20261116/数据/rebuilt_private.rds')$data
a$race6<-factor(a$RIDRETH3,levels=c(1,2,3,4,6,7),labels=c('Mexican American','Other Hispanic','NH White','NH Black','NH Asian','Other/multiracial'))
a$duration<-factor(ifelse(a$birth=='NonUS'&a$DMDYRSUS%in%1:9,ifelse(a$DMDYRSUS<4,'<10 y',ifelse(a$DMDYRSUS<6,'10-<20 y','20+ y')),NA))
nut<-read.csv(file.path(out,'Source_data/18_nutrition_four_state_descriptive.csv'));nut$stat<-'mean (95% CI)';nut$q25<-nut$q75<-NA_real_
nut[nut$variable=='risk',c('value','low','high')]<-100*nut[nut$variable=='risk',c('value','low','high')]
dem<-list()
for(v in c('LBXSCU','LBXSSE','LBXSZN','LBXSIR','URXUIO','race6','educ','PIR','duration')){
 w<-if(v%in%c('LBXSCU','LBXSSE','LBXSZN','URXUIO'))'WTSA2YR'else 'WTMEC2YR'
 b<-data.frame(a[c('age','group','SDMVPSU','SDMVSTRA')],val=a[[v]],wt=a[[w]])
 b<-subset(b,wt>0);d<-svydesign(~SDMVPSU,strata=~SDMVSTRA,weights=~wt,data=b,nest=T)
 for(k in c('Q0T0','Q1T0','Q0T1','Q1T1')){
   ds<-subset(d,age>=20&group==k&!is.na(val));n<-nrow(ds$variables)
   if(v%in%nut$variable){
    idx<-which(nut$variable==v&nut$group==k);stopifnot(n==nut$n[idx])
    if(v=='URXUIO'){
     q<-as.numeric(coef(svyquantile(~val,ds,c(.25,.5,.75),ci=F))); ci<-confint(svyquantile(~val,ds,.5,ci=T,df=17))
     nut[idx,c('value','low','high')]<-c(q[2],ci);nut$stat[idx]<-'median (95% CI)';nut$q25[idx]<-q[1];nut$q75[idx]<-q[3]
    }else{
     m<-svymean(~log(val),ds);nut[idx,c('value','low','high')]<-exp(c(coef(m),confint(m,df=17)));nut$stat[idx]<-'geometric mean (95% CI)'
    }
   }else{
    m<-svymean(~val,ds);p<-if(is.factor(b$val))100 else 1;ci<-confint(m,df=17)
    dem[[paste(v,k)]]<-data.frame(variable=v,category=sub('^val','',names(coef(m))),group=k,n=n,value=p*coef(m),low=p*ci[,1],high=p*ci[,2])
   }
 }
}
saveout(nut,'nutrition');saveout(do.call(rbind,dem),'demographics')
nut$state<-factor(nut$group,levels=c('Q0T0','Q1T0','Q0T1','Q1T1'),labels=c('Q−T−','Q+T−','Q−T+','Q+T+'))
levels<-unique(nut$variable);nut$facet<-factor(nut$variable,levels=levels,labels=paste0(nut$label[match(levels,nut$variable)],' — ',ifelse(levels=='URXUIO','median',ifelse(grepl('^LBX',levels),'geometric mean',ifelse(levels=='risk','prevalence','mean')))))
g<-ggplot(nut,aes(value,state))+geom_segment(aes(x=low,xend=high,yend=state),col='#4C7884')+geom_point(size=1.3,col='#16465D')+facet_wrap(~facet,scales='free_x',ncol=2,labeller=label_wrap_gen(26))+labs(x=NULL,y=NULL)+theme_bw(base_size=7)+theme(strip.text=element_text(size=7),panel.grid.minor=element_blank())
saveplot(g,'Figure1_C_nutrition',5.1,5.8)
