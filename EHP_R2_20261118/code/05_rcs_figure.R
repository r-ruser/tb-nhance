.libPaths(c('TB_metabolome_exposome_20261115_2130/analysis_MWAS_ExWAS/Rlib',.libPaths()))
library(ggplot2)
root<-'IJHEH_Validation_20261118'
d<-read.csv(file.path(root,'Validation/07_RCS_curve_coordinates.csv'))
cols<-c(LBXBCD='Blood cadmium',LBXBPB='Blood lead',LBXTHG='Total mercury',LBXBSE='Blood selenium',LBXSIR='Serum iron',LBXSSE='Serum selenium',LBXSZN='Serum zinc',LBXSCU='Serum copper')
units<-c(LBXBCD='µg/L',LBXBPB='µg/dL',LBXTHG='µg/L',LBXBSE='µg/L',LBXSIR='µg/dL',LBXSSE='µg/L',LBXSZN='µg/dL',LBXSCU='µg/dL')
for(k in c('coded_T10','QFT_antigen','TST_induration')){
 z<-subset(d,outcome==k)
 labels<-paste0(cols[z$variable], ' (', units[z$variable], ')\nn=', z$n)
 z$marker<-factor(labels,levels=unique(labels[order(match(z$variable,names(cols)))]))
 if(k=='coded_T10'){z$beta<-exp(z$beta);z$low<-exp(z$low);z$high<-exp(z$high)}
 ref<-if(k=='coded_T10')1 else 0
 yl<-switch(k,coded_T10='Discordance odds ratio',QFT_antigen='Difference in transformed QFT response',TST_induration='Difference in log(1 + induration)')
 p<-ggplot(z,aes(x,beta))+geom_hline(yintercept=ref,colour='#707B85',linetype=2,linewidth=.35)+
 geom_ribbon(aes(ymin=low,ymax=high),fill='#0072B2',alpha=.15)+geom_line(colour='#0072B2',linewidth=.55)+
 facet_wrap(~marker,ncol=4,scales='free_x')+labs(x='Exposure concentration (log2 original units)',y=yl)+
 theme_classic(base_family='Arial',base_size=8)+theme(strip.background=element_blank(),strip.text=element_text(face='bold'))
 if(k=='coded_T10')p<-p+scale_y_log10(labels=scales::label_number())
 stem<-file.path(root,'Figures',paste0('FigureS2_',k))
 ggsave(paste0(stem,'.png'),p,width=7.2,height=4.5,dpi=600,bg='white')
 ggsave(paste0(stem,'.svg'),p,width=7.2,height=4.5,bg='white')
 ggsave(paste0(stem,'.pdf'),p,width=7.2,height=4.5,bg='white',device=function(filename,width,height,...)quartz(type='pdf',file=filename,width=width,height=height,family='Arial'))
}
