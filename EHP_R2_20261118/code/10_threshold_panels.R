# Run from the new manuscript directory with its sibling NHANES audit directory available.
library(grid)
library(svglite)
a <- readRDS("../NHANES_TB_Stage3_Inference_Audit_20261116/数据/rebuilt_private.rds")$data
a <- subset(a, age >= 20 & WTMEC2YR > 0 & !is.na(group))
stopifnot(nrow(a)==4215, !anyDuplicated(a$SEQN))
m <- read.csv("Validation/03_outcome_sensitivity.csv")
reference <- subset(m, scenario=="coded_T10" & status=="ok",c(variable,beta))
counts <- head(read.csv("Validation/01_threshold_counts.csv"),9)
rows <- list()
for (i in seq_len(nrow(counts))) {
 k <- counts$scenario[i]
 q <- as.numeric(substr(k,2,5))
 t <- as.numeric(sub(".*T","",k))
 d <- a$LBXTBA-a$LBXTBN
 mit <- a$LBXTBM-a$LBXTBN
 positive <- d>=q-1e-10 & d>=0.25*a$LBXTBN-1e-10
 eligible <- is.finite(d) & is.finite(mit) & a$LBXTBN<=8 & (positive | mit>=0.5)
 state <- ifelse(eligible,paste0("Q",as.integer(positive),"T",as.integer(a$TBDRUIND>=t)),NA_character_)
 fits <- merge(reference,subset(m,scenario==k & status=="ok",c(variable,beta)),by="variable")
 rows[[i]] <- data.frame(q,t,classified=sum(eligible),q_only=sum(state=="Q1T0",na.rm=T),
 t_only=sum(state=="Q0T1",na.rm=T),changed=sum(state!=a$group,na.rm=T),
 comparable=nrow(fits),same=sum(sign(fits$beta.x)==sign(fits$beta.y)))
}
z <- do.call(rbind,rows)
stopifnot(all(z$q_only==counts$qft_only),all(z$t_only==counts$tst_only))
write.csv(z,"Source_data/Figure1D_threshold_summary.csv",row.names=FALSE)
panel <- function(){
 grid.newpage()
 label <- function(x,y,s,b=FALSE,size=10,col="#142F40") grid.text(s,x,y,gp=gpar(fontfamily="Helvetica",fontsize=size,fontface=if(b)2 else 1,col=col))
 label(.04,.94,"D",TRUE,14)
 label(.48,.94,"QFT antigen minus nil (IU/mL)",TRUE)
 for(j in 1:3) label(.25+(j-1)*.29,.83,c("0.20","0.35","0.70")[j],TRUE)
 label(.07,.83,"TST (mm)",TRUE,9)
 for(i in 1:9){
  col <- (i-1)%%3; row <- (i-1)%/%3; x <- .25+col*.29; y <- .66-row*.235
  grid.rect(x,y,width=.274,height=.219,gp=gpar(fill=if(i==5)"#DCEBEA" else "#F1F5F6",col=NA))
  if(col==0)label(.07,y,z$t[i],TRUE)
  label(x,y+.055,sprintf("%d / %d",z$q_only[i],z$t_only[i]),TRUE,12)
  label(x,y-.015,sprintf("Changed: %d",z$changed[i]),size=8.7)
  label(x,y-.069,sprintf("Same sign: %d/%d",z$same[i],z$comparable[i]),size=8.7)
 }
 label(.49,.027,"Bold counts: QFT-only / TST-only",size=9)
}
export <- function(stem,h,fun){
 svglite(paste0(stem,".svg"),width=6.6,height=h);fun();dev.off()
 pdf(paste0(stem,".pdf"),width=6.6,height=h);fun();dev.off()
 png(paste0(stem,".png"),width=6.6,height=h,units="in",res=350,type="quartz");fun();dev.off()
}
export("Figures/Figure1D_thresholds",3.15,panel)
f <- subset(read.csv("Source_data/Figure3_all_34_discordance.csv"),layer=="module")

export("Figures/Figure3B_modules",2.1,function(){
 grid.newpage()
 x <- function(v) .39 + (.72-.39)*(log(v)-log(.4))/(log(2)-log(.4))
 label <- function(x,y,s,bold=FALSE,just="centre",size=9){
   grid.text(s,x,y,just=just,gp=gpar(fontfamily="Helvetica",fontsize=size,fontface=if(bold)2 else 1,col="#203441"))
 }
 segment <- function(x1,y1,x2,y2,...) grid.lines(c(x1,x2),c(y1,y2),gp=gpar(...))
 y <- c(.71,.56,.41,.26)
 label(.012,.9,"B",TRUE,size=13)
 label(.055,.9,"Fatty acid module",TRUE,"left")
 label(.315,.9,"n",TRUE)
 label(.55,.9,"Adjusted odds ratio",TRUE)
 label(.86,.9,"OR (95% CI)",TRUE)
 segment(x(1),.17,x(1),.80,col="#586877",lty=2)
 for(i in seq_along(y)){
   label(.055,y[i],c("Saturated","Monounsaturated","Omega-3","Omega-6")[i],just="left")
   label(.315,y[i],f$n[i])
   lo<-x(f$OR_low[i]);hi<-x(f$OR_high[i])
   segment(lo,y[i],hi,y[i],col="#256478",lwd=1.5)
   for(end in c(lo,hi)) segment(end,y[i]-.018,end,y[i]+.018,col="#256478")
   grid.points(x(f$OR[i]),y[i],default.units="npc",pch=15,size=unit(2.7,"mm"),gp=gpar(col="#174B5D"))
   label(.86,y[i],sprintf("%.2f (%.2f to %.2f)",f$OR[i],f$OR_low[i],f$OR_high[i]))
 }
 segment(x(.4),.14,x(2),.14,col="#84939B")
 for(tick in c(.5,1,2)){
   segment(x(tick),.14,x(tick),.115,col="#73838C")
   label(x(tick),.052,format(tick),size=8)
 }
})
print(z)
