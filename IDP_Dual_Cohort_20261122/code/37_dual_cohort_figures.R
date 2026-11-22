# Run from the submission directory. All inputs contain aggregate estimates.
library(grid)
read <- function(x) read.csv(file.path("B_results",paste0(x,".csv")))
rates <- read("external_completion")
groups <- c("NHANES no documented reaction","NHANES documented no reaction","TBESC enrolled")
gap <- subset(rates, grepl("minus",contrast))
gap <- gap[match(groups,gap$cohort),]
models <- subset(read("external_common_models"),term=="foreign")
models <- models[match(c("NHANES","TBESC"),models$cohort),]
ipw <- read("completion_ipw188")
both <- subset(ipw,status=="ok" & primary_status=="ok")
same <- sum(sign(log(both$OR))==sign(log(both$primary_OR)))
lip <- with(subset(read("bridge6392"),direction!="unavailable"),table(layer,direction))
bmi <- subset(read("external_common_models"),term=="bmi5")
bmi <- bmi[match(models$cohort,bmi$cohort),]
stopifnot(nrow(ipw)==188,nrow(both)==156,same==154,all(round(gap$estimate,1)==c(-6.3,.2,6.2)))
navy <- "#193748"; blue <- "#28658A"; teal <- "#007465"; grey <- "#53656D"; pale <- "#D3DEE1"
tx <- function(x,y,s,size=10,col=navy,bold=FALSE,j="left")grid.text(s,x,y,just=j,gp=gpar(fontfamily="Arial",fontsize=size,col=col,fontface=if(bold)"bold" else "plain",lineheight=1.03))
ln <- function(x,y,xx,yy,col=pale,w=1,style=1)grid.segments(x,y,xx,yy,gp=gpar(col=col,lwd=w,lty=style))
panel <- function(x,y,w,h)pushViewport(viewport(x=x,y=y,width=w,height=h,just=c("left","bottom")))
headpanel <- function(s){tx(0,.96,s,12,navy,TRUE);ln(0,.88,1,.88)}
signed <- function(x)sub("-","−",sprintf("%+.1f",x),fixed=TRUE)
fmt <- function(x)sub("-","−",sprintf("%.1f",x),fixed=TRUE)
ci <- function(r)paste0(signed(r$estimate),"  (",fmt(r$low),", ",fmt(r$high),")")
interval <- function(x,y,lo,hi,col,pch){ln(lo,y,hi,y,col,2.4);grid.points(x,y,pch=pch,size=unit(2.6,"mm"),gp=gpar(col=col,fill=col))}
forest <- function(kind){
 isgap <- kind=="gap"
 d <- if(isgap)gap else models
 labs <- if(isgap)c("NHANES · broad eligibility","NHANES · history documented","TBESC · enrolled") else c("NHANES","TBESC")
 ys <- if(isgap)c(.65,.43,.21) else c(.60,.30)
 map <- if(isgap) function(x).07+.86*(x+12)/24 else function(x).07+.86*(log(x)-log(.07))/(log(1.2)-log(.07))
 ref<-map(if(isgap)0 else 1)
 for(y in ys)ln(ref,y-.055,ref,y+.052,grey,1,2)
 for(i in seq_along(ys)){
  y<-ys[i];r<-d[i,];co<-if(grepl("TBESC",labs[i]))teal else blue
  tx(0,y+.09,labs[i],9,co,TRUE)
  interval(map(if(isgap)r$estimate else r$OR),y,map(r$low),map(r$high),co,if(co==teal)15 else 16)
  label<-if(isgap)ci(r) else sprintf("%.2f  (%.2f, %.2f)",r$OR,r$low,r$high)
  tx(1,y+.09,label,9,navy,FALSE,"right")
 }
 ticks<-if(isgap)c(-12,0,12) else c(.1,.25,.5,1)
 for(k in ticks){xx<-map(k);tx(xx,.045,format(k,trim=TRUE),8,grey,j="centre")}
 if(!isgap)tx(.5,.142,"Both estimates favor TST-only in foreign-born adults",8.5,navy,TRUE,"centre")
}
flow <- read("completion_flow")
flowtext <- list(
 N=c(paste0(format(flow$n[1],big.mark=",")," examined\nage ≥20"),"5,235 eligible*\n4,215 paired","188 exposures\n30 fatty acids + 4 modules"),
 T=c("22,020 enrolled\n18 clinics","14,233 eligible\nage ≥25","13,546 paired\n3,378 discordant"))
# Write the same directed paths and node labels used in the figure as editable Mermaid.
mmd<-c("flowchart LR")
for(k in names(flowtext)){
 labels<-flowtext[[k]]
 for(i in seq_along(labels))mmd<-c(mmd,sprintf('  %s%d["%s"]',k,i,gsub("\n","<br/>",labels[i])))
 mmd<-c(mmd,paste0("  ",paste0(k,1:3,collapse=" --> ")))
}
writeLines(mmd,"Figures/Figure1_flow.mmd")
flowpanel <- function(){
 headpanel("A   Who contributed to each analysis?")
 for(i in 1:2){
  k<-names(flowtext)[i];y<-c(.63,.28)[i];co<-c(blue,teal)[i]
  tx(0,y,c("NHANES","TBESC")[i],11,co,TRUE)
  for(j in 1:3){
   x<-c(.26,.55,.85)[j]
   grid.roundrect(x,y,.252,.245,r=unit(1.5,"mm"),gp=gpar(fill="white",col=co))
   tx(x,y,flowtext[[k]][j],10,navy,j="centre")
   if(j<3)grid.segments(x+.127,y,x+.159,y,arrow=arrow(length=unit(1.5,"mm"),type="closed"),gp=gpar(col=grey))
  }
 }
 tx(0,.055,"* No documented severe prior TST reaction. Between-cohort comparisons use age ≥25 in both cohorts.",9,grey)
}
notes <- function(){
 headpanel("C   Environmental and host profiles")
 tx(0,.79,"NHANES exposure screen",11,blue,TRUE)
 tx(0,.68,sprintf("%d screened  ·  %d estimable",nrow(ipw),nrow(both)),11)
 tx(0,.58,sprintf("%d retained direction after response weighting",same),9.3)
 tx(0,.48,"BMI per 5 kg/m² · adjusted discordance OR",10,navy,TRUE)
 for(i in 1:2)tx(0,.40-(i-1)*.075,sprintf("%s   %.2f (%.2f, %.2f)",models$cohort[i],bmi$OR[i],bmi$low[i],bmi$high[i]),10)
 ln(0,.25,1,.25)
 tx(0,.19,"Lipid sign agreement · NHANES",10,navy,TRUE)
 tx(0,.09,sprintf("Fatty acids %.0f%%   ·   modules %.0f%%",100*lip[1,"aligned"]/sum(lip[1,]),100*lip[2,"aligned"]/sum(lip[2,])),11,teal,TRUE)
}
draw <- function(ga){
 grid.newpage()
 if(ga){
  for(i in 1:3){panel((i-1)/3+.014,.035,.305,.93)
   if(i==1){headpanel("A   Paired-test completion");tx(0,.82,"Foreign-born minus US-born, pp (95% CI)",9,grey);forest("gap")}
   if(i==2){headpanel("B   Direction of discordance");tx(0,.82,"QFT-only vs TST-only · adjusted OR (95% CI)",9,grey);forest("or")}
   if(i==3)notes()
   popViewport()}
 }else{
  panel(.04,.66,.92,.32);flowpanel();popViewport()
  panel(.04,.31,.44,.32);headpanel("B   Completion by birthplace");tx(0,.82,"Foreign-born minus US-born, pp (95% CI)",9,grey);forest("gap");popViewport()
  panel(.52,.31,.44,.32);headpanel("C   Discordance by birthplace");tx(0,.82,"QFT-only vs TST-only · adjusted OR (95% CI)",9,grey);forest("or");popViewport()
  panel(.04,.025,.92,.255);headpanel("D   Host characteristics and robustness of the exposure screen")
  tx(0,.73,"BMI and discordance",11,blue,TRUE);tx(.52,.73,"NHANES exposure and lipid profiles",11,blue,TRUE)
  tx(0,.52,"Adjusted OR per 5 kg/m² (95% CI)",10)
  for(i in 1:2)tx(0,.32-(i-1)*.175,sprintf("%s   %.2f (%.2f, %.2f)",models$cohort[i],bmi$OR[i],bmi$low[i],bmi$high[i]),11)
  tx(.52,.52,sprintf("%d of %d estimable in both screens",nrow(both),nrow(ipw)),11)
  tx(.52,.32,sprintf("%d of %d retained direction after weighting",same,nrow(both)),10)
  tx(.52,.145,sprintf("Lipid sign agreement: %.0f%% fatty acids · %.0f%% modules",100*lip[1,"aligned"]/sum(lip[1,]),100*lip[2,"aligned"]/sum(lip[2,])),10)
  popViewport()
 }
}
for(ga in c(FALSE,TRUE)){
 path<-if(ga)"Graphical_Abstract" else "Figures/Figure1_complete"
 w<-if(ga)12 else 8.2;h<-if(ga)3.85 else 7.2
 svglite::svglite(paste0(path,".svg"),width=w,height=h);draw(ga);dev.off()
 quartz(type="pdf",file=paste0(path,".pdf"),width=w,height=h);draw(ga);dev.off()
 for(dpi in c(350,100)){
  png(paste0(path,if(dpi==100)"_preview.png" else ".png"),type="quartz",width=w,height=h,units="in",res=dpi,bg="white");draw(ga);dev.off()
 }
}
