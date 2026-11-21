# From the release directory: Rscript code/33_graphical_abstract.R
library(grid)
read_result <- function(name) read.csv(file.path("B_results",paste0(name,".csv")))
flow <- read_result("completion_flow")
nativity <- function(name) subset(read_result(name), variable=="nativity" & level=="Foreign-born")
gap <- rbind(nativity("completion_groups"),nativity("completion_groups_documented"))
ipw <- read_result("completion_ipw188")
primary <- read_result("primary188")
both <- subset(ipw,primary_status=="ok" & status=="ok")
hit <- primary$variable[primary$q_188<.05]
stopifnot(identical(flow$n,c(5319L,5235L,4410L)),all(flow$paired==4215),
          nrow(ipw)==188,nrow(both)==156,setequal(hit,ipw$variable[ipw$q188<.05]))
same <- sum(sign(log(both$OR))==sign(log(both$primary_OR)))
bridge <- read_result("bridge6392")
lip <- with(subset(bridge,direction!="unavailable"),table(layer,direction))
stopifnot(same==154,length(hit)==3,sum(lip)==5270)
ink <- "#193748"; blue <- "#196C92"; green <- "#087867"; grey <- "#586B74"; rule <- "#D3DEE1"
txt <- function(x,y,s,size=11,col=ink,bold=FALSE,j="left"){
 grid.text(s,x,y,just=j,gp=gpar(fontfamily="Arial",fontsize=size,col=col,fontface=if(bold)"bold" else "plain",lineheight=1.05))
}
line <- function(x,y,xx,yy,col=rule,width=1,dash=1){
 grid.segments(x,y,xx,yy,gp=gpar(col=col,lwd=width,lty=dash,lineend="butt"))
}
panel <- function(i,name){
 pushViewport(viewport(x=(i-.5)/3,width=.31,height=.95))
 txt(0,.96,LETTERS[i],13,blue,TRUE);txt(.09,.96,name,14,ink,TRUE)
 line(0,.91,1,.91)
}
comma <- function(x)format(x,big.mark=",",trim=TRUE)
draw <- function(){
 grid.newpage()
 panel(1,"Paired QFT/TST completion")
 txt(0,.85,"US adults  ·  NHANES 2011–2012",10,grey)
 labs <- c("All examined adults","No documented severe TST reaction","Explicitly denied a prior severe reaction")
 for(i in 1:3){
  y<-c(.73,.49,.25)[i]
  txt(0,y,labs[i],11,ink,TRUE)
  txt(0,y-.07,paste(comma(flow$paired[i]),"/",comma(flow$n[i])),12)
  txt(1,y-.07,sprintf("%.1f%%",flow$pct[i]),15,blue,TRUE,"right")
  line(0,y-.12,1,y-.12,rule,4);line(0,y-.12,flow$pct[i]/100,y-.12,blue,4)
 }
 txt(0,.04,"Unweighted counts  ·  survey-weighted rates",9,grey);popViewport()
 panel(2,"Completion by birthplace")
 txt(0,.85,"Foreign-born minus US-born, percentage points",10,grey)
 trans <- function(x).04+(x+12)/16*.90
 line(trans(0),.25,trans(0),.70,grey,1,2)
 for(i in 1:2){
  y<-c(.65,.37)[i];r<-gap[i,]
  txt(0,y+.105,c("No documented severe reaction","Explicitly denied a severe reaction")[i],11,ink,TRUE)
  line(trans(r$diff_low),y,trans(r$diff_high),y,blue,2.5)
  grid.points(trans(r$difference),y,pch=16,size=unit(3,"mm"),gp=gpar(col=blue))
  txt(0,y-.065,sprintf("%+.1f  (%.1f to %.1f)",r$difference,r$diff_low,r$diff_high),11,blue,TRUE)
 }
 line(trans(-12),.19,trans(4),.19)
 for(k in seq(-12,4,4)){line(trans(k),.19,trans(k),.175,grey);txt(trans(k),.13,k,9,grey,j="centre")}
 txt(0,.04,"Points and lines show differences and 95% CIs.",9,grey);popViewport()
 panel(3,"Exposure and lipid findings")
 txt(0,.85,sprintf("%d biomarkers screened  ·  %d estimable",nrow(ipw),nrow(both)),10,grey)
 # A compact identity plot shows all paired, estimable marker-specific ORs.
 box <- c(.15,.56,.445,.765)
 coord <- function(x,a,b)a+(log10(x)+2.05)/3.2*(b-a)
 X<-function(x)coord(x,box[1],box[2]);Y<-function(y)coord(y,box[3],box[4])
 line(box[1],box[3],box[2],box[4],grey,1)
 line(box[1],box[3],box[2],box[3]);line(box[1],box[3],box[1],box[4])
 for(k in c(.01,.1,1,10)){txt(X(k),box[3]-.035,k,8,grey,j="centre");txt(box[1]-.02,Y(k),k,8,grey,j="right")}
 grid.points(X(both$primary_OR),Y(both$OR),pch=16,size=unit(1.25,"mm"),gp=gpar(col=adjustcolor(blue,.55)))
 txt(.15,.795,"Response-weighted OR",9,grey)
 txt(.15,.355,"Primary OR (log scales)",9,grey)
 txt(.64,.69,paste0(same,"/",nrow(both)),17,blue,TRUE)
 txt(.64,.60,"same direction",10)
 txt(.64,.50,paste0(length(hit)," of ",nrow(ipw)),13,blue,TRUE)
 txt(.64,.415,"passed FDR in\nboth analyses",9)
 line(0,.30,1,.30)
 txt(0,.245,"Alignment with discordance direction",11,ink,TRUE)
 for(i in 1:2){
  y<-c(.15,.065)[i];a<-lip[i,"aligned"];n<-sum(lip[i,])
  txt(0,y,c("Fatty acids","Lipid modules")[i],10)
  txt(.43,y,paste(comma(a),"/",comma(n)),10,grey)
  txt(1,y,sprintf("%.0f%%",100*a/n),13,green,TRUE,"right")
 };popViewport()
}
w<-12;h<-4.45
svglite::svglite("Graphical_Abstract.svg",width=w,height=h);draw();dev.off()
grDevices::quartz(type="pdf",file="Graphical_Abstract.pdf",width=w,height=h);draw();dev.off()
for(dpi in c(350,100)){
 name<-if(dpi==350)"Graphical_Abstract.png" else "Graphical_Abstract_preview.png"
 png(name,type="quartz",width=w,height=h,units="in",res=dpi,bg="white");draw();dev.off()
}
