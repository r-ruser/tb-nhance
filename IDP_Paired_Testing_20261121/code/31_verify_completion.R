# Independent Taylor variance checks from public XPT data without survey estimators.
.libPaths(c("TB_metabolome_exposome_20261115_2130/analysis_MWAS_ExWAS/Rlib",.libPaths()))
library(haven)
raw <- Sys.glob("NHANES离线数据库/*/2011-2012")[1]
files <- c("Demographics/demo_g.xpt","Laboratory/tb_g.xpt","Examination/tbx_g.xpt")
data <- Reduce(function(x,y)merge(x,y,by="SEQN",all.x=TRUE),lapply(file.path(raw,files),read_xpt))
data <- subset(data,WTMEC2YR>0)
domain <- with(data,RIDAGEYR>=20 & !(TBQ070 %in% 1))
y <- with(data,as.numeric(LBXTBIN %in% 1:2 & !is.na(TBDRUIND)))
linear <- function(mask){
 mask[is.na(mask)] <- FALSE
 w <- data$WTMEC2YR*mask
 p <- sum(w*y)/sum(w)
 list(p=p,lin=w*(y-p)/sum(w))
}
variance <- function(lin){
 psu <- aggregate(lin,list(stratum=data$SDMVSTRA,psu=data$SDMVPSU),sum)
 sum(vapply(split(psu$x,psu$stratum),function(z)length(z)*var(z),numeric(1)))
}
out <- "IDP_Paired_Testing_20261121/B_results"
g <- read.csv(file.path(out,"completion_groups.csv"))
nat <- function(k)linear(domain & data$DMDBORN4==k)
us <- nat(1); foreign<-nat(2); contrast<-100*(foreign$p-us$p)
ci<-contrast+c(-1,1)*qt(.975,17)*100*sqrt(variance(foreign$lin-us$lin))
target<-subset(g,variable=="nativity" & level=="Foreign-born")
stopifnot(abs(contrast-target$difference)<1e-8,max(abs(ci-c(target$diff_low,target$diff_high)))<1e-8)
domains<-list(data$RIDAGEYR>=20,domain,domain & data$TBQ070 %in% 2)
f<-read.csv(file.path(out,"completion_flow.csv"))
for(i in seq_along(domains)){
 z<-linear(domains[[i]]); ci<-100*plogis(qlogis(z$p)+c(-1,1)*qt(.975,17)*sqrt(variance(z$lin))/(z$p*(1-z$p)))
 stopifnot(max(abs(c(100*z$p,ci)-unlist(f[i,c("pct","low","high")])))<1e-8)
}
cat("Independent XPT linkage, three weighted completion rates and intervals, and foreign-born minus US-born contrast and interval agree within 1e-8 percentage points.\n")
