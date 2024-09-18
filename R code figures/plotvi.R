library(ggplot2)
 cc <- read.table("histtndvig2_1",header=FALSE,sep = " ",dec=".")
 cc2 <- read.table("histtndvig2_2",header=FALSE,sep = " ",dec=".")
 cc3 <- read.table("histtndvig2_3",header=FALSE,sep = " ",dec=".")
 cc4 <- read.table("histtndvig2_4",header=FALSE,sep = " ",dec=".")
cc5 <- read.table("histtndvig2_5",header=FALSE,sep = " ",dec=".")
data<-(cbind(cc$V1,cc$V2,cc$V2-cc$V3,cc$V2+cc$V3,"Crop"))
data<-rbind(data,cbind(cc2$V1,cc2$V2,cc2$V2-cc2$V3,cc2$V2+cc2$V3,"Grassland"))
data<-rbind(data,cbind(cc3$V1,cc3$V2,cc3$V2-cc3$V3,cc3$V2+cc3$V3,"Broadleaf ever."))
data<-rbind(data,cbind(cc4$V1,cc4$V2,cc4$V2-cc4$V3,cc4$V2+cc4$V3,"Broadleaf dec."))
data<-rbind(data,cbind(cc5$V1,cc5$V2,cc5$V2-cc5$V3,cc5$V2+cc5$V3,"Needleleaf"))
data<-data.frame(data)
names(data)<-c("temperature","vi_delta","deltamin","deltamax","types")
data$temperature=as.numeric(data$temperature)
data$deltamin=as.numeric(data$deltamin)
data$deltamax=as.numeric(data$deltamax)
data$vi_delta=as.numeric(data$vi_delta)
data[data==0]<-NaN

myplot<-ggplot(data, aes(temperature, vi_delta, colour=types, fill=types)) +
  scale_colour_manual(values = c("Crop"="red", "Grassland"="orange" ,"Broadleaf ever."= "springgreen4","Broadleaf dec."="yellowgreen","Needleaf"="aquamarine4")) +geom_line(show.legend = FALSE)+
  geom_ribbon(aes(x=temperature, ymax=deltamax, ymin=deltamin, fill=types),colour = NA,alpha=0.2) +xlim(35,52)+ylim(-0.015,0.01)+theme_minimal()+scale_fill_discrete(name = "Veg. type") + xlab("maximum surface temperature (°C)")+ ylab(expression(delta~"HS_ref (NDVI)"))+geom_hline(yintercept=0, linetype="dashed", color = "black")+theme_bw()+theme(legend.position=c(0.2,0.2))
png("vigraph.png")
print(myplot)
dev.off()




