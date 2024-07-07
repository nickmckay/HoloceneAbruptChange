#Code to reproduce Figure 1 of "The 4.2 ka event is not remarkable in the context of Holocene climate variability"

#Author: Chris Hancock

#Load packages
library(tidyverse)
library(egg)

RAWtheme <- function(base_size = 8, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(size = base_size),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size),
      plot.title = element_text(size = base_size),
      plot.subtitle = element_text(size = base_size),
      plot.caption = element_text(size = base_size)
    )
}

#Load Data  ----
df <- read.csv('metadata/Lit Review_Metanalysis - MetaAnalysis - NEW.csv')
#Set up data to find
res <- 0.1 #ka
bins<-seq(3.4,5,res)
binscount <- rep(0,length(bins))
ends <-begs <- c()
#Find info
for (i in 1:nrow(df)){
  beg <- as.numeric(df[i,'Excursion.start..ka.'])
  end <- as.numeric(df[i,'Excursion.end..ka.'])
  #If no dates, continue
  if ((is.na(end)) & (is.na(beg))){next}
  #If one date is missing, assume end=beginning
  if ((is.na(end)==TRUE) & (is.na(beg)==FALSE)){end<-beg}
  if ((is.na(end)==FALSE) & (is.na(beg)==TRUE)){beg<-end}
  #Make sure all values in ka
  if (beg > 1000){beg <- beg/1000}
  if (end > 1000){end <- end/1000}
  #Save begining and end in vector
  begs <- c(begs,beg)
  ends <- c(ends,end)
  #Id which timebins the 4.2 event covers
  for (t in 1:length(bins)){
    if (((bins[t]+res/2)>=end) & ((bins[t]-res/2)<=beg)){
      binscount[t] <- binscount[t]+1
    }
  }
}

#Plot Fig1 (histogram) ----
Fig1a <- ggplot() +
  #geom_vline(aes(xintercept = 4.2))+
  geom_histogram(aes(x=begs,fill=paste0('beginning (',median(begs),' ka)')),alpha=0.6,binwidth=0.1)+
  geom_vline(aes(xintercept=median(begs)),color='#1b9e77',size=1)+
  geom_histogram(aes(x=ends,fill=paste0('end (',median(ends),'.0 ka)')),alpha=0.6,binwidth=0.1)+
  geom_vline(aes(xintercept=median(ends)),color='#d95f02',size=1)+
  scale_x_reverse('age (ka)',limits=c(4.65,3.35),breaks=seq(3.4,4.6,0.1),expand=c(0,0),
                  labels=c(c('3.4','','','3.7','','','4.0','','','4.3','','','4.6')))+
  scale_y_continuous('count (# of papers)',limits=c(0,15),expand=c(0,0))+
  labs(title="4.2 ka event timing")+
  scale_fill_manual(values = c("#1b9e77", "#d95f02"))+
  RAWtheme()+
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = c(0.77,0.77),
    legend.title = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    legend.background = element_rect(fill='grey99',color=NA)
    )
Fig1a

Fig1b <- ggplot() +
  #geom_vline(aes(xintercept = 4.2))+
  geom_histogram(aes(x=(df$Duration..kyr)*1000),fill='#7570b3',
                 alpha=0.6,binwidth=100)+
  geom_vline(aes(xintercept=median(df$Duration..kyr,na.rm=T)*1000),color='#7570b3',size=1)+
  scale_x_continuous('duration (years)',limits=c(50,750),breaks=seq(100,700,100),expand=c(0,0),
                  labels=c(c('100','','300','','500','','700')))+
  scale_y_continuous('count (# of papers)',limits=c(0,15),expand=c(0,0))+
  labs(title="4.2 ka event timing")+
  RAWtheme()+
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    )
Fig1b

#Save
Fig1 <- egg::ggarrange(plots=list(Fig1a,Fig1b),ncol=2, labels=c('a)','b)'),
          label.args = list(gp = grid::gpar(fontsize = 9, fontfamily = "", fontface = "bold"), hjust = -0.5, vjust=1.5))
ggsave('Fig1.pdf',Fig1, dpi = 400, height = 2.5, width = 6.5)

#Plot FigS1 (map) -----
#dataS1a
dataS1a <- df[(df$Source %in% c("Literature search","Railsback et al. 2018","Search & Railsback et al. 2018")),]
dataS1a <- dataS1a[!is.na(dataS1a$Lat....),]
dataS1a <- dataS1a[!is.na(dataS1a$Lon....),]
dataS1a <- dataS1a %>% sf::st_as_sf(coords = c("Lon....", "Lat...."), crs = 4326)
dataS1a$Event.Type <- ifelse(dataS1a$Change.type == 'Excursion','excursion','other')
dataS1a$AuthorInterp1 <- NA
dataS1a$AuthorInterp2 <- NA
for (i in 1:nrow(dataS1a)){
  interp <- tolower(dataS1a$Climate.signal[i])
  if(interp %in% c("no event",'cold','warm','wet','dry')){
    dataS1a$AuthorInterp1[i] <- dataS1a$AuthorInterp2[i] <- interp
  }else if(grepl(' & ',interp)){
    dataS1a$AuthorInterp1[i] <- strsplit(interp,' & ')[[1]][1]
    dataS1a$AuthorInterp2[i] <- strsplit(interp,' & ')[[1]][2]
  }else{
    dataS1a$AuthorInterp1[i] <- dataS1a$AuthorInterp2[i] <- 'other'
  }
}
levs<-c('wet','dry','cold','warm','other','no event')
dataS1a$AuthorInterp1 <- factor(dataS1a$AuthorInterp1,levels=levs)
dataS1a$AuthorInterp2 <- factor(dataS1a$AuthorInterp2,levels=levs)
dataS1a$Event.Type <- factor(dataS1a$Event.Type,levels=c('excursion','other'))

#dataS2a
dataS2a <- df[(df$Source %in% c("Wang et al. 2016","Marchant & Hooghiemstra 2004")),]
dataS2a <- dataS2a[!is.na(dataS2a$Lat....),]
dataS2a <- dataS2a[!is.na(dataS2a$Lon....),]
dataS2a <- dataS2a %>% sf::st_as_sf(coords = c("Lon....", "Lat...."), crs = 4326)
dataS2a$Event.Type <- ifelse(dataS2a$Change.type == 'Excursion','excursion','other')
dataS2a$AuthorInterp1 <- NA
dataS2a$AuthorInterp2 <- NA
for (i in 1:nrow(dataS2a)){
  interp <- tolower(dataS2a$Climate.signal[i])
  if(interp %in% c("no event",'cold','warm','wet','dry')){
    dataS2a$AuthorInterp1[i] <- dataS2a$AuthorInterp2[i] <- interp
  }else if(grepl(' & ',interp)){
    dataS2a$AuthorInterp1[i] <- strsplit(interp,' & ')[[1]][1]
    dataS2a$AuthorInterp2[i] <- strsplit(interp,' & ')[[1]][2]
  }else{
    dataS2a$AuthorInterp1[i] <- dataS2a$AuthorInterp2[i] <- 'other'
  }
}
levs<-c('wet','dry','cold','warm','other','no event')
dataS2a$AuthorInterp1 <- factor(dataS2a$AuthorInterp1,levels=levs)
dataS2a$AuthorInterp2 <- factor(dataS2a$AuthorInterp2,levels=levs)
dataS2a$Event.Type <- factor(dataS2a$Event.Type,levels=c('excursion','other'))

#base +
# Plot the map using ggplot2
base <- ggplot(data = sf::st_as_sf(rworldmap::getMap(resolution = "low"))) +
  geom_sf(fill = "grey95",color='grey') +
  coord_sf(crs = "+proj=robin") +
  scale_y_continuous(breaks = c(-89.9,seq(-60, 60, by = 30),89.9)) +
  RAWtheme() +
  scale_shape_manual(values=c(21,24),drop=FALSE)+
  scale_fill_manual(values=c('forestgreen','goldenrod','dodgerblue','firebrick','grey','black'),drop=FALSE)+
  scale_color_manual(values=c('forestgreen','goldenrod','dodgerblue','firebrick','grey','black'),drop=FALSE)+
  guides(fill="none",  color=guide_legend(nrow=2, title.position = "top"), shape=guide_legend(nrow=2, title.position = "top")) +
  labs(color="author interpretation", shape="event type") +
  theme(panel.ontop = TRUE,
        panel.background = element_rect(fill = NA,color=NA),
        panel.border = element_rect(fill = NA,color=NA),
        panel.grid.major = element_line(color='grey80'),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
  )

FigS1a <- base+
  geom_sf(data=dataS1a,aes(fill=AuthorInterp1,color=AuthorInterp2,shape=Event.Type),
          stroke=0.8,size=0.9,alpha=0.7)+
  labs(title='Metanalysis (literature review and Railsback et al. (2018))')+
  coord_sf(crs = "+proj=robin")

FigS1b <- base+
  geom_sf(data=dataS2a,aes(fill=AuthorInterp1,color=AuthorInterp2,shape=Event.Type),
          stroke=0.9,size=0.9,alpha=0.7)+
  labs(title='Marchant and Hooghiemstra (2004) and Wang et al. (2016)')+
  coord_sf(crs = "+proj=robin")


#Save
FigS1 <- egg::ggarrange(plots=list(FigS1a+theme(legend.position = "none"),FigS1b),
                                   ncol=1, labels=c('a)','b)'),
                       label.args = list(gp = grid::gpar(fontsize = 9, fontfamily = "", fontface = "bold"),
                                         hjust = -0.5, vjust=1.5))
ggsave('FigS1.pdf',FigS1, dpi = 400, height = 6.5, width = 6.5)
#Map


#View(df)
#bm <-
