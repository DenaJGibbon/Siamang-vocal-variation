library(stringr)
library(dplyr)
library(ggplot2)
library(lme4)
library(bbmle)
library(plyr)
library(flextable)
library(apcluster)
library(randomForest)
library(cluster)
library(clValid)

# V1 focus only on US2s that were annotated in the same way

files <- list.files('combine',recursive = T,
                    full.names = T,pattern = '.txt')
short.files <- list.files('combine',recursive = T,
                          full.names = F,pattern = '.txt')


# Code to read in selection tables
combined.US2.df <- data.frame()
for(a in 1:length(files)) { tryCatch({
  print(a) 
  temp.table <- read.delim2(files[a],stringsAsFactors = F)
  temp.table <- subset(temp.table,View=="Spectrogram 1")
  
  temp.table <- temp.table[order( as.numeric(temp.table$Begin.Time..s.)),]
  
  group <- str_split_fixed(short.files[a],pattern = 'Song',n=2)[,1]
  group.label <- str_split_fixed(short.files[a],pattern = '.txt',n=2)[,1]
  new.temp.table <-cbind.data.frame(temp.table,group,group.label)
  
  NoteIntervalList <- list()
  for(b in 1: (nrow(new.temp.table)-1 )){
    FirstRow <- new.temp.table[b,]
    SecondRow <- new.temp.table[b+1,]
    InterNoteInterval <- as.numeric(SecondRow$Begin.Time..s.) - as.numeric(FirstRow$End.Time..s.)
    NoteIntervalList[[b]] <- InterNoteInterval
  }
  
  NoteIntervalList[[ nrow(new.temp.table)]] <- 'NA'
  
  new.temp.table$IOI <- unlist(NoteIntervalList)
  
  combined.US2.df <- rbind.data.frame(combined.US2.df,new.temp.table)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



combined.US2.df$IOI <- as.numeric(combined.US2.df$IOI)
hist( as.numeric(combined.US2.df$IOI))

combined.US2.df <- droplevels(subset(combined.US2.df,IOI <= 5))
hist( as.numeric(combined.US2.df$IOI))
range( as.numeric(combined.US2.df$IOI))



TempUS <- subset(combined.US2.df, Sequence='us')
hist( as.numeric(TempUS$IOI))

# solution
combined.US2.df[,4:13] <- 
  combined.US2.df[,4:13] %>% mutate_if(is.character,as.numeric)

pair.index <- unique(combined.US2.df$group.label)

siamangUS2df <- data.frame()
for(b in 1:length(pair.index)){
  
     temp.table <- subset(combined.US2.df, group.label==pair.index[b])
  
      individual <- temp.table$group[1]
      note.dur <- mean(temp.table$End.Time..s. - temp.table$Begin.Time..s.)  
      call.dur  <- max(temp.table$End.Time..s.) - min(temp.table$Begin.Time..s.)
      call.id <- temp.table$group.label[1]
      nnotes <- nrow(temp.table)
      
       
        mean95 <-  mean(temp.table$Freq.95...Hz.)
        max95 <-  max(temp.table$Freq.95...Hz.)
        min95 <-  min(temp.table$Freq.95...Hz.)
        minbw <-  min(temp.table$BW.90...Hz.)
        maxbw <-  max(temp.table$BW.90...Hz.)
        meanbw <- mean(temp.table$BW.90...Hz.)
        mindurnote <-  min(temp.table$Dur.90...s.)
        maxdurnote <- max(temp.table$Dur.90...s.)
        noterate <- nrow(temp.table)/call.dur
        note1dur <- temp.table[1,]$Dur.90...s.

        note1maxfreq <- temp.table[1,]$Freq.95...Hz.
        note2dur <- temp.table[2,]$Dur.90...s.

        note2maxfreq <- temp.table[2,]$Freq.95...Hz.
        range.bw <- maxbw-minbw
        rest.dur <- call.dur - sum(temp.table$Dur.90...s.)
        lastnotedur <- temp.table[nrow(temp.table),]$Dur.90...s.
        lastnoteminfreq <- temp.table[nrow(temp.table),]$Freq.5...Hz.
        lastnotemaxfreq <- temp.table[nrow(temp.table),]$Freq.95...Hz.
        
        us2 <- temp.table[ which(temp.table$Sequence=='us'),]
        
        if( nrow(us2) >0){
          
        us2 <- us2[,c(10,11,13)]
        colnames(us2) <- c("US2Dur.90...s.","US2Freq.95...Hz.", "US2BW.90...Hz.")
        
        n.bk <- length(which(temp.table$Sequence=='bk'))
        n.bm <- length(which(temp.table$Sequence=='bm'))
        
        temp.US2.df <- cbind.data.frame(individual,call.id,call.dur,nnotes,
                                         min95, minbw,maxbw,mean95,max95,meanbw,mindurnote,
                                         maxdurnote,noterate,note1dur,note1maxfreq,
                                         note2dur,note2maxfreq,
                                         range.bw,rest.dur,
                                         lastnotedur,lastnotemaxfreq,
                                        us2,n.bk,n.bm)
        
        siamangUS2df <- rbind.data.frame(siamangUS2df,temp.US2.df)
        }
      }

siamangUS2df$individual <- as.factor(siamangUS2df$individual)

table(siamangUS2df$individual)

siamangUS2df$sequence <-
  as.factor(str_split_fixed(siamangUS2df$call.id, pattern = '[.]', n=2)[,2])
siamangUS2df$song <-
  as.factor(str_split_fixed(siamangUS2df$call.id, pattern = 'US', n=2)[,1])

siamangUS2df$site <- recode(siamangUS2df$individual,
                            Group1='Sikundar',Group2='Sikundar',
                            Group3='Kutapanjang',Group4='Kutapanjang',
                            Group5='Kutapanjang',Group6='Kutapanjang',
                            Group7='Kutapanjang',Group8='Kutapanjang',
                            Group9='Soraya',Group10='Soraya')


table(siamangUS2df$site)

siamangUS2df <- droplevels(subset(siamangUS2df,individual !='IndvB' ))

plot(siamangUS2df$call.dur ~ as.numeric(siamangUS2df$sequence))

# Supervised clustering using spectrogram features by site
ml.model.specfeatures.svm <- e1071::svm(siamangUS2df[,-c(1:3,27,28,29)], 
                                             siamangUS2df$individual, kernel = "radial", 
                                             cross = nrow(siamangUS2df))


ml.model.specfeatures.svm$tot.accuracy

ml.model.specfeatures.site <- e1071::svm(siamangUS2df[,-c(1:3,27,28,29)], 
                                        siamangUS2df$site, kernel = "radial", 
                                        cross = nrow(siamangUS2df))


ml.model.specfeatures.site$tot.accuracy

ml.model.specfeatures.site$accuracies
mean(ml.model.specfeatures.site$accuracies)
sd(ml.model.specfeatures.site$accuracies)


ml.model.specfeatures.sequence <- e1071::svm(siamangUS2df[,-c(1:3,27,28,29)], 
                                         siamangUS2df$sequence, kernel = "radial", 
                                         cross = nrow(siamangUS2df))


ml.model.specfeatures.sequence$tot.accuracy

intern <- clValid(siamangUS2df[,-c(1:3,27,28,29)], 2:8, clMethods=c("hierarchical","kmeans","pam"),
                  maxitems= nrow(siamangUS2df),validation="internal")
intern@labels
summary(intern)
optimalScores(intern)
plot(intern)

pca.vals <- princomp(siamangUS2df[,-c(1:2,25:28,29)])
temp.pca <- cbind.data.frame(pca.vals$scores[,1],pca.vals$scores[,2],intern@labels)
colnames(temp.pca) <- c('PCA1','PCA2','Cluster')
temp.pca$Cluster <- as.factor(temp.pca$Cluster)
ggpubr::ggscatter(data=temp.pca,
                  x='PCA1',y='PCA2',color='Cluster')+ theme(legend.position = "none")

## Part 2b. Using UMAP to visualize differences in individual males
male.individual.umap <- 
  umap::umap(siamangUS2df[,-c(1:2,25:28,29)],labels=as.numeric(siamangUS2df$individual),
             controlscale=TRUE,scale=3)

plot.for.male.individuals <-
  cbind.data.frame(male.individual.umap$layout[,1:2],
                   as.factor(siamangUS2df$individual))

colnames(plot.for.male.individuals) <-
  c("Dim.1", "Dim.2", "individual")

my_plot_umap_specfeatures_individ <-
  ggplot(data = plot.for.male.individuals, aes(
    x = Dim.1,
    y = Dim.2,
    colour = individual
  )) +
  geom_point(size = 3) +
  #geom_text(label=plot.for.male.individuals$individual) +
  #scale_color_continuous(values = matlab::jet.colors (length(unique(plot.for.male.individuals$individual)))) +
  theme_bw() + ggtitle('Spectrogram features') + xlab('UMAP: Dim 1')+ylab('UMAP: Dim 2')+
  labs(color='Pair')#+ theme(legend.position = "none")

my_plot_umap_specfeatures_individ

# Do the cluster in the US2 phrases occur at different points in the song?
# Run random forest
rf2.male <- randomForest(x = siamangUS2df[,-c(1:2,17:29)], 
                         ntree = 10000, proximity = TRUE)
rf2.male

# Run iteratively over multiple cluster solutions
prox <- rf2.male$proximity
n.clusters <- seq(2,60,5)
Male.sil.df <- data.frame()
for(a in 1:length(n.clusters)){
  pam.rf <- pam(prox, n.clusters[a])
  
  
  sil <-
    cluster::silhouette(x = pam.rf$clustering,
                        dist = dist(siamangUS2df[,-c(1:2,25:28,29)]))
  
  sil.val <- (summary(sil)$avg.width)
  temp.sil.df <-  cbind.data.frame(sil.val,n.clusters[a])
  Male.sil.df <- rbind.data.frame(Male.sil.df,temp.sil.df)
}

pam.rf <- pam(prox, Male.sil.df[which.max(Male.sil.df$sil.val),]$`n.clusters[a]`)

siamangUS2df.clustering <- 
  siamangUS2df

siamangUS2df.clustering$RF <- 
  pam.rf$clustering



# Run the UMAP algorithm which is used for data visualization 
hornbill.umap.affinity <- 
  umap::umap(siamangUS2df[,-c(1:2,25:28)],
             controlscale=TRUE,scale=3, labels=as.factor(affinity.df@idx))

# Create a new dataframe that has the UMAP coordinates 
hornbill.umap.affinity <- as.data.frame(hornbill.umap.affinity$layout)


# Change the cluster ID to a factor instead of a number
hornbill.umap.affinity$AffinityID <- 
  as.factor(affinity.df@idx)


# Add more informative column names
colnames(hornbill.umap.affinity) <-
  c("Dim.1", "Dim.2", "AffinityID")

# Plot our results

my_plot_hornbill.ClusterID <- 
  ggpubr::ggscatter(data = hornbill.umap.affinity,
                    x = 'Dim.1',
                    y = 'Dim.2',
                    color = 'AffinityID',size=3) + 
  scale_color_manual(values = matlab::jet.colors (length(unique(hornbill.umap.affinity$AffinityID)))) +
  theme_bw() + ggtitle('Siamang US2 (Affinity)') + 
  xlab('UMAP: Dim 1')+ylab('UMAP: Dim 2')+theme(legend.position = "none") 




qvals <- seq(0,1,0.1)
Male.sil.df.aff <- data.frame()
for(a in 1:length(qvals)){
  cluster.df <- apcluster::apcluster(
    negDistMat(r = 2),
    siamangUS2df[,-c(1:2,25:28)],
    q=qvals[a]
  )
  
  sil <-
    cluster::silhouette(x = cluster.df@idx,
                        dist = dist(siamangUS2df[,-c(1:2,25:28)]))
  
  sil.val <- (summary(sil)$avg.width)
  temp.sil.df <-  cbind.data.frame(sil.val,qvals[a])
  Male.sil.df.aff <- rbind.data.frame(Male.sil.df.aff,temp.sil.df)
}

affinity.df <- apcluster::apcluster(
  negDistMat(r = 2),
  siamangUS2df[,-c(1:2,25:28)], 
  q=qvals[which.max(Male.sil.df.aff$sil.val)]
)

length(unique(affinity.df@exemplars))

subset(cbind.data.frame(siamangUS2df,affinity.df@idx),individual=='Group1')
affinty <- affinity.df@idx
temp.df <- cbind.data.frame(siamangUS2df,affinty)

ggpubr::ggboxplot(data=temp.df,x="affinty",y="call.dur")
ggpubr::ggboxplot(data=temp.df,x="affinty",y="US2Freq.95...Hz.")

# Run the UMAP algorithm which is used for data visualization 
hornbill.umap.affinity <- 
  umap::umap(siamangUS2df[,-c(1:2,25:29)],
             controlscale=TRUE,scale=3, labels=as.factor(affinity.df@idx))

# Create a new dataframe that has the UMAP coordinates 
hornbill.umap.affinity <- as.data.frame(hornbill.umap.affinity$layout)


# Change the cluster ID to a factor instead of a number
hornbill.umap.affinity$AffinityID <- 
  as.factor(affinity.df@idx)


# Add more informative column names
colnames(hornbill.umap.affinity) <-
  c("Dim.1", "Dim.2", "AffinityID")

# Plot our results

my_plot_hornbill.ClusterID <- 
  ggpubr::ggscatter(data = hornbill.umap.affinity,
                    x = 'Dim.1',
                    y = 'Dim.2',
                    color = 'AffinityID',size=3) + 
  scale_color_manual(values = matlab::jet.colors (length(unique(hornbill.umap.affinity$AffinityID)))) +
  theme_bw() + ggtitle('Siamang US2 (Affinity)') + 
  xlab('UMAP: Dim 1')+ylab('UMAP: Dim 2')+theme(legend.position = "none") 

my_plot_hornbill.ClusterID

siamangUS2df$sequence <- as.numeric(siamangUS2df$sequence)

ggpubr::ggscatter(data=siamangUS2df,
                  y= 'call.dur', x='sequence',facet.by = 'individual')
