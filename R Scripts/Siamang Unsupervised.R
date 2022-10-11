library(randomForest)
library(cluster)

# Load data ---------------------------------------------------------------


# Male unsupervised clustering --------------------------------------------


nrow(siamangUS2df)

siamangUS2df[,-c(1:3,27:30)] <-
  siamangUS2df[,-c(1:3,27:30)] %>% mutate_if(is.character,as.numeric)


# 
# # Run affinity prop clustering
# q.val.seq <- seq(from=0.1,to=0.9,by=0.1)
# 
# 
# AcousticSignal.sil.df <- data.frame()
# for(a in 1:length(q.val.seq)){
#   print(a)
#   AcousticSignalsAP <-
#     apcluster::apcluster(negDistMat(r=2),q=q.val.seq[a],
#                          siamangUS2df[,-c(1:3,27:30)],
#                          maxits=100000,convits=10000)
#   
#   
#   sil <-
#     cluster::silhouette(x = AcousticSignalsAP@idx,
#                         dist = dist(siamangUS2df[,-c(1:3,27:30)]))
#   
#   sil.val <- (summary(sil)$avg.width)
#   temp.sil.df <-  cbind.data.frame(sil.val,q.val.seq[a])
#   AcousticSignal.sil.df <- rbind.data.frame(AcousticSignal.sil.df,temp.sil.df)
# }
# 
# MaxSil <- which.max(AcousticSignal.sil.df$sil.val)
# 

AcousticSignalsAP <-
  apcluster::apcluster(negDistMat(r=2),q=0.1,# q.val.seq[MaxSil],
                       siamangUS2df[,-c(1:3,27:30)],
                       maxits=100000,convits=10000)

#print(q.val.seq[MaxSil])
print(paste('N clusters=', length(AcousticSignalsAP@exemplars)))
#AcousticSignalsMFCCs$class <- as.factor(AcousticSignalsAP@idx)

SiamangUS2.umap <-
  umap::umap(siamangUS2df[,-c(1:3,27:30)],labels= AcousticSignalsAP@idx,
             controlscale=TRUE,scale=3)

plot.for.SiamangUS2s <-
  cbind.data.frame(SiamangUS2.umap$layout[,1:2],
                   siamangUS2df$site)

colnames(plot.for.SiamangUS2s) <-
  c("Dim.1", "Dim.2", "site")


# intern <- clValid(siamangUS2df[,-c(1:3,27,28,29)], 2:8, clMethods=c("hierarchical","kmeans","pam"),
#                   maxitems= nrow(siamangUS2df),validation="internal")
# intern@clusterObjs$hierarchical$order
# summary(intern)
# optimalScores(intern)


# clusterhierarchical <- clusters(intern,"pam")
# two <- cutree(clusterhierarchical,3)


plot.for.SiamangUS2s$Clusters <- as.factor(AcousticSignalsAP@idx)
site <- siamangUS2df$site
pair <- siamangUS2df$individual
Male.site <- ggplot(plot.for.SiamangUS2s, aes(x = Dim.1,
                                                   y = Dim.2, col = site)) +
  geom_point(size = 3) +
  scale_color_manual(values = matlab::jet.colors (length(unique(plot.for.SiamangUS2s$site))))+labs(color='Site')+
  ggtitle('A. US-II sites') +theme_bw()+theme(legend.position = "none") +xlim(-12,12)+ylim(-10,10)

Male.site

Male.cluster <- ggplot(plot.for.SiamangUS2s, aes(x = Dim.1,
                                                      y = Dim.2, col = Clusters)) +
  geom_point(size = 3) +
  scale_color_manual(values = matlab::jet.colors (length(unique(AcousticSignalsAP@idx))))+
  ggtitle('C. US-II unsupervised clustering') +
  theme_bw()+theme(legend.position = "none") +xlim(-12,12)+ylim(-10,10)

Male.cluster

Male.pair <- ggplot(plot.for.SiamangUS2s, aes(x = Dim.1,
                                              y = Dim.2, col = pair)) +
  geom_point(size = 3) +
  scale_color_manual(values = matlab::jet.colors (length(unique(pair))))+labs(color='Pair')+
  ggtitle('B. US-II pairs') +theme_bw()+theme(legend.position = "none") +xlim(-12,12)+ylim(-10,10)

Male.pair

cowplot::plot_grid(Male.site,Male.pair,Male.cluster,nrow=2)

Extremes.umap <- siamangUS2df[c(which.max(plot.for.SiamangUS2s$Dim.1),
               which.min(plot.for.SiamangUS2s$Dim.1),
               which.max(plot.for.SiamangUS2s$Dim.2),
               which.min(plot.for.SiamangUS2s$Dim.2)),]

AcousticSignalsAP@idx[c(which.max(plot.for.SiamangUS2s$Dim.1),
                        which.min(plot.for.SiamangUS2s$Dim.1),
                        which.max(plot.for.SiamangUS2s$Dim.2),
                        which.min(plot.for.SiamangUS2s$Dim.2))]

