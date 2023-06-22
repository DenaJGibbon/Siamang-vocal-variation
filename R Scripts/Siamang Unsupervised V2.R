library(randomForest)
library(cluster)
library(ggpubr)

# Load data ---------------------------------------------------------------


# Siamang unsupervised clustering --------------------------------------------


nrow(siamangUS2df)

siamangUS2df[,-c(1:3,27:30)] <-
  siamangUS2df[,-c(1:3,27:30)] %>% mutate_if(is.character,as.numeric)



intern <- clValid(siamangUS2df[,-c(1:3,27:30)], 3:4, clMethods=c("hierarchical","kmeans","pam","som","Mclust"),
                  maxitems= nrow(siamangUS2df),validation="internal")
intern@clusterObjs$hierarchical$order
summary(intern)
intern@measures
optimalScores(intern)



clusterhierarchical <- clusters(intern,"kmeans")

two <-clusterhierarchical$`3`$cluster

SiamangUS2.umap <-
  umap::umap(siamangUS2df[,-c(1:3,27:30)],labels= two,
             controlscale=TRUE,scale=3)

plot.for.SiamangUS2s <-
  cbind.data.frame(SiamangUS2.umap$layout[,1:2],
                   siamangUS2df$site)

colnames(plot.for.SiamangUS2s) <-
  c("Dim.1", "Dim.2", "site")



plot.for.SiamangUS2s$Clusters <- as.factor(two)
site <- siamangUS2df$site
pair <- siamangUS2df$individual
Siamang.site <- ggplot(plot.for.SiamangUS2s, aes(x = Dim.1,
                                                   y = Dim.2, col = site)) +
  geom_point(size = 3) +
  scale_color_manual(values = matlab::jet.colors (length(unique(plot.for.SiamangUS2s$site))))+labs(color='Site')+
  ggtitle('A. US-II sites') +theme_bw()+theme(legend.position = "none") #+xlim(-12,12)+ylim(-10,10)

Siamang.site

Siamang.cluster <- ggplot(plot.for.SiamangUS2s, aes(x = Dim.1,
                                                      y = Dim.2, col = Clusters)) +
  geom_point(size = 3) +
  scale_color_manual(values = matlab::jet.colors (length(unique(plot.for.SiamangUS2s$Clusters))))+
  ggtitle('C. US-II unsupervised clustering') +
  theme_bw()+theme(legend.position = "none") #+xlim(-12,12)+ylim(-10,10)

Siamang.cluster

Siamang.pair <- ggplot(plot.for.SiamangUS2s, aes(x = Dim.1,
                                              y = Dim.2, col = pair)) +
  geom_point(size = 3) +
  scale_color_manual(values = matlab::jet.colors (length(unique(pair))))+labs(color='Pair')+
  ggtitle('B. US-II pairs') +theme_bw()+theme(legend.position = "none")# +xlim(-12,12)+ylim(-10,10)

Siamang.pair

cowplot::plot_grid(Siamang.site,Siamang.pair,Siamang.cluster,nrow=2)

siamangUS2dfunsup <-siamangUS2df
siamangUS2dfunsup$cluster <- two

colnames(siamangUS2df[,-c(1:3,27:30)])
ggpubr::gghistogram(x= 'US2Dur.90...s.',data=siamangUS2dfunsup, facet.by = 'cluster')


