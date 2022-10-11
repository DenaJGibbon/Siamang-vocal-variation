# R code for Chapter LD
# Sound files available upon request
# Part 1 load libraries, read files, list files, get file names. (lines 8-22)
# Part 2 combine data frame, check structure, prepare classification (lines 80-93)
# Part 3 SVM and Check accuracy of SVM classification  (lines 32-47)
# Part 4 Convert to string and matrix (lines 96-154)
# Part 5 Finding means, correlation plots, and statistical analysis (lines 156-240)
library(stringdist)
library(stringr)
library(e1071)
library(matrixStats)
library(corrplot)

# Set directory location for all files
input.dir <- '/Users/Lomez/Desktop/combine'

# List files
ch4phrases <- list.files(input.dir,full.names = T, recursive = T,
                         pattern = '.txt')
# Get file names
ch4phrasesnames <- list.files(input.dir,full.names = F,recursive = T,
                         pattern = '.txt')

# Loop to combine all .txt into dataframe; if there is an error it will print out the index
alltextdf <- data.frame()
for(a in 1:length(ch4phrases)){tryCatch({
 print(a)
 temp.text <- ch4phrases[[a]]
 names <- ch4phrasesnames[a]
 #names <- stringr::str_split_fixed(names,pattern='/',n=3)[,3]
 names <- (stringr::str_split_fixed(names,pattern='.txt',n=2)[,1])
 tempdf <- read.delim(temp.text,stringsAsFactors = F)

 tempdf<- cbind.data.frame(tempdf,names)
 alltextdf <- rbind.data.frame(alltextdf,tempdf)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

}

# Check structure of dataframe
str(alltextdf)

# Create group id column
alltextdf$group <- str_split_fixed(alltextdf$names, pattern = 'Song',n=2)[,1]

# Isolate features for classification
call.index <- unique(alltextdf$names)
#length of call.index
length(call.index)
# Loop to estimate features from each phrase
feature.df <- data.frame()
for(i in 1:length(call.index)){
 print( paste(i, 'out of', length(call.index)))
 temp.call.df <-  subset(alltextdf, names==call.index[i])
 nnotes <- nrow(temp.call.df)
 minnotedur <- min(temp.call.df$Dur.90...s.)
 maxnotedur <- max(temp.call.df$Dur.90...s.)
 minfreq <- min(temp.call.df$Freq.95...Hz.)
 maxfreq <- max(temp.call.df$Freq.95...Hz.)
 maxbw <- max(temp.call.df$BW.90...Hz.)
 minbw <- min(temp.call.df$BW.90...Hz.)
 calldur <- max(temp.call.df$End.Time..s.) - min(temp.call.df$Begin.Time..s.)
 Nbm <- length( which(temp.call.df$Sequence=='bm')=='TRUE') 
 Nbk <-  length( which(temp.call.df$Sequence=='bk')=='TRUE') 
 Nus <-  temp.call.df[which(temp.call.df$Sequence=='us'),]$Dur.90...s. 
 names <- temp.call.df$names[1]
 group <- temp.call.df$group[1]
 features <- cbind.data.frame(names,group,nnotes,minnotedur,maxnotedur,minfreq,maxfreq,maxbw,minbw,calldur,
                              Nbm, Nbk,Nus)
 feature.df <- rbind.data.frame(feature.df,features)
}

apply(feature.df, 2,class)
nnum<-names(feature.df[,-c(1,2)])
feature.df$group<-as.factor(feature.df$group)
feature.df<-data.table::setDT(feature.df)[,(nnum):=lapply(.SD,as.numeric),.SDcols=nnum]
sapply(feature.df,class)

# SVM Model
svm.model <-
  svm(
    feature.df[,-c(1,2)],
    feature.df$group,
    kernel = "radial",
    cross = 10
  )

# Check accuracy of SVM classification
# Chance would be 1/12 and high is over 90% - possible factors are aggregating data, having 11 levels, possible predictors were left out. 
svm.model$tot.accuracy
summary(svm.model)
# Identify all unique names
call.index <- (unique(alltextdf$names))

# Convert each sequence for US2 to string for use in LD distance
sequence.list <- list()

for(b in 1:length(call.index)){
  tempseqdf <- subset(alltextdf,names==call.index[b])

  sequence <- c(tempseqdf$Sequence)
  call <- call.index[b]


  sequence.list[[b]] <- paste( unlist(sequence), collapse='|')
  }



# Convert to matrix
matrix.output <- stringdistmatrix(sequence.list)
matrix.output <- as.matrix(matrix.output, labels=TRUE)
row.names(matrix.output) <- call.index
colnames(matrix.output) <- call.index

# Start of code from Gamba's group
M <- matrix.output
#create the dimnames for the matrix
dimnames(M)<-list(NULL,NULL)
id<- str_split_fixed(call.index, pattern = 'Song',n=2)[,1]
dimnames(M)<-list(id,id)

# this function creates a matrix containing the mean values of an original matrix
AverageMatValsFast <- function(mat) {
  uniRow <- unique(rownames(mat))
  uniCol <- unique(colnames(mat))
  newmat1 <- matrix(numeric(0), nrow=length(uniRow), ncol=ncol(mat))
  rownames(newmat1) <- uniRow
  colnames(newmat1) <- colnames(mat)
  
  newmat <- matrix(numeric(0), nrow=length(uniRow), ncol=length(uniRow))
  rownames(newmat) <- uniRow
  colnames(newmat) <- uniCol
  
  for (i in 1:nrow(newmat1)) {
    rowMatch <- which(rownames(mat)==uniRow[i])
    if (length(rowMatch)>1) {
      newmat1[i,] <- colMeans(mat[rowMatch,])
    } else {
      newmat1[i,] <- mat[rowMatch,]
    }
  }
  
  for (j in 1:ncol(newmat)) {
    colMatch <- which(colnames(mat)==uniCol[j])
    
    if (length(rowMatch)>1) {
      newmat[,j] <- round(rowMeans(newmat1[,colMatch]),3)
    } else {
      newmat[,j] <- round(newmat1[,colMatch], 3)
    }
  }
  newmat
}
# calculate the mean for our matrix
av_M_sing_D<-AverageMatValsFast(M)

# Create a corrplot
corrplot::corrplot (av_M_sing_D,is.corr = F,cl.lim = c(4,66))

# this function creates a matrix containing the standard deviation 
SDMatValsFast <- function(mat) {
  uniRow <- unique(rownames(mat))
  uniCol <- unique(colnames(mat))
  newmat1 <- matrix(numeric(0), nrow=length(uniRow), ncol=ncol(mat))
  rownames(newmat1) <- uniRow
  colnames(newmat1) <- colnames(mat)
  
  newmat <- matrix(numeric(0), nrow=length(uniRow), ncol=length(uniRow))
  rownames(newmat) <- uniRow
  colnames(newmat) <- uniCol
  
  for (i in 1:nrow(newmat1)) {
    rowMatch <- which(rownames(mat)==uniRow[i])
    if (length(rowMatch)>1) {
      newmat1[i,] <- colSds(mat[rowMatch,])
    } else {
      newmat1[i,] <- mat[rowMatch,]
    }
  }
  
  for (j in 1:ncol(newmat)) {
    colMatch <- which(colnames(mat)==uniCol[j])
    
    if (length(rowMatch)>1) {
      newmat[,j] <- round(rowSds (newmat1[,colMatch]),3)
    } else {
      newmat[,j] <- round(newmat1[,colMatch], 3)
    }
  }
  newmat
}

# calculate the sd for our matrix
av_SD_sing<-SDMatValsFast(M)

corrplot::corrplot (av_SD_sing,is.corr = F)

summary(feature.df$calldur)
Test<-av_M_sing_D
rownames(Test)<-NULL
Test<-data.frame(Test)

#Organize LD values for statistical analysis

bbb <- function(x){
  x_col <- unique(colnames(x))
  x_row <- unique(rownames(x))
  x_col1<-paste0(x_col,"_")
  
  mm <- vector(mode = "list", length = length(x_col))
  
  for(i in 1:length(x_col)){
    mm[[i]] <- x[rownames(x) %in% x_row[i],colnames(x) %in% x_col[i]]
    }  
  names(mm) <- x_col1
  mm
  }

vv <- bbb(x=M)
vb <- lapply(vv, c)
yy<-tibble::enframe(unlist(vb, recursive = FALSE))

yy$Group<-str_extract(yy$name,pattern = (".+(?=_)"))
unique(yy$Group)
library(ggplot2)
ggplot(data=yy)+
  geom_histogram(mapping = aes(x=value))+
  facet_wrap(~Group)
#Welch one way ANOVA not assuming equal variances - because the box plot shows unequal variances. 
oneway.test(value~Group,data = yy)
boxplot(value~Group,data = yy)

#Paired test to see which groups are different.

yy$Group<-as.factor(yy$Group)
write.csv(yy,"US2.csv")
install.packages("haven")
library(haven)
write_sav(yy,"US2Phrase.sav")

###############################################################
#Graphs of LSI scores within and between

BetweenGroups<-c(31.234,28.249,24.843,26.594,40.019,15.691,24.712,25.618,24.46,36.672,4.91)
Mean_BG<-mean(BetweenGroups)
SD_BG<-sd(BetweenGroups)
WithinGroups<-c(30.712,	30.857,	29.782,	36.645,	31.522,	29.314,	30.913,	32.115,	37.857,	49.415,	27.543,	29.736,	34.878,	26.201,	27.152,	27.62,	28.07,	39.298,	42.187,	30.709,	34.61,	22.183,	26.54,	25.619,	25.212,	41.441,	34.189,	35.211,	32.727,	27.563,	30.389,	32.282,	35.268,	52.723,	33.649,	33.137,	34.181,	34.868,	43.638,	50.724,	25.121,	22.585,	20.916,	43.962,	26.351,	26.124,	26.715,	37.548,	42.759,	25.478,	40.93,	37.175,	42.899,	32.95,	65.318)
Mean_WG<-mean(WithinGroups)
SD_WG<-sd(WithinGroups)
DenaPlot<-data.frame(Mean_Group=c(Mean_BG,Mean_WG),sd_Group=c(SD_BG,SD_WG),group=c("Between","Within"))
DenaPlot$lower<-DenaPlot$Mean_Group-DenaPlot$sd_Group
DenaPlot$upper<-DenaPlot$Mean_Group+DenaPlot$sd_Group
library(ggplot2)
ggplot(data = DenaPlot,mapping=aes(x=group,y=Mean_Group))+
  geom_point(size=4)+
  geom_errorbar(mapping=aes(ymin=lower,ymax=upper),width=.05)+
  xlab("group(s)")+
  ylab("LSI")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
