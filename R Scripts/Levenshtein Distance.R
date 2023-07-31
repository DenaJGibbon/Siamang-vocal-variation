# Prepare for LD
library(stringdist)
library(stringr)
library(dplyr)
library(ggpubr)

# List of the Raven selection tables
files <- list.files(
  'combine',
  recursive = T,
  full.names = T,
  pattern = '.txt'
)
short.files <- list.files(
  'combine',
  recursive = T,
  full.names = F,
  pattern = '.txt'
)


# Prepare the data --------------------------------------------------------
# Code to read in selection tables
combined.US2.df <- data.frame()
for (a in 1:length(files)) {
  tryCatch({
    print(a)
    temp.table <- read.delim2(files[a], stringsAsFactors = F)
    temp.table <- subset(temp.table, View == "Spectrogram 1")

    temp.table <-
      temp.table[order(as.numeric(temp.table$Begin.Time..s.)), ]

    group <- str_split_fixed(short.files[a], pattern = 'Song', n = 2)[, 1]
    group.label <-
      str_split_fixed(short.files[a], pattern = '.txt', n = 2)[, 1]
    new.temp.table <- cbind.data.frame(temp.table, group, group.label)

    NoteIntervalList <- list()
    for (b in 1:(nrow(new.temp.table) - 1)) {
      FirstRow <- new.temp.table[b, ]
      SecondRow <- new.temp.table[b + 1, ]
      InterNoteInterval <-
        as.numeric(SecondRow$Begin.Time..s.) - as.numeric(FirstRow$End.Time..s.)
      NoteIntervalList[[b]] <- InterNoteInterval
    }

    NoteIntervalList[[nrow(new.temp.table)]] <- 'NA'

    new.temp.table$IOI <- unlist(NoteIntervalList)

    combined.US2.df <-
      rbind.data.frame(combined.US2.df, new.temp.table)
  }, error = function(e) {
    cat("ERROR :", conditionMessage(e), "\n")
  })
}



combined.US2.df$IOI <- as.numeric(combined.US2.df$IOI)
hist(as.numeric(combined.US2.df$IOI))

combined.US2.df <- droplevels(subset(combined.US2.df, IOI <= 5))
hist(as.numeric(combined.US2.df$IOI))
range(as.numeric(combined.US2.df$IOI))

# solution
combined.US2.df[, 4:13] <-
  combined.US2.df[, 4:13] %>% mutate_if(is.character, as.numeric)

pair.index <- unique(combined.US2.df$group.label)

siamangLDdf <- data.frame()
for (b in 1:length(pair.index)) {
  temp.table <- subset(combined.US2.df, group.label == pair.index[b])
  group <- temp.table$group[1]
  nnotes <- nrow(temp.table)
  LDSequence <- paste(temp.table$Sequence, collapse = ' ')
  TempRow <- cbind.data.frame(group, LDSequence, nnotes)
  siamangLDdf <- rbind.data.frame(siamangLDdf, TempRow)
}

siamangLDdf <-
  droplevels(subset(siamangLDdf, group != "IndvB"))

siamangLDdf$site <- recode(
  siamangLDdf$group,
  Group1 = 'Sikundar',
  Group2 = 'Sikundar',
  Group3 = 'Kutapanjang',
  Group4 = 'Kutapanjang',
  Group5 = 'Kutapanjang',
  Group6 = 'Kutapanjang',
  Group7 = 'Kutapanjang',
  Group8 = 'Kutapanjang',
  Group9 = 'Soraya',
  Group10 = 'Soraya'
)


siamangLDdf$group <-
  factor(
    siamangLDdf$group,
    levels = c(
      "Group1",
      "Group2",
      "Group3",
      "Group4",
      "Group5",
      "Group6",
      "Group7",
      "Group8",
      "Group9",
      "Group10"
    )
  )

# LD on dataset -----------------------------------------------------------

id <- siamangLDdf$group
# Convert to matrix
matrix.output <- stringdistmatrix(siamangLDdf$LDSequence,
                                  method = 'lv')
matrix.output <- as.matrix(matrix.output, labels = TRUE)
row.names(matrix.output) <- id
colnames(matrix.output) <- id

# Start of code from Gamba's group
M <- matrix.output
#create the dimnames for the matrix
dimnames(M) <- list(NULL, NULL)
dimnames(M) <- list(id, id)

# this function creates a matrix containing the mean values of an original matrix
AverageMatValsFast <- function(mat) {
  uniRow <- unique(rownames(mat))
  uniCol <- unique(colnames(mat))
  newmat1 <- matrix(numeric(0), nrow = length(uniRow), ncol = ncol(mat))
  rownames(newmat1) <- uniRow
  colnames(newmat1) <- colnames(mat)

  newmat <-
    matrix(numeric(0),
           nrow = length(uniRow),
           ncol = length(uniRow))
  rownames(newmat) <- uniRow
  colnames(newmat) <- uniCol

  for (i in 1:nrow(newmat1)) {
    rowMatch <- which(rownames(mat) == uniRow[i])
    if (length(rowMatch) > 1) {
      newmat1[i, ] <- colMeans(mat[rowMatch, ])
    } else {
      newmat1[i, ] <- mat[rowMatch, ]
    }
  }

  for (j in 1:ncol(newmat)) {
    colMatch <- which(colnames(mat) == uniCol[j])

    if (length(rowMatch) > 1) {
      newmat[, j] <- round(rowMeans(newmat1[, colMatch]), 3)
    } else {
      newmat[, j] <- round(newmat1[, colMatch], 3)
    }
  }
  newmat
}

# calculate the mean for our matrix
SiamangLDGroupActual <- AverageMatValsFast(M)

SiamangLDGroupActualMean <- colMeans(SiamangLDGroupActual)
SiamangLDGroupActualSD <-
  sapply(as.data.frame(SiamangLDGroupActual), sd)
Group <- unique(id)
Type <- rep('Actual', length(Group))

ActualLDDF <- cbind.data.frame(SiamangLDGroupActualMean,
                               SiamangLDGroupActualSD, Group, Type)
# Create a corrplot
corrplot::corrplot (SiamangLDGroupActual,
                    is.corr = F,
                    cl.align.text = 'l')

# Create a corrplot
corrplot::corrplot (SiamangLDGroupActualSD,
                    is.corr = F,
                    cl.align.text = 'l')

# LD on random data -------------------------------------------------------
# Randomization
LDRandomDF <- data.frame()
for (a in 1:25) {
  set.seed(a)
  print(a)
  vocs.vec <- c('bm', 'bk', 'us')

  siamangLDdfrandom <- data.frame()
  for (b in 1:nrow(siamangLDdf)) {
    TempRow <- siamangLDdf[b, ]
    TempRow$LDSequence.ran <-
      paste(sample(vocs.vec, TempRow$nnotes, replace = T),
            collapse = ' ')
    siamangLDdfrandom <- rbind.data.frame(siamangLDdfrandom, TempRow)
  }

  id <- siamangLDdfrandom$group

  # Convert to matrix
  matrix.output.ran <-
    stringdistmatrix(siamangLDdfrandom$LDSequence.ran,
                     method = 'lv')
  matrix.output.ran <- as.matrix(matrix.output.ran, labels = TRUE)
  row.names(matrix.output.ran) <- id
  colnames(matrix.output.ran) <- id

  # Start of code from Gamba's group
  M <- matrix.output.ran
  #create the dimnames for the matrix
  dimnames(M) <- list(NULL, NULL)
  dimnames(M) <- list(id, id)


  # calculate the mean for our matrix
  SiamangLDGroupRandom <- AverageMatValsFast(M)

  MeanLD <- colMeans(SiamangLDGroupRandom)
  SDLD <- sapply(as.data.frame(SiamangLDGroupRandom), sd)
  Group <- attr(MeanLD, "names")
  TempDF <- cbind.data.frame(MeanLD, SDLD, Group)
  LDRandomDF <- rbind.data.frame(LDRandomDF, TempDF)

}


# Create a corrplot
corrplot::corrplot (SiamangLDGroupRandom, is.corr = F)

LDRandomDF$Type <- rep('Random', nrow(LDRandomDF))

colnames(ActualLDDF) <- colnames(LDRandomDF)

CombinedLDDF <- rbind.data.frame(LDRandomDF, ActualLDDF)
CombinedLDDF$Group <-
  factor(
    CombinedLDDF$Group,
    levels = c(
      "Group1",
      "Group2",
      "Group3",
      "Group4",
      "Group5",
      "Group6",
      "Group7",
      "Group8",
      "Group9",
      "Group10"
    )
  )


ggboxplot(
  data = CombinedLDDF,
  x = 'Group',
  y = 'MeanLD',
  color = 'Type',
  fill = 'Type',
  outlier.shape = NA
) +
  scale_color_manual(values = c('black', 'red')) +
  scale_fill_manual(values = c('black', 'red')) +
  theme(legend.title = element_blank()) + xlab('')
