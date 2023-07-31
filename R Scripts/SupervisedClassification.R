# Load required libraries
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


# Data Preparation --------------------------------------------------------
# Link to Raven selection tables
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


# Initialize an empty data frame
combined.US2.df <- data.frame()

# Loop over the length of the files
for (a in 1:length(files)) {
  # Read in the data from a file in the 'files' list
  temp.table <- read.delim2(files[a], stringsAsFactors = F)

  # Subset the data for only those rows where View is 'Spectrogram 1'
  temp.table <- subset(temp.table, View == "Spectrogram 1")

  # Order the data by the Begin Time in ascending order
  temp.table <-
    temp.table[order(as.numeric(temp.table$Begin.Time..s.)), ]

  # Split the file names to create group labels
  group <- str_split_fixed(short.files[a], pattern = 'Song', n = 2)[, 1]
  group.label <-
    str_split_fixed(short.files[a], pattern = '.txt', n = 2)[, 1]

  # Add the group labels to the data frame
  new.temp.table <- cbind.data.frame(temp.table, group, group.label)

  # Initialize an empty list for note intervals
  NoteIntervalList <- list()

  # Loop over the rows of the data frame to calculate inter-note intervals
  for (b in 1:(nrow(new.temp.table) - 1)) {
    FirstRow <- new.temp.table[b, ]
    SecondRow <- new.temp.table[b + 1, ]
    InterNoteInterval <-
      as.numeric(SecondRow$Begin.Time..s.) - as.numeric(FirstRow$End.Time..s.)
    NoteIntervalList[[b]] <- InterNoteInterval
  }

  # Add 'NA' to the end of the list
  NoteIntervalList[[nrow(new.temp.table)]] <- 'NA'

  # Add the Inter-Note Interval list as a column to the data frame
  new.temp.table$IOI <- unlist(NoteIntervalList)

  # Combine the new temporary table with the main combined data frame
  combined.US2.df <-
    rbind.data.frame(combined.US2.df, new.temp.table)
}

# Convert the IOI column to numeric
combined.US2.df$IOI <- as.numeric(combined.US2.df$IOI)

# Plot a histogram of the IOI column
hist(as.numeric(combined.US2.df$IOI))

# Subset data where IOI <= 5, and drop unused factor levels
#combined.US2.df <- droplevels(subset(combined.US2.df,IOI <= 5))

# Plot a histogram of the updated IOI column
hist(as.numeric(combined.US2.df$IOI))

# Display the range of the IOI column
range(as.numeric(combined.US2.df$IOI))

# Display a frequency table of the Sequence column
table(combined.US2.df$Sequence)

# Subset data where Sequence is 'us'
TempUS <- subset(combined.US2.df, Sequence = 'us')

# Plot a histogram of the IOI column for the subsetted data
hist(as.numeric(TempUS$IOI))

# Convert character columns in columns 4 to 13 to numeric
combined.US2.df[, 4:13] <-
  combined.US2.df[, 4:13] %>% mutate_if(is.character, as.numeric)


# Loop to estimate features for each phrase -------------------------------

pair.index <- unique(combined.US2.df$group.label)

siamangUS2df <- data.frame()
for (b in 1:length(pair.index)) {
  temp.table <- subset(combined.US2.df, group.label == pair.index[b])

  individual <- temp.table$group[1]
  note.dur <-
    mean(temp.table$End.Time..s. - temp.table$Begin.Time..s.)
  call.dur  <-
    max(temp.table$End.Time..s.) - min(temp.table$Begin.Time..s.)
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
  noterate <- nrow(temp.table) / call.dur
  note1dur <- temp.table[1, ]$Dur.90...s.

  note1maxfreq <- temp.table[1, ]$Freq.95...Hz.
  note2dur <- temp.table[2, ]$Dur.90...s.

  note2maxfreq <- temp.table[2, ]$Freq.95...Hz.
  range.bw <- maxbw - minbw
  rest.dur <- call.dur - sum(temp.table$Dur.90...s.)
  lastnotedur <- temp.table[nrow(temp.table), ]$Dur.90...s.
  lastnoteminfreq <-
    temp.table[nrow(temp.table), ]$Freq.5...Hz.
  lastnotemaxfreq <-
    temp.table[nrow(temp.table), ]$Freq.95...Hz.

  us2 <- temp.table[which(temp.table$Sequence == 'us'), ]

  if (nrow(us2)  > 0) {
    us2 <- us2[, c(10, 11, 13)]
    colnames(us2) <-
      c("US2Dur.90...s.", "US2Freq.95...Hz.", "US2BW.90...Hz.")

    n.bk <- length(which(temp.table$Sequence == 'bk'))
    n.bm <- length(which(temp.table$Sequence == 'bm'))

    temp.US2.df <-
      cbind.data.frame(
        individual,
        call.id,
        call.dur,
        nnotes,
        min95,
        minbw,
        maxbw,
        mean95,
        max95,
        meanbw,
        mindurnote,
        maxdurnote,
        noterate,
        note1dur,
        note1maxfreq,
        note2dur,
        note2maxfreq,
        range.bw,
        rest.dur,
        lastnotedur,
        lastnotemaxfreq,
        us2,
        n.bk,
        n.bm
      )

    temp.US2.df <- temp.US2.df[1, ]
    siamangUS2df <- rbind.data.frame(siamangUS2df, temp.US2.df)
  }
}

siamangUS2df$individual <- as.factor(siamangUS2df$individual)

table(siamangUS2df$individual)

siamangUS2df$sequence <-
  as.factor(str_split_fixed(siamangUS2df$call.id, pattern = '[.]', n = 2)[, 2])
siamangUS2df$song <-
  as.factor(str_split_fixed(siamangUS2df$call.id, pattern = 'US', n = 2)[, 1])

siamangUS2df$site <- recode(
  siamangUS2df$individual,
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


table(siamangUS2df$site)

siamangUS2df <-
  droplevels(subset(siamangUS2df, individual != 'IndvB'))

# Check number of rows
nrow(siamangUS2df)
sum(siamangUS2df$nnotes)


# Supervised Clustering ---------------------------------------------------

# Supervised clustering using spectrogram features by site
ml.model.specfeatures.svm <-
  e1071::svm(
    siamangUS2df[, -c(1:3, 18, 27, 28, 29)],
    siamangUS2df$individual,
    kernel = "radial",
    cross = nrow(siamangUS2df)
  )

levels(ml.model.specfeatures.svm$fitted) <-
  c(
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

levels(siamangUS2df$individual) <-
  c(
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

CaretConfusion <-
  caret::confusionMatrix(ml.model.specfeatures.svm$fitted, siamangUS2df$individual)
CaretConfusion$table

write.csv(CaretConfusion$table,
          'simangconfusionmatrix.csv')


ml.model.specfeatures.svm$tot.accuracy


ml.model.specfeatures.site <-
  e1071::svm(
    siamangUS2df[, -c(1:3, 18, 27, 28, 29)],
    siamangUS2df$site,
    kernel = "radial",
    cross = nrow(siamangUS2df)
  )


ml.model.specfeatures.site$tot.accuracy

ml.model.specfeatures.site$accuracies
mean(ml.model.specfeatures.site$accuracies)
sd(ml.model.specfeatures.site$accuracies)


ml.model.specfeatures.sequence <-
  e1071::svm(
    siamangUS2df[, -c(1:3, 18, 27, 28, 29)],
    siamangUS2df$sequence,
    kernel = "radial",
    cross = nrow(siamangUS2df)
  )


ml.model.specfeatures.sequence$tot.accuracy


# Table summarizing features ------------------------------------
table(siamangUS2df$individual)

min.max <-
  data.frame(
    mean = sapply(siamangUS2df[, -c(1:3, 18, 27, 28, 29)], mean),
    sd = sapply(siamangUS2df[, -c(1:3, 18, 27, 28, 29)], sd),
    min = sapply(siamangUS2df[, -c(1:3, 18, 27, 28, 29)], min),
    max = sapply(siamangUS2df[, -c(1:3, 18, 27, 28, 29)], max)
  )

min.max$se <- min.max$sd / sqrt(nrow(min.max))

min.max <- round(cbind.data.frame(min.max), 2)

min.max$mean.se <- paste(min.max$mean, '±', min.max$se)
min.max$range <- paste(min.max$min, '-', min.max$max, sep = '')

NewMinMax <-
  min.max[, c('mean.se', 'range')] # paste(min.max$mean.se , '\n', min.max$range )


Feature <- c(
  'Number of notes',
  'Minimum low frequency (Hz)',
  'Minimum bandwidth (Hz)',
  'Maximum bandwidth (Hz)',
  'Mean maximum frequency (Hz)',
  'Maximum high frequency (Hz)',
  'Mean bandwidth (Hz)',
  'Minimum note duration (s)',
  'Maximum note duration (s)',
  'Note rate (number of notes / duration)',
  'Note 1 duration (s)',
  'Note 1 maximum frequency (Hz)',
  'Note 2 duration (s)',
  'Note 2 maximum frequency (Hz)',
  'Rest duration (s)',
  'Last note duration (s)',
  'Last note maximum frequency (Hz)',
  'US2 Duration (s)',
  'US2 maximum frequency (Hz)',
  'US2 bandwidth (Hz)',
  'Number of bark notes',
  'Number of boom notes'
)

table.with.features <- cbind.data.frame(Feature, NewMinMax)
colnames(table.with.features) <- c('Feature', 'Mean ± SEM', 'Range')

myft <- flextable((table.with.features))
myft <- width(myft, width = 1)
myft <- bold(myft, part = "header")
myft

save_as_docx(myft, path = 'Table 2.Siamang Features Summary.docx')
