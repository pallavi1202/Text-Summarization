rm(list=ls()); cat("\014")
library(tm)
library(textmineR); library(igraph) 
library(ggplot2)
library(stringdist)
#install.packages("dendextend")
library(dendextend)

pdf.loc <- file.path("C:/Users/palla/OneDrive/Desktop/Movies")
myFiles <- normalizePath(list.files(path = pdf.loc, pattern = "pdf",  full.names = TRUE))
Docs.corpus <- Corpus(URISource(myFiles), readerControl = list(reader = readPDF(engine = "xpdf")))

#Most frequent Terms
temp.corpus <- tm_map(Docs.corpus, removeWords, words = stopwords())
temp.corpus <- tm_map(temp.corpus, removePunctuation)
dtm <- DocumentTermMatrix(temp.corpus)
term.freq <- colSums(as.matrix(dtm))
filtered <- term.freq[!(names(term.freq) %in% stopwords("SMART"))]
sorted <- sort(filtered[filtered >= 20], decreasing = TRUE)
tf_df = data.frame(word=names(sorted), frequency=sorted)

ggplot(tf_df, mapping = aes(word, frequency, fill=word)) +
  geom_col(show.legend = FALSE) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
  ggtitle("Term Frequency", subtitle = "Accross All TED Talks")

library(wordcloud)
palette <- brewer.pal(8,"Dark2")

set.seed(137)
wordcloud(words=names(sorted), 
          freq=sorted, 
          min.freq=3, 
          random.order=F,
          colors=palette)


# Extract content of each document in a format needed by textmineR to find text summary
corpora <- data.frame()
for (ff in 1:length(Docs.corpus)) {
  doc.content <- Docs.corpus[[ff]]$content
  temp <- paste(doc.content, collapse = ""); 
  temp1 <- data.frame(content=temp, stringsAsFactors = FALSE)
  corpora <- rbind(corpora, temp1)
}

corpora$id <- seq_along(corpora$content) # enumerate documents with id
corpora$name <- rownames(dtm)



# Find summary for each Doc. 
## Doc 1
embedding(corpora[1,], 10, 5) 

## Doc 2 
embedding(corpora[2,], 10, 15) 

## Doc 3 
embedding(corpora[3,], 10, 7) 

## Doc 4 
embedding(corpora[4,], 15, 9) 

## Doc 5 
embedding(corpora[5,], 8, 3) 

## Doc 6 
embedding(corpora[6,], 10, 5)

## Doc 7 
embedding(corpora[7,], 12, 7)

## Doc 8 
embedding(corpora[8,], 9, 5) 

## Doc 9 
embedding(corpora[9,], 10, 4) 

## Doc 10 
embedding(corpora[10,], 15, 8) 


# Perform hierarchical clustering on the documents

tfidfdtm <- DocumentTermMatrix(corpora, control = list(bounds = list(global = c(2, Inf)),weighting = function(x) weightTfIdf(x, normalize = FALSE)))
d <- dist(dtm, method="euclidian")
hc <- hclust(d, "ward.D")
dend <- as.dendrogram(hc)
labels(dend) <- labels.txt[order.dendrogram(dend)]
dend <- set(dend, "labels_cex", 0.6)
dend <- color_branches(dend, k = 6)
par(mar = c(3,0,2,6))
plot(dend,horiz = T,main = "Clustering Dendrogram of Ted Talks")






# Embedding & Summarization Functions

embedding <- function(corpora, Number.of.Sentences, Number.of.Topics)  { # corpora$id and corpora$content
  
  # 1) Create a TCM (term-co-occurrence matrix) with co-occurrence weights between the words, using skip grams, and a 5-word window
  # 
  tcm <- CreateTcm(doc_vec = corpora$content,
                   skipgram_window = 10, # 5-word window on each side
                   verbose = FALSE,
                   cpus = 2)
  
  
  # 2) Create embedding model. Assign probability to each word using Latent Dirichlet Allocation(LDA) and embed words in k topics. 
  # This function uses Gibbs sampling creating different embeding each run if you don't use set.seed!
  set.seed(123) # Keep this value to create simulated values that are reproducible.
  embeddings <- FitLdaModel(dtm = tcm,
                            k = Number.of.Topics, # Integer number of topics
                            iterations = 200,
                            burnin = 180,
                            alpha = 0.1,
                            beta = 0.05,
                            optimize_alpha = TRUE,
                            calc_likelihood = FALSE,
                            calc_coherence = FALSE,
                            calc_r2 = FALSE,
                            cpus = 2)
  
  doc <- corpora$content
  names(doc) <- corpora$id
  gamma = embeddings$gamma # Generate gamma - matrix of probability for each word (columns) to belong to one of the k topics (rows).
  sums <- summarizer(doc, gamma, Number.of.Sentences) # 3) Get Summary 
  
  return(sums)
}


summarizer <- function(doc, gamma, Number.of.Sentences) {
  
  if (length(doc) > 1 ) # Recursively manage to handle multiple docs at once
    return(sapply(doc, function(d) try(summarizer(d, gamma, Number.of.Sentences)))) # use a try statement to catch any problems that may arise
  
  # 1) Parse doc into sentences
  sentences <- stringi::stri_split_boundaries(doc, type = "sentence")[[ 1 ]]
  names(sentences) <- seq_along(sentences) # Enumerate Sentences
  
  # 2) Find Dtm (each sentence is a document – a row in Dtm), All possible words are the columns. 
  dtm <- CreateDtm(sentences, ngram_window = c(1,1), verbose = FALSE, cpus = 2) # Convert sentences to a document term matrix.
  # Rows - number of sentences, Columns - all the worsd
  
  dtm <- dtm[ rowSums(dtm) > 2 , ] # Remove any documents with 2 or fewer words
  vocab <- intersect(colnames(dtm), colnames(gamma)) # 
  dtm <- dtm / rowSums(dtm) # Scale word frequency between [0,1] 
  # Dtm X Gamma - Matrix product of each scaled word frequency (dtm) with the embedded values (gamma). 
  # Transforms the dtm matrix to number of sentences times the number of topics.
  dtm_topic <- dtm[ , vocab ] %*% t(gamma[ , vocab ]) 
  dtm_topic <- as.matrix(dtm_topic) # Each sentence (row) is represented by a vector of k topics.
  
  # 3) Get the pairwise distances between each document (embedded sentence)
  e_dist <- CalcHellingerDist(dtm_topic) # e_dist is a square matrix of size as the Num of documents
  
  e_measure <- (1 - e_dist) * 100 # Turn distance into a 0 to 100 measure
  diag(e_measure) <- 0 # set diagonal elements to zero (sentences connected to themselves)
  
  # Turn into a nearest-neighbor graph (Select "Closest.Neighbors" highest values)
  # Keep connections only to the top "Closest.Neighbors" most similar documents (sentences). This creates non-symmetric matrix. 
  Closest.Neighbors = 5
  e_measure1 <- apply(e_measure, 1, function(x){
    x[ x < sort(x, decreasing = TRUE)[ Closest.Neighbors ] ] <- 0 # Set to zero all elements except the max "Closest.Neighbors" values per ROW
    x
  })
  # NOTE:  e_measure1 becomes transposed!
  
  e_measure2 <- pmax(e_measure1, t(e_measure1)) # Take elementwise max, to turn e_measure into symmetric adjacency matrix. 
  e_measure3 <- graph.adjacency(e_measure2, mode = "undirected", weighted = TRUE) # Create a graph from adjacency matrix
  
  ev <- evcent(e_measure3) # Calculate eigenvector centrality
  
  # Format the result
  # Order ev.vector from largest to smallest and take top "Number.of.Sentences"
  result <- sentences[ names(ev$vector)[ order(ev$vector, decreasing = TRUE)[ 1:Number.of.Sentences ] ] ] 
  result <- result[ order(as.numeric(names(result))) ] # Order sentences by index of appearence
  result <- paste(result, collapse = " ")
  return(result)
}




