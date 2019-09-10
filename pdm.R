## Poisson Dirichlet Model
## Yanshan Wang (email: wang.yanshan@mayo.edu)
##
## Partial codes are based on https://mollermara.com/blog/lda/

library(rjags)
library(plyr)
library(coda)
Rdir <- "path_to_workspace"
setwd(Rdir)

# load data
# dxcount: codes count, patient by CCS
# etable: population estimated data
load(paste(data_dir, "path_to_data", sep=""))

# check whether dxcount and etable are aligned
mtdm = t(dxcount)
E = etable

# number of topics
K = 20
words <- colnames(E)

genPDM <-
  function(mtdm,
           words,
           K,
           E,
           alpha.Words = 0.1,
           alpha.Topics = 0.1) {
    word <- do.call(rbind.fill.matrix,
                    lapply(1:ncol(mtdm), function(i)
                      t(rep(
                        1:length(mtdm[, i]), mtdm[, i]
                      ))))
    wordCount <- t(mtdm)
    N <- ncol(mtdm)                 #Number of documents
    Nwords <- length(words)         #Number of terms
    alphaTopics <- rep(alpha.Topics, K)      #Hyperprior on topics
    alphaWords <- rep(alpha.Words, Nwords)  #Hyperprior on words
    wordtopic <- matrix(NA, nrow(word), Nwords)
    ## Length of documents needed for indexing in jags
    doclengths <- rowSums(!is.na(word))
    ## How much we believe each document belongs to each of K topics
    topicdist <- matrix(NA, N, K)
    ## How much we believe each word belongs to each of K topics
    worddist <- matrix(NA, K, Nwords)
    ## Initialize gamma
    gamma <- t(rep(NA, N))
    ## All the parameters to be passed to jags
    dataList <- list(
      alphaTopics = alphaTopics,
      alphaWords = alphaWords,
      topicdist = topicdist,
      wordtopic = wordtopic,
      wordCount = wordCount,
      Ndocs = N,
      Ktopics = K,
      Nwords = Nwords,
      worddist = worddist,
      E = E,
      gamma = gamma
    )
    
    jags.model(
      'pdmmodel.jag',
      data = dataList,
      n.chains = 2,
      n.adapt = 100
    )
  }
## parameters
n.iter = 1000
n.burn = 500
thin = 1
jags <- genPDM(mtdm, words, K, E)
update(jags, n.burn)
## sampling
samples <-
  jags.samples(jags, c('worddist', 'topicdist', 'wordtopic'), n.iter, thin)
saveRDS(samples,
        paste(data_dir, "path_to_save",disease_id, sep=""))
