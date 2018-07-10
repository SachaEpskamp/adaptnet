adaptnet <- function(
  network, # Estimated network model
  thresholds, # Estimated thresholds
  observed, # Vector with observed scores (NA for missing scores)
  dataBank, # Databank or string indicating where to load the databank. If missing, databank is generated
  dataBank_minN = 500,
  dataBank_N = 5000,
  ... # Arguments sent to IsingSampler
){
  # Check if network is a matrix:
  if (!is.matrix(network) || nrow(network) != ncol(network)){
    stop("'network' must be a symmetrical matrix")
  }

  # Number of variables:
  nVar <- ncol(network)

  # Check if thresholds is vector:
  if (!is.vector(thresholds) || !is.numeric(thresholds) || length(thresholds) != nVar){
    stop("'thresholds' must be a numeric vector")
  }

  # If observed is missing, create:
  if (missing(observed)){
    observed <- rep(NA, nVar)
  }

  # Currently known:
  curKnown <- which(!is.na(observed))

  # Currently unknown:
  curUnknown <- which(is.na(observed))

  # If the databank is missing, simulate a databank:
  if (missing(dataBank)){
    newCases <- matrix(observed,dataBank_N,length(observed),byrow=TRUE)
    dataBank <- as.data.frame(IsingSampler::IsingSampler(dataBank_N, network, thresholds, ..., constrain = newCases))
    names(dataBank) <- colnames(network)
  }

  # Subset the databank:
  subBank <- dataBank
  if (any(!is.na(observed))){
    for (i in which(!is.na(observed))){
      subBank <- subBank[subBank[,i] == observed[i],]
    }
  }

  # If the databank is too small, replenish:
  if (nrow(subBank) < dataBank_minN){
    newCases <- matrix(observed,dataBank_N,length(observed),byrow=TRUE)
    subBank <- as.data.frame(IsingSampler::IsingSampler(dataBank_N, network, thresholds, ..., constrain = newCases))
    names(subBank) <- names(dataBank)
  }

  # Compute entropy:
  Ent <- infotheo::entropy(subBank, method = "mm")

  # Compute conditional entropies:
  condEnt <- sapply(curUnknown,function(i){
    infotheo::condentropy(subBank[,-c(curKnown,i),drop=FALSE],subBank[,c(curKnown,i),drop=FALSE])
  })

  # REsults:
  Results <- list(
    mostInformative = curUnknown[which.min(condEnt)],
    entropy = Ent,
    condentropy = data.frame(
      variable = curUnknown,
      condentropy = condEnt
    ),
    dataBank = subBank
  )
  class(Results) <- c("adaptnet","list")
  return(Results)
}
