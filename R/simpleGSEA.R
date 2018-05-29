get_sd <- function(x){
  sqrt(rowSums((x - Matrix::rowMeans(x))^2)/(dim(x)[2] - 1))
}

get_s2n_ranking <- function(x, data, index, phen, permute = 0){
  if ((x-101) %% 100  == 0){
    if ((x+99) <= permute){
      print(paste0("Permuted phenotypes.......permutations: ", x, "-", x+99))
    }else{
      print(paste0("Permuted phenotypes.......permutations: ", x, "-", permute))
    }

  }

  if (permute == 0){
    label1 <- phen
    print(phen)
  }else{
    label1 <- sample(index, size = length(phen))
  }

  print(class(data[,label1]))
  print(dim(data[,label1]))

  label2 <- index[!index %in% label1]

  mean1 <- Matrix::rowMeans(data[,label1])
  mean2 <- Matrix::rowMeans(data[,label2])

  sd1 <- get_sd(data[,label1])
  sd2 <- get_sd(data[,label2])
  sd1 <- ifelse(sd1 > 0.2*abs(mean1), sd1, 0.2*abs(mean1))
  sd1 <- ifelse(sd1 == 0, 0.2, sd1)
  sd2 <- ifelse(sd2 > 0.2*abs(mean2), sd2, 0.2*abs(mean2))
  sd2 <- ifelse(sd2 == 0, 0.2, sd2)

  s2n <- (mean1 - mean2)/(sd1 + sd2)
}


#' Do a simplified version of GSEA with most default parameters
#'
#' `simple_gsea` is used to output result of GSEA. You just need to import expression matrix and gene set as element files.
#' Then you have to provide the group info of samples. See this examples for more information
#'
#'
#' @param dataexpr Matrix or dataframe of expression data (rows are genes and columns are samples)
#' @param group A vector which distinguishs some samples to others(Just support two phenotypes)
#' @param geneset Gene set databases in Molecular Signatures Database (MSigDB)
#' @param gsminsize Minimum size (in genes) for database gene sets to be considered (default: 15)
#' @param gsmaxsize Maximum size (in genes) for database gene sets to be considered (default: 500)
#' @param nperm Number of random permutations (default: 1000)
#'
#' @return gseaResult object
#' @export
#'
#' @examples
#' data <- read.table(file = "expression.txt", sep = "\t", header = T, row.names = 1, stringsAsFactors = F, quote = "")
#' group <- c(rep("m", 15), rep("f", 17))
#' class <- factor(group)
#' gene_set <- readLines("C1.gmt")
#' res <- simple_gsea(dataexpr = data, group = class, geneset = gene_set, gsminsize = 15, gsmaxsize = 500, nperm = 1000)
#'
simple_gsea <- function(dataexpr, group, geneset, gsminsize=15, gsmaxsize=500, nperm=1000){
  design <- model.matrix(~group)
  colnames(design) <- levels(group)
  rownames(design) <- colnames(dataexpr)

  phen1 <- colnames(design)[1]
  phen2 <- colnames(design)[2]

  # Filter gene sets
  # gsminsize = 15
  # gsmaxsize = 500
  set_list <- list()
  set_name <- c()
  for(line in geneset){
    tmp <- strsplit(line, "\t")[[1]]
    gset_name <- tmp[1]
    gset_id <- tmp[-c(1,2)]
    gset_id <- gset_id[gset_id %in% rownames(dataexpr)]
    if (length(gset_id) >= gsminsize & length(gset_id) <= gsmaxsize){
      set_list <- c(set_list, list(gset_id))
      set_name <- c(set_name, gset_name)
    }
  }
  names(set_list) <- set_name

  N <- nrow(dataexpr)
  Ng <- length(set_list)
  Ns <- ncol(dataexpr)
  max.Ng <- length(geneset)
  set_size <- unlist(lapply(set_list, length))

  print(c("Number of genes:", N))
  print(c("Number of Gene Sets:", Ng))
  print(c("Number of samples:", Ns))
  print(c("Original number of Gene Sets:", max.Ng))

  dataexpr <- dataexpr + 0.00000001

  sample_number <- nrow(design)
  obs.phen1_index <- as.numeric(which((design[,1]-design[,2]) == 1))
  # obs.phen2_index <- as.numeric(which((design[,1]-design[,2]) == 0))

  # nperm <- 1000

  sigma.correction <- "GeneCluster"

  data <- Matrix(as.matrix(dataexpr))

  obs.s2n <- as.numeric(get_s2n_ranking(x = 0, data = data, index = 1:sample_number, phen = obs.phen1_index, permute = 0))
  obs.s2n_index <- order(obs.s2n, decreasing = T)
  obs.s2n <- obs.s2n[obs.s2n_index]

  s2n_unsorted <- sapply(1:nperm, get_s2n_ranking, data = data, index = 1:sample_number, phen = obs.phen1_index, permute = nperm, simplify = "matrix")
  s2n_index <- apply(s2n_unsorted, 2, function(x){
    order(x, decreasing = T)
  })

  # compute observed ES
  obs.ES <- c()
  obs.ES_index <- c()
  obs.RES <- matrix(nrow= Ng, ncol= N)
  obs.indicator <- matrix(nrow = Ng, ncol = N)

  for (i in 1:Ng){
    set <- set_list[[i]]
    obs.gene_label <- row.names(dataexpr)[obs.s2n_index]
    obs.tag.indicator <- sign(obs.gene_label %in% set)
    obs.indicator[i,] <- obs.tag.indicator
    obs.no.tag.indicator <- 1 - obs.tag.indicator

    obs.correl.vector <- abs(obs.s2n)
    obs.sum.correl.tag <- sum(obs.correl.vector[obs.tag.indicator == 1])
    obs.ES_cumsum <- cumsum(obs.tag.indicator*obs.correl.vector/obs.sum.correl.tag - obs.no.tag.indicator/(length(obs.gene_label) - length(set)))
    obs.RES[i,] <- obs.ES_cumsum
    obs.max.ES <- max(obs.ES_cumsum)
    obs.min.ES <- min(obs.ES_cumsum)
    if (obs.max.ES > -obs.min.ES){
      obs.ES[i] <- signif(obs.max.ES, digits = 5)
      obs.ES_index[i] <- which.max(obs.ES_cumsum)
    }else{
      obs.ES[i] <- signif(obs.min.ES, digits = 5)
      obs.ES_index[i] <- which.min(obs.ES_cumsum)
    }
  }

  # compute permutations ES
  ES <- matrix(0, nrow = Ng, ncol = nperm)
  for (i in 1:Ng){
    set <- set_list[[i]]
    print(paste("Computing random permutations' enrichment for gene set:", i, names(set_list)[i], sep=" "))

    for (r in 1:nperm){
      s2n_col_sorted <- s2n_unsorted[,r][s2n_index[,r]]
      gene_label <- row.names(dataexpr)[s2n_index[,r]]
      tag.indicator <- sign(gene_label %in% set)
      no.tag.indicator <- 1 - tag.indicator

      correl.vector <- abs(s2n_col_sorted)
      sum.correl.tag <- sum(correl.vector[tag.indicator == 1])
      ES_cumsum <- cumsum(tag.indicator*correl.vector/sum.correl.tag - no.tag.indicator/(length(gene_label) - length(set)))
      max.ES <- max(ES_cumsum)
      min.ES <- min(ES_cumsum)
      ES[i,r] <- signif(ifelse(max.ES > -min.ES, max.ES, min.ES), digits = 5)
    }
  }

  # Compute 3 types of p-values

  # Find nominal p-values
  print("Computing nominal p-values...")

  nom_pval <- c()
  for (i in 1:Ng){
    pos.ES <- c()
    neg.ES <- c()
    for(j in 1:nperm){
      if (ES[i,j] >= 0){
        pos.ES <- c(pos.ES, ES[i,j])
      }else{
        neg.ES <- c(neg.ES, ES[i,j])
      }
    }
    if (obs.ES[i] >= 0){
      nom_pval[i] <- signif(sum(pos.ES >= obs.ES[i])/length(pos.ES), digits = 5)
    }else{
      nom_pval[i] <- signif(sum(neg.ES <= obs.ES[i])/length(neg.ES), digits = 5)
    }
  }

  # Rescaling normalization for each gene set null
  print("Computing rescaling normalization for each gene set null...")

  ES_norm <- matrix(0, nrow = Ng, ncol = nperm)
  obs.ES_norm <- c()
  for (i in 1:Ng){
    pos.ES <- c()
    neg.ES <- c()
    for(j in 1:nperm){
      if (ES[i,j] >= 0){
        pos.ES <- c(pos.ES, ES[i,j])
      }else{
        neg.ES <- c(neg.ES, ES[i,j])
      }
    }

    pos.m <- mean(pos.ES)
    neg.m <- mean(abs(neg.ES))
    pos.ES <- pos.ES/pos.m
    neg.ES <- neg.ES/neg.m

    for (j in 1:nperm){
      if (ES[i,j] >= 0){
        ES_norm[i,j] <- ES[i,j]/pos.m
      }else{
        ES_norm[i,j] <- ES[i,j]/neg.m
      }
    }

    if (obs.ES[i] > 0){
      obs.ES_norm[i] <- obs.ES[i]/pos.m
    }else{
      obs.ES_norm[i] <- obs.ES[i]/neg.m
    }
  }

  # Compute FWER p-vals
  print("Computing FWER p-values...")

  max.ES_set <- c()
  min.ES_set <- c()
  for (j in 1:nperm){
    pos.ES <- c()
    neg.ES <- c()
    for (i in 1:Ng){
      if (ES_norm[i,j] >= 0){
        pos.ES <- c(pos.ES, ES_norm[i,j])
      }else{
        neg.ES <- c(neg.ES, ES_norm[i,j])
      }
    }
    if (length(pos.ES) > 0){
      max.ES_set <- c(max.ES_set, max(pos.ES))
    }
    if (length(neg.ES) > 0){
      min.ES_set <- c(min.ES_set, min(neg.ES))
    }
  }

  fwer_pval <- c()
  for (i in 1:Ng){
    if (obs.ES_norm[i] >= 0){
      fwer_pval[i] <- signif(sum(max.ES_set >= obs.ES_norm[i])/length(max.ES_set), digits = 5)
    }else{
      fwer_pval[i] <- signif(sum(min.ES_set <= obs.ES_norm[i])/length(min.ES_set), digits = 5)
    }
  }

  # Compute FDRs
  print("Computing FDR q-values...")

  tmp.ES.index <- order(obs.ES_norm, decreasing = T)
  ori.index <- order(tmp.ES.index)
  obs.ES_norm.sorted <- obs.ES_norm[tmp.ES.index]
  # set_name.sorted <- set_name[tmp.ES.index]

  NES <- c()
  ES_norm.mean <- c()
  obs.ES_norm.mean <- c()
  FDR.mean <- c()
  for (k in 1:Ng){
    NES[k] <- obs.ES_norm.sorted[k]

    count.col <- c()
    obs.count.col <- c()
    for (i in 1:nperm){
      ES_norm_vec <- ES_norm[,i]
      # obs.ES_norm_vec <- obs.ES_norm[i]
      if (NES[k] >= 0){
        count.col.norm <- sum(ES_norm_vec >= 0)
        obs.count.col.norm <- sum(obs.ES_norm >= 0)
        count.col[i] <- ifelse(count.col.norm > 0, sum(ES_norm_vec >= NES[k])/count.col.norm, 0)
        obs.count.col[i] <- ifelse(obs.count.col.norm > 0, sum(obs.ES_norm >= NES[k])/obs.count.col.norm, 0)
      }else{
        count.col.norm <- sum(ES_norm_vec < 0)
        obs.count.col.norm <- sum(obs.ES_norm < 0)
        count.col[i] <- ifelse(count.col.norm > 0, sum(ES_norm_vec <= NES[k])/count.col.norm, 0)
        obs.count.col[i] <- ifelse(obs.count.col.norm > 0, sum(obs.ES_norm <= NES[k])/obs.count.col.norm, 0)
      }
    }

    ES_norm.mean[k] <- mean(count.col)
    obs.ES_norm.mean[k] <- mean(obs.count.col)
    FDR.mean[k] <- ifelse(ES_norm.mean[k]/obs.ES_norm.mean[k] < 1, ES_norm.mean[k]/obs.ES_norm.mean[k], 1)
  }
  FDR.mean.sorted <- FDR.mean[ori.index]

  # Produce results report
  print("Producing result tables and plots...")

  obs.ES <- signif(obs.ES, digits = 5)
  obs.ES_norm <- signif(obs.ES_norm, digits = 5)
  nom_pval <- signif(nom_pval, digits = 4)
  FDR.mean.sorted <- signif(FDR.mean.sorted, digits = 5)
  fwer_pval <- signif(fwer_pval, digits = 3)

  report <- data.frame(GS_NAME = set_name,
                       SIZE = set_size,
                       ES = obs.ES,
                       NES = obs.ES_norm,
                       `NOM p-val` = nom_pval,
                       `FDR q-val` = FDR.mean.sorted,
                       `FWER p-val` = fwer_pval,
                       `RANK AT MAX` = obs.ES_index)

  result <- list(reprot = report,
                 RES = obs.RES,
                 s2n = obs.s2n,
                 pickindex = obs.ES_index,
                 indicator = obs.indicator,
                 group = c(phen1, phen2)
  )
}


