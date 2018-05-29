plot_func <- function(gset, RES, ES, s2n, pick, ind, gset_size, N, group){
  min.RES <- min(RES)
  max.RES <- max(RES)
  max.RES <- ifelse(max.RES < 0.3, 0.3, max.RES)
  min.RES <- ifelse(min.RES > -0.3, -0.3, min.RES)

  delta <- (max.RES - min.RES)*0.50
  min.plot <- min.RES - 2*delta
  max.plot <- max.RES
  max.corr <- max(s2n)
  min.corr <- min(s2n)
  obs.correl.vector.norm <- (s2n - min.corr)/(max.corr - min.corr)*1.25*delta + min.plot
  zero.corr.line <- (- min.corr/(max.corr - min.corr))*1.25*delta + min.plot
  col <- ifelse(ES > 0, 2, 4)

  phen1 <- group[1]
  phen2 <- group[2]

  sub_title <- paste("Number of genes: ", N, " (in list), ", gset_size, " (in gene set)", sep = "", collapse="")
  main_title <- paste("Gene Set ", ":", gset)

  plot(x = 1:N, y = RES, main = main_title, sub = sub_title,
       xlab = "Gene List Index", ylab = "Running Enrichment Score (RES)",
       xlim=c(1, N), ylim=c(min.plot, max.plot), type = "l",
       lwd = 2, cex = 1, col = col)
  for (j in seq(1, N, 20)) {
    lines(c(j, j), c(zero.corr.line, obs.correl.vector.norm[j]), lwd = 1,
          cex = 1, col = colors()[12]) # shading of correlation plot
  }
  lines(c(1, N), c(0, 0), lwd = 1, lty = 2, cex = 1, col = 1) # zero RES line
  lines(c(pick, pick), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = col) # max enrichment vertical line
  for (j in 1:N) {
    if (ind[j] == 1) {
      lines(c(j, j), c(min.plot + 1.25*delta, min.plot + 1.75*delta), lwd = 1, lty = 1, cex = 1, col = 1)  # enrichment tags
    }
  }
  lines(1:N, obs.correl.vector.norm, type = "l", lwd = 1, cex = 1, col = 1)
  lines(c(1, N), c(zero.corr.line, zero.corr.line), lwd = 1, lty = 1, cex = 1, col = 1) # zero correlation horizontal line
  arg.correl <- which(abs(s2n) == min(abs(s2n)))
  lines(c(arg.correl, arg.correl), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = 3) # zero crossing correlation vertical line

  leg.txt <- paste0("  ", phen1, "  ")
  text(x=1, y=min.plot, adj = c(0, 0), labels=leg.txt, cex = 1.0)

  leg.txt <- paste0("  ", phen2, "  ")
  text(x=N, y=min.plot, adj = c(1, 0), labels=leg.txt, cex = 1.0)

  adjx <- ifelse(ES > 0, 0, 1)

  leg.txt <- paste("Peak at ", pick, sep="", collapse="")
  text(x=pick, y=min.plot + 1.8*delta, adj = c(adjx, 0), labels=leg.txt, cex = 1.0)

  leg.txt <- paste("Zero crossing at ", arg.correl, sep="", collapse="")
  text(x=arg.correl, y=min.plot + 1.95*delta, adj = c(adjx, 0), labels=leg.txt, cex = 1.0)
}

# gseaplot(expr = data_set, list = res, gset = "chr10q11")

#' Plot GSEA graph
#'
#' This GSEA graph is nearly the same with GSEA official plot with very minimal change to the source code.
#'
#' @param expr Matrix or dataframe of expression data (rows are genes and columns are samples)
#' @param list Result of simple_gsea function
#' @param gset You can choose one gene set as input or Null which means those satisfied gene sets will be plot
#' @param nompval Significance threshold for nominal p-vals for gene sets (default: 0.05)
#' @param fdr Significance threshold for FDR q-vals for gene sets (default: 0.25)
#'
#' @return GSEA graph
#' @export
#'
#' @examples
#' gseaplot <- function(expr = data, list = res, gset = NULL, nompval = 0.05, fdr = 0.25)
#' gseaplot <- function(expr = data, list = res, gset = "chr10q11")
#'
gseaplot <- function(expr, list, gset=NULL, nompval=0.05, fdr=0.25){
  obs.s2n <- list$s2n
  gene_num <- nrow(expr)
  group <- list$group

  if (length(gset) == 1){
    local <- match(gset, list$reprot$GS_NAME)
    obs.RES <- list$RES[local,]
    pick_index <- list$pickindex[local]
    indicator <- list$indicator[local,]
    obs.ES <- list$reprot$ES[local]
    gset_size <- list$reprot$SIZE[local]

    plot_func(gset = gset, RES = obs.RES, ES = obs.ES, s2n = obs.s2n,
              pick = pick_index, ind = indicator, gset_size = gset_size,
              N = gene_num, group = group)
  }else if (length(gset) > 1){
    stop("The number of gene set over the limit, only one gene set can be showing")
  }else{
    print(paste0("Pval < ", nompval, " & FDR < ", fdr, " geneset plot will be save in pwd !"))
    result <- list$reprot
    result <- subset(result, NOM.p.val <= nompval & FDR.q.val <= fdr)

    if (nrow(result) == 0){
      stop("No GSEA enrichment above the threshold !")
    }else{
      for (i in 1:nrow(result)){
        local <- match(result$GS_NAME[i], list$reprot$GS_NAME)
        obs.RES <- list$RES[local,]
        pick_index <- list$pickindex[local]
        indicator <- list$indicator[local,]
        obs.ES <- list$reprot$ES[local]
        gset_size <- list$reprot$SIZE[local]

        png(filename = paste0(result$GS_NAME[i], "_", local, ".png"), height = 15, width = 25, units = "cm", res = 300)

        plot_func(gset = result$GS_NAME[i], RES = obs.RES, ES = obs.ES, s2n = obs.s2n,
                  pick = pick_index, ind = indicator, gset_size = gset_size,
                  N = gene_num, group = group)

        dev.off()
      }
    }
  }
}
