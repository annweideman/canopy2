#' Internal function
#'
#' Modify Canopy function to remove cutoff value for posterior probabilities.
#'
#' @keywords internal

canopy_post_v2 <- function(sampchain, projectname, K, numchain, burnin, thin,
                           optK, C = NULL, post.config.cutoff = NULL)
{
  if (is.null(C)) {
    C <- diag(nrow(sampchain[[1]][[1]][[1]]$cna))
    colnames(C) <- rownames(C) <- rownames(sampchain[[1]][[1]][[1]]$cna)
  }
  sampchaink <- sampchain[[which(K == optK)]]
  numchain <- length(sampchaink)
  samptreenew <- sampchaink[[1]][(burnin + 1):length(sampchaink[[1]])]
  numpostburn <- length(samptreenew)
  temp <- thin * c(1:(numpostburn/thin))
  samptreethin <- samptreenew[temp]
  length(samptreethin)
  for (numi in 2:numchain) {
    samptreenew <- sampchaink[[numi]][(burnin + 1):length(sampchaink[[numi]])]
    numpostburn <- length(samptreenew)
    temp <- thin * c(1:(numpostburn/thin))
    samptreethin <- c(samptreethin, samptreenew[temp])
  }
  samptreethin.lik <- rep(NA, length(samptreethin))
  for (treei in 1:length(samptreethin)) {
    samptreethin.lik[treei] <- samptreethin[[treei]]$likelihood
  }
  samptreethin <- samptreethin[which((rank(-1 * samptreethin.lik,
                                           ties.method = "first")) <= 5 * (length(samptreethin)/numchain))]
  samptreethin.lik <- rep(NA, length(samptreethin))
  for (treei in 1:length(samptreethin)) {
    samptreethin.lik[treei] <- samptreethin[[treei]]$likelihood
  }
  if (!is.null(sampchain[[1]][[1]][[1]]$cna)) {
    for (i in 1:length(samptreethin)) {
      samptreethin[[i]] <- Canopy::sortcna(samptreethin[[i]], C)
    }
  }
  for (i in 1:length(samptreethin)) {
    samptreethin[[i]]$clonalmut <- Canopy::getclonalcomposition(samptreethin[[i]])
  }
  config <- rep(NA, length(samptreethin))
  config[1] <- 1
  categ <- 1
  for (i in 2:length(samptreethin)) {
    for (categi in 1:categ) {
      list.a <- samptreethin[[i]]$clonalmut
      list.b <- samptreethin[[which(config == categi)[1]]]$clonalmut
      if ((sum(is.element(list.a, list.b)) == optK) & (sum(is.element(list.b,
                                                                      list.a)) == optK)) {
        config[i] <- categi
      }
    }
    if (is.na(config[i])) {
      config[i] <- categ + 1
      categ <- categ + 1
    }
  }
  z.temp <- (samptreethin.lik - mean(samptreethin.lik))/stats::sd(samptreethin.lik)
  samptreethin <- samptreethin[z.temp <= 1.5 & z.temp >= -1.5]
  samptreethin.lik <- samptreethin.lik[z.temp <= 1.5 & z.temp >=
                                         -1.5]
  config <- config[z.temp <= 1.5 & z.temp >= -1.5]
  config.summary <- matrix(nrow = length(unique(config)), ncol = 3)
  colnames(config.summary) <- c("Configuration", "Post_prob",
                                "Mean_post_lik")
  config.summary[, 1] <- unique(config)
  for (i in 1:nrow(config.summary)) {
    configi <- config.summary[i, 1]
    configi.temp <- which(config == configi)
    config.summary[i, 2] <- round(length(configi.temp)/length(config),
                                  3)
    config.summary[i, 3] <- round(max(samptreethin.lik[which(config ==
                                                               configi)]), 2)
  }

  for (treei in 1:length(samptreethin)) {
    output.tree <- samptreethin[[treei]]
    output.tree.Z <- output.tree$Z[, 2:ncol(output.tree$Z),
                                   drop = FALSE]
    output.tree.P <- apply(output.tree$P[2:nrow(output.tree$P),
                                         , drop = FALSE], 2, function(x) {
                                           x/sum(x)
                                         })
    output.tree$CCF <- output.tree.Z %*% output.tree.P
    samptreethin[[treei]] <- output.tree
  }
  return(list(samptreethin, samptreethin.lik, config, config.summary))
}
