#' Probability Density Function for 3-parameter twopiece distributions
#' @param x: vector of quantiles.
#' @param p: vector of probabilities.
#' @param n: number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param  mu: location parameter.
#' @param par1: scale parameter 1.
#' @param par2: scale parameter 2.
#' @param FUN: a symmetric density f.
#' @param param: parameterizations used.
#' @param log, log.p: logical; if TRUE, probabilities p are given as log(p).
#' @return
#' @export

dtp3 <- function (x, mu, par1, par2, FUN, param = "tp", log = FALSE)
{
  param = match.arg(param, choices = c("tp", "eps", "isf"))
  if (!is.logical(log)) {
    stop("log.p must be a boolean")
  }
  # FIX (item 2): replaced ifelse()-based validation with proper if/stop guards
  if (param == "tp") {
    if (!(par1 > 0 & par2 > 0)) {
      stop("invalid arguments: par1 or/and par2 is/are not positive in the parametrization tp")
    }
    logPDF <- log(2) + ifelse(x < mu, FUN((x - mu)/par1, log = T),
                              FUN((x - mu)/par2, log = T)) - log(par1 + par2)
  }
  if (param == "eps") {
    sigma = par1
    gamma = par2
    if (!(sigma > 0 & abs(gamma) < 1)) {
      stop("invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps")
    }
    logPDF <- ifelse(x < mu, FUN((x - mu)/(sigma * (1 + gamma)), log = T),
                     FUN((x - mu)/(sigma * (1 - gamma)), log = T)) - log(sigma)
  }
  if (param == "isf") {
    sigma = par1
    gamma = par2
    if (!(sigma > 0 & gamma > 0)) {
      stop("invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf")
    }
    logPDF <- log(2) + ifelse(x < mu, FUN((x - mu)/(sigma * gamma), log = T),
                              FUN((x - mu)/(sigma/gamma), log = T)) - log(sigma * (gamma + 1/gamma))
  }
  if (log) return(logPDF) else return(exp(logPDF))
}

#' Probability Density Function for 4-parameter twopiece distributions
#' @param x: vector of quantiles.
#' @param p: vector of probabilities.
#' @param n: number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param mu: location parameter.
#' @param par1: scale parameter 1.
#' @param par2: scale parameter 2.
#' @param delta: shape parameter.
#' @param FUN: a symmetric density f.
#' @param param: parameterizations used.
#' @param log, log.p: logical; if TRUE, probabilities p are given as log(p).
#' @return
#' @export

dtp4 <- function (x, mu, par1, par2, delta, FUN, param = "tp", log = FALSE)
{
  param = match.arg(param, choices = c("tp", "eps", "isf"))
  # FIX (item 2): replaced ifelse()-based validation with proper if/stop guards
  if (param == "tp") {
    if (!(par1 > 0 & par2 > 0 & delta > 0)) {
      stop("invalid arguments: par1 or/and par2 or/and delta is/are no positive in the parametrization tp")
    }
    logPDF <- log(2) + ifelse(x < mu, FUN((x - mu)/par1, delta, log = T),
                              FUN((x - mu)/par2, delta, log = T)) - log(par1 + par2)
  }
  if (param == "eps") {
    sigma = par1
    gamma = par2
    if (!(sigma > 0 & abs(gamma) < 1 & delta > 0)) {
      stop("invalid arguments: par1 or/and delta is/are no positive or/and abs(par2) is no less that 1 in the parametrization eps")
    }
    logPDF <- ifelse(x < mu, FUN((x - mu)/(sigma * (1 + gamma)), delta, log = T),
                     FUN((x - mu)/(sigma * (1 - gamma)), delta, log = T)) - log(sigma)
  }
  if (param == "isf") {
    sigma = par1
    gamma = par2
    if (!(sigma > 0 & gamma > 0 & delta > 0)) {
      stop("invalid arguments: par1 or/and par2 or/and delta is/are no positive in the parametrization isf")
    }
    logPDF <- log(2) + ifelse(x < mu, FUN((x - mu)/(sigma * gamma), delta, log = T),
                              FUN((x - mu)/(sigma/gamma), delta, log = T)) - log(sigma * (gamma + 1/gamma))
  }
  if (log) return(logPDF) else return(exp(logPDF))
}



#' Cumulative Probability Function for 3-parameter twopiece distributions
#' @param x: vector of quantiles.
#' @param p: vector of probabilities.
#' @param n: number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param  mu: location parameter.
#' @param par1: scale parameter 1.
#' @param par2: scale parameter 2.
#' @param FUN: a symmetric density f.
#' @param param: parameterizations used.
#' @param log, log.p: logical; if TRUE, probabilities p are given as log(p).
#' @return
#' @export
ptp3 <- function (x, mu, par1, par2, FUN, param = "tp", log.p = FALSE)
{
  param = match.arg(param, choices = c("tp", "eps", "isf"))
  if (!is.logical(log.p)) {
    stop("log.p must be a boolean")
  }
  if (param == "tp") {
    if (!(par1 > 0 & par2 > 0)) {
      stop("invalid arguments: par1 or/and par2 is/are not positive in the parametrization tp")
    }
    CDF <- ifelse(x < mu, 2 * par1 * FUN((x - mu)/par1, log.p = F)/(par1 +
                                                                      par2), (par1 + par2 * (2 * FUN((x - mu)/par2, log.p = F) -
                                                                                               1))/(par1 + par2))
  }
  if (param == "eps") {
    sigma = par1
    gamma = par2
    if (!(sigma > 0 & abs(gamma) < 1)) {
      stop("invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps")
    }
    CDF <- ifelse(x < mu, (1 + gamma) * FUN((x - mu)/(sigma *
                                                        (1 + gamma)), log.p = F), gamma + (1 - gamma) * FUN((x -
                                                                                                               mu)/(sigma * (1 - gamma)), log.p = F))
  }
  # FIX (item 1): moved sigma/gamma assignment to BEFORE the validity check
  if (param == "isf") {
    sigma = par1
    gamma = par2
    if (!(sigma > 0 & gamma > 0)) {
      stop("invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf")
    }
    CDF <- ifelse(x < mu, 2 * gamma^2 * FUN((x - mu)/(sigma *
                                                        gamma), log.p = F)/(1 + gamma^2), (gamma^2 - 1 +
                                                                                             2 * FUN((x - mu)/(sigma/gamma), log.p = F))/(1 +
                                                                                                                                            gamma^2))
  }
  ifelse(log.p, return(log(CDF)), return(CDF))
}

#' Cumulative Probability Function for 4-parameter twopiece distributions
#' @param x: vector of quantiles.
#' @param p: vector of probabilities.
#' @param n: number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param mu: location parameter.
#' @param par1: scale parameter 1.
#' @param par2: scale parameter 2.
#' @param delta: shape parameter.
#' @param FUN: a symmetric density f.
#' @param param: parameterizations used.
#' @param log, log.p: logical; if TRUE, probabilities p are given as log(p).
#' @return
#' @export
ptp4 <- function (x, mu, par1, par2, delta, FUN, param = "tp", log.p = FALSE)
{
  param = match.arg(param, choices = c("tp", "eps", "isf"))
  if (!is.logical(log.p)) {
    stop("log.p must be a boolean")
  }
  # FIX (item 2): replaced ifelse()-based validation with proper if/stop guards
  if (param == "tp") {
    if (!(par1 > 0 & par2 > 0)) {
      stop("invalid arguments: par1 or/and par2 is/are not positive in the parametrization tp")
    }
    CDF <- ifelse(x < mu, 2 * par1 * FUN((x - mu)/par1, delta, log.p = F)/(par1 + par2),
                  (par1 + par2 * (2 * FUN((x - mu)/par2, delta, log.p = F) - 1))/(par1 + par2))
  }
  if (param == "eps") {
    sigma = par1
    gamma = par2
    if (!(sigma > 0 & abs(gamma) < 1)) {
      stop("invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps")
    }
    CDF <- ifelse(x < mu, (1 + gamma) * FUN((x - mu)/(sigma * (1 + gamma)), delta, log.p = F),
                  gamma + (1 - gamma) * FUN((x - mu)/(sigma * (1 - gamma)), delta, log.p = F))
  }
  if (param == "isf") {
    sigma = par1
    gamma = par2
    if (!(sigma > 0 & gamma > 0)) {
      stop("invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf")
    }
    CDF <- ifelse(x < mu, 2 * gamma^2 * FUN((x - mu)/(sigma * gamma), delta, log.p = F)/(1 + gamma^2),
                  (gamma^2 - 1 + 2 * FUN((x - mu)/(sigma/gamma), delta, log.p = F))/(1 + gamma^2))
  }
  if (log.p) return(log(CDF)) else return(CDF)
}

#' Quantile Function for 3-parameter twopiece distributions
#' @param x: vector of quantiles.
#' @param p: vector of probabilities.
#' @param n: number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param  mu: location parameter.
#' @param par1: scale parameter 1.
#' @param par2: scale parameter 2.
#' @param FUN: a symmetric density f.
#' @param param: parameterizations used.
#' @param log, log.p: logical; if TRUE, probabilities p are given as log(p).
#' @return
#' @export
qtp3 <- function (p, mu, par1, par2, FUN, param = "tp")
{
  param = match.arg(param, choices = c("tp", "eps", "isf"))
  # FIX (item 2): replaced ifelse()-based validation with proper if/stop guards
  if (param == "tp") {
    if (!(par1 > 0 & par2 > 0)) {
      stop("invalid arguments: par1 or/and par2 is/are not positive in the parametrization tp")
    }
    Q <- ifelse(p < par1/(par1 + par2),
                mu + par1 * FUN(0.5 * p * (par1 + par2)/par1),
                mu + par2 * FUN(0.5 * ((par1 + par2) * (1 + p) - 2 * par1)/par2))
  }
  if (param == "eps") {
    sigma = par1
    gamma = par2
    if (!(sigma > 0 & abs(gamma) < 1)) {
      stop("invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps")
    }
    Q <- ifelse(p < 0.5 * (1 + gamma),
                mu + sigma * (1 + gamma) * FUN(p/(1 + gamma)),
                mu + sigma * (1 - gamma) * FUN((p - gamma)/(1 - gamma)))
  }
  if (param == "isf") {
    sigma = par1
    gamma = par2
    if (!(sigma > 0 & gamma > 0)) {
      stop("invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf")
    }
    Q <- ifelse(p < gamma^2/(1 + gamma^2),
                mu + sigma * gamma * FUN(0.5 * p * (1 + gamma^2)/gamma^2),
                mu + sigma * FUN(0.5 * (p * (1 + gamma^2) + 1 - gamma^2))/gamma)
  }
  return(Q)
}

#' Quantile Function for 4-parameter twopiece distributions
#' @param x: vector of quantiles.
#' @param p: vector of probabilities.
#' @param n: number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param mu: location parameter.
#' @param par1: scale parameter 1.
#' @param par2: scale parameter 2.
#' @param delta: shape parameter.
#' @param FUN: a symmetric density f.
#' @param param: parameterizations used.
#' @param log, log.p: logical; if TRUE, probabilities p are given as log(p).
#' @return
#' @export
qtp4 <- function (p, mu, par1, par2, delta, FUN, param = "tp")
{
  param = match.arg(param, choices = c("tp", "eps", "isf"))
  # FIX (item 2): replaced ifelse()-based validation with proper if/stop guards
  if (param == "tp") {
    if (!(par1 > 0 & par2 > 0)) {
      stop("invalid arguments: par1 or/and par2 is/are not positive in the parametrization tp")
    }
    Q <- ifelse(p < par1/(par1 + par2),
                mu + par1 * FUN(0.5 * p * (par1 + par2)/par1, delta),
                mu + par2 * FUN(0.5 * ((par1 + par2) * (1 + p) - 2 * par1)/par2, delta))
  }
  if (param == "eps") {
    sigma = par1
    gamma = par2
    if (!(sigma > 0 & abs(gamma) < 1)) {
      stop("invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps")
    }
    Q <- ifelse(p < 0.5 * (1 + gamma),
                mu + sigma * (1 + gamma) * FUN(p/(1 + gamma), delta),
                mu + sigma * (1 - gamma) * FUN((p - gamma)/(1 - gamma), delta))
  }
  if (param == "isf") {
    sigma = par1
    gamma = par2
    if (!(sigma > 0 & gamma > 0)) {
      stop("invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf")
    }
    Q <- ifelse(p < gamma^2/(1 + gamma^2),
                mu + sigma * gamma * FUN(0.5 * p * (1 + gamma^2)/gamma^2, delta),
                mu + sigma * FUN(0.5 * (p * (1 + gamma^2) + 1 - gamma^2), delta)/gamma)
  }
  return(Q)
}

#' Random Number Generation Function for 3-parameter twopiece distributions
#' @param x: vector of quantiles.
#' @param p: vector of probabilities.
#' @param n: number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param  mu: location parameter.
#' @param par1: scale parameter 1.
#' @param par2: scale parameter 2.
#' @param FUN: a symmetric density f.
#' @param param: parameterizations used.
#' @param log, log.p: logical; if TRUE, probabilities p are given as log(p).
#' @return
#' @export
rtp3 <- function (n, mu, par1, par2, FUN, param = "tp")
{
  param = match.arg(param, choices = c("tp", "eps", "isf"))
  # FIX (item 3): generate a single shared pool of n draws so each observation
  # comes from one call to FUN, not two independent calls evaluated by ifelse.
  if (param == "tp") {
    if (!(par1 > 0 & par2 > 0)) {
      stop("invalid arguments: par1 or/and par2 is/are not positive in the parametrization tp")
    }
    u <- runif(n)
    z <- abs(FUN(n))
    sample <- ifelse(u < par1/(par1 + par2), mu - par1 * z, mu + par2 * z)
  }
  if (param == "eps") {
    sigma = par1
    gamma = par2
    if (!(sigma > 0 & abs(gamma) < 1)) {
      stop("invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps")
    }
    u <- runif(n)
    z <- abs(FUN(n))
    sample <- ifelse(u < 0.5 * (1 + gamma), mu - sigma * (1 + gamma) * z,
                     mu + sigma * (1 - gamma) * z)
  }
  if (param == "isf") {
    sigma = par1
    gamma = par2
    if (!(sigma > 0 & gamma > 0)) {
      stop("invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf")
    }
    u <- runif(n)
    z <- abs(FUN(n))
    sample <- ifelse(u < gamma^2/(1 + gamma^2), mu - sigma * gamma * z,
                     mu + sigma * z/gamma)
  }
  return(sample)
}

#' Random Number Generation Function for 4-parameter twopiece distributions
#' @param x: vector of quantiles.
#' @param p: vector of probabilities.
#' @param n: number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param mu: location parameter.
#' @param par1: scale parameter 1.
#' @param par2: scale parameter 2.
#' @param delta: shape parameter.
#' @param FUN: a symmetric density f.
#' @param param: parameterizations used.
#' @param log, log.p: logical; if TRUE, probabilities p are given as log(p).
#' @return
#' @export
rtp4 <- function (n, mu, par1, par2, delta, FUN, param = "tp")
{
  param = match.arg(param, choices = c("tp", "eps", "isf"))
  # FIX (item 3): generate a single shared pool of n draws so each observation
  # comes from one call to FUN, not two independent calls evaluated by ifelse.
  if (param == "tp") {
    if (!(par1 > 0 & par2 > 0 & delta > 0)) {
      stop("invalid arguments: par1 or/and par2 is/are not positive in the parametrization tp")
    }
    u <- runif(n)
    z <- abs(FUN(n, delta))
    sample <- ifelse(u < par1/(par1 + par2), mu - par1 * z, mu + par2 * z)
  }
  if (param == "eps") {
    sigma = par1
    gamma = par2
    if (!(sigma > 0 & abs(gamma) < 1 & delta > 0)) {
      stop("invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps")
    }
    u <- runif(n)
    z <- abs(FUN(n, delta))
    sample <- ifelse(u < 0.5 * (1 + gamma), mu - sigma * (1 + gamma) * z,
                     mu + sigma * (1 - gamma) * z)
  }
  if (param == "isf") {
    sigma = par1
    gamma = par2
    if (!(sigma > 0 & gamma > 0 & delta > 0)) {
      stop("invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf")
    }
    u <- runif(n)
    z <- abs(FUN(n, delta))
    sample <- ifelse(u < gamma^2/(1 + gamma^2), mu - sigma * gamma * z,
                     mu + sigma * z/gamma)
  }
  return(sample)
}
