# Probability Density Functions

dtp3 <- function (x, mu, par1, par2, FUN, param = "tp", log = FALSE) 
{
  param = match.arg(param, choices = c("tp", "eps", "isf"))
  if (!is.logical(log)) {
    stop("log.p must be a boolean")
  }
  if (param == "tp") {
    ifelse(par1 > 0 & par2 > 0, logPDF <- log(2) + ifelse(x < 
                                                            mu, FUN((x - mu)/par1, log = T), FUN((x - mu)/par2, 
                                                                                                 log = T)) - log(par1 + par2), logPDF <- "invalid arguments: par1 or/and par2 is/are not positive in the parametrization tp")
  }
  if (param == "eps") {
    sigma = par1
    gamma = par2
    ifelse(sigma > 0 & abs(gamma) < 1, logPDF <- ifelse(x < 
                                                          mu, FUN((x - mu)/(sigma * (1 + gamma)), log = T), 
                                                        FUN((x - mu)/(sigma * (1 - gamma)), log = T)) - log(sigma), 
           logPDF <- "invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps")
  }
  if (param == "isf") {
    sigma = par1
    gamma = par2
    ifelse(sigma > 0 & gamma > 0, logPDF <- log(2) + ifelse(x < 
                                                              mu, FUN((x - mu)/(sigma * gamma), log = T), FUN((x - 
                                                                                                                 mu)/(sigma/gamma), log = T)) - log(sigma * (gamma + 
                                                                                                                                                               1/gamma)), logPDF <- "invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf")
  }
  ifelse(is.numeric(logPDF), ifelse(log, return(logPDF), return(exp(logPDF))), 
         logPDF)
}


dtp4 <- function (x, mu, par1, par2, delta, FUN, param = "tp", log = FALSE) 
{
  param = match.arg(param, choices = c("tp", "eps", "isf"))
  if (param == "tp") {
    ifelse(par1 > 0 & par2 > 0 & delta > 0, logPDF <- log(2) + 
             ifelse(x < mu, FUN((x - mu)/par1, delta, log = T), 
                    FUN((x - mu)/par2, delta, log = T)) - log(par1 + 
                                                                par2), logPDF <- "invalid arguments: par1 or/and par2 or/and delta is/are no positive in the parametrization tp")
  }
  if (param == "eps") {
    sigma = par1
    gamma = par2
    ifelse(sigma > 0 & abs(gamma) < 1 & delta > 0, logPDF <- ifelse(x < 
                                                                      mu, FUN((x - mu)/(sigma * (1 + gamma)), delta, log = T), 
                                                                    FUN((x - mu)/(sigma * (1 - gamma)), delta, log = T)) - 
             log(sigma), logPDF <- "invalid arguments: par1 or/and delta is/are no positive or/and abs(par2) is no less that 1 in the parametrization eps")
  }
  if (param == "isf") {
    sigma = par1
    gamma = par2
    ifelse(sigma > 0 & gamma > 0 & delta > 0, logPDF <- log(2) + 
             ifelse(x < mu, FUN((x - mu)/(sigma * gamma), delta, 
                                log = T), FUN((x - mu)/(sigma/gamma), delta, 
                                              log = T)) - log(sigma * (gamma + 1/gamma)), logPDF <- "invalid arguments: par1 or/and par2 or/and delta is/are no positive in the parametrization isf")
  }
  ifelse(is.numeric(logPDF), ifelse(log, return(logPDF), return(exp(logPDF))), 
         logPDF)
}



# Cumulative Probability Functions

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
  if (param == "isf") {
    if (!(sigma > 0 & gamma > 0)) {
      stop("invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf")
    }
    sigma = par1
    gamma = par2
    CDF <- ifelse(x < mu, 2 * gamma^2 * FUN((x - mu)/(sigma * 
                                                        gamma), log.p = F)/(1 + gamma^2), (gamma^2 - 1 + 
                                                                                             2 * FUN((x - mu)/(sigma/gamma), log.p = F))/(1 + 
                                                                                                                                            gamma^2))
  }
  ifelse(log.p, return(log(CDF)), return(CDF))
}


ptp4 <- function (x, mu, par1, par2, delta, FUN, param = "tp", log.p = FALSE) 
{
  param = match.arg(param, choices = c("tp", "eps", "isf"))
  if (!is.logical(log.p)) {
    stop("log.p must be a boolean")
  }
  if (param == "tp") {
    ifelse(par1 > 0 & par2 > 0, CDF <- ifelse(x < mu, 2 * 
                                                par1 * FUN((x - mu)/par1, delta, log.p = F)/(par1 + 
                                                                                               par2), (par1 + par2 * (2 * FUN((x - mu)/par2, delta, 
                                                                                                                              log.p = F) - 1))/(par1 + par2)), CDF <- "invalid arguments: par1 or/and par2 is/are not positive in the parametrization tp")
  }
  if (param == "eps") {
    sigma = par1
    gamma = par2
    ifelse(sigma > 0 & abs(gamma) < 1, CDF <- ifelse(x < 
                                                       mu, (1 + gamma) * FUN((x - mu)/(sigma * (1 + gamma)), 
                                                                             delta, log.p = F), gamma + (1 - gamma) * FUN((x - 
                                                                                                                             mu)/(sigma * (1 - gamma)), delta, log.p = F)), CDF <- "invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps")
  }
  if (param == "isf") {
    sigma = par1
    gamma = par2
    ifelse(sigma > 0 & gamma > 0, CDF <- ifelse(x < mu, 2 * 
                                                  gamma^2 * FUN((x - mu)/(sigma * gamma), delta, log.p = F)/(1 + 
                                                                                                               gamma^2), (gamma^2 - 1 + 2 * FUN((x - mu)/(sigma/gamma), 
                                                                                                                                                delta, log.p = F))/(1 + gamma^2)), CDF <- "invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf")
  }
  ifelse(is.numeric(CDF), ifelse(log.p, return(log(CDF)), return(CDF)), 
         CDF)
}

# Quantile Functions

qtp3 <- function (p, mu, par1, par2, FUN, param = "tp") 
{
  param = match.arg(param, choices = c("tp", "eps", "isf"))
  if (param == "tp") {
    ifelse(par1 > 0 & par2 > 0, Q <- ifelse(p < par1/(par1 + 
                                                        par2), mu + par1 * FUN(0.5 * p * (par1 + par2)/par1), 
                                            mu + par2 * FUN(0.5 * ((par1 + par2) * (1 + p) - 
                                                                     2 * par1)/par2)), Q <- "invalid arguments: par1 or/and par2 is/are not positive in the parametrization tp")
  }
  if (param == "eps") {
    sigma = par1
    gamma = par2
    ifelse(sigma > 0 & abs(gamma) < 1, Q <- ifelse(p < 0.5 * 
                                                     (1 + gamma), mu + sigma * (1 + gamma) * FUN(p/(1 + 
                                                                                                      gamma)), mu + sigma * (1 - gamma) * FUN((p - gamma)/(1 - 
                                                                                                                                                             gamma))), Q <- "invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps")
  }
  if (param == "isf") {
    sigma = par1
    gamma = par2
    ifelse(sigma > 0 & gamma > 0, Q <- ifelse(p < gamma^2/(1 + 
                                                             gamma^2), mu + sigma * gamma * FUN(0.5 * p * (1 + 
                                                                                                             gamma^2)/gamma^2), mu + sigma * FUN(0.5 * (p * (1 + 
                                                                                                                                                               gamma^2) + 1 - gamma^2))/gamma), Q <- "invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf")
  }
  return(Q)
}


qtp4 <- function (p, mu, par1, par2, delta, FUN, param = "tp") 
{
  param = match.arg(param, choices = c("tp", "eps", "isf"))
  if (param == "tp") {
    ifelse(par1 > 0 & par2 > 0, Q <- ifelse(p < par1/(par1 + 
                                                        par2), mu + par1 * FUN(0.5 * p * (par1 + par2)/par1, 
                                                                               delta), mu + par2 * FUN(0.5 * ((par1 + par2) * (1 + 
                                                                                                                                 p) - 2 * par1)/par2, delta)), Q <- "invalid arguments: par1 or/and par2 is/are not positive in the parametrization tp")
  }
  if (param == "eps") {
    sigma = par1
    gamma = par2
    ifelse(sigma > 0 & abs(gamma) < 1, Q <- ifelse(p < 0.5 * 
                                                     (1 + gamma), mu + sigma * (1 + gamma) * FUN(p/(1 + 
                                                                                                      gamma), delta), mu + sigma * (1 - gamma) * FUN((p - 
                                                                                                                                                        gamma)/(1 - gamma), delta)), Q <- "invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps")
  }
  if (param == "isf") {
    sigma = par1
    gamma = par2
    ifelse(sigma > 0 & gamma > 0, Q <- ifelse(p < gamma^2/(1 + 
                                                             gamma^2), mu + sigma * gamma * FUN(0.5 * p * (1 + 
                                                                                                             gamma^2)/gamma^2, delta), mu + sigma * FUN(0.5 * 
                                                                                                                                                          (p * (1 + gamma^2) + 1 - gamma^2), delta)/gamma), 
           Q <- "invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf")
  }
  return(Q)
}

# Random Number Generation Functions

rtp3 <- function (n, mu, par1, par2, FUN, param = "tp") 
{
  param = match.arg(param, choices = c("tp", "eps", "isf"))
  if (param == "tp") {
    ifelse(par1 > 0 & par2 > 0, sample <- ifelse(runif(n) < 
                                                   par1/(par1 + par2), mu - par1 * abs(FUN(n)), mu + 
                                                   par2 * abs(FUN(n))), sample <- "invalid arguments: par1 or/and par2 is/are not positive in the parametrization tp")
  }
  if (param == "eps") {
    sigma = par1
    gamma = par2
    ifelse(sigma > 0 & abs(gamma) < 1, sample <- ifelse(runif(n) < 
                                                          0.5 * (1 + gamma), mu - sigma * (1 + gamma) * abs(FUN(n)), 
                                                        mu + sigma * (1 - gamma) * abs(FUN(n))), sample <- "invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps")
  }
  if (param == "isf") {
    sigma = par1
    gamma = par2
    ifelse(sigma > 0 & gamma > 0, sample <- ifelse(runif(n) < 
                                                     gamma^2/(1 + gamma^2), mu - sigma * gamma * abs(FUN(n)), 
                                                   mu + sigma * abs(FUN(n))/gamma), sample <- "invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf")
  }
  return(sample)
}


rtp4 <- function (n, mu, par1, par2, delta, FUN, param = "tp") 
{
  param = match.arg(param, choices = c("tp", "eps", "isf"))
  if (param == "tp") {
    ifelse(par1 > 0 & par2 > 0 & delta > 0, sample <- ifelse(runif(n) < 
                                                               par1/(par1 + par2), mu - par1 * abs(FUN(n, delta)), 
                                                             mu + par2 * abs(FUN(n, delta))), sample <- "invalid arguments: par1 or/and par2 is/are not positive in the parametrization tp")
  }
  if (param == "eps") {
    sigma = par1
    gamma = par2
    ifelse(sigma > 0 & abs(gamma) < 1 & delta > 0, sample <- ifelse(runif(n) < 
                                                                      0.5 * (1 + gamma), mu - sigma * (1 + gamma) * abs(FUN(n, 
                                                                                                                            delta)), mu + sigma * (1 - gamma) * abs(FUN(n, delta))), 
           sample <- "invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps")
  }
  if (param == "isf") {
    sigma = par1
    gamma = par2
    ifelse(sigma > 0 & gamma > 0 & delta > 0, sample <- ifelse(runif(n) < 
                                                                 gamma^2/(1 + gamma^2), mu - sigma * gamma * abs(FUN(n, 
                                                                                                                     delta)), mu + sigma * abs(FUN(n, delta))/gamma), 
           sample <- "invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf")
  }
  return(sample)
}

