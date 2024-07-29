########## quasi-binary endpoints ########
UNITED.oc = function(ttox,
                    teff,
                    wt_tox,
                    wt_eff,
                    cohortsize = 3,
                    ncohort = 10,
                    ntrial = 5000,
                    phiT = 0.3,
                    muT = 0.9,
                    phiE = 0.5,
                    muE = 0.85,
                    alphaT=0.5,
                    betaT=0.5,
                    alphaE=0.5,
                    betaE=0.5) {
  library("Iso")
  
  ###find admissible set of toxicity
  adm_tox <- function(n, ytox) {
    nn = n[n != 0]
    yytox = ytox[which(n != 0)]
    at = alphaT + yytox
    bt = betaT + nn - yytox
    Tox_prob = 1 - pbeta(phiT, at, bt)
    AT_naive = which(Tox_prob <= muT)
    if (length(AT_naive)==0){
      AT=AT_naive
    }
    else{
      full_seq = seq(min(AT_naive),max(AT_naive),1)
      if (length(setdiff(full_seq, AT_naive)) == 0) {
        AT = AT_naive
      }
      else{
        AT = AT_naive[AT_naive < min(setdiff(full_seq, AT_naive))]
      }
      AT = 1:max(AT)
    }
    return(AT)
  }
  
  ###find admissible set of efficacy
  adm_eff <- function(n, yeff) {
    nn = n[n != 0]
    yyeff = yeff[which(n != 0)]
    ae = alphaE + yyeff
    be = betaE + nn - yyeff
    Eff_prob = pbeta(phiE, ae, be)
    AR_naive = which(Eff_prob < muE)
    if (length(AR_naive)==0){
      AR=AR_naive
    }
    else{
      full_seq = seq(min(AR_naive),max(AR_naive),1)
      if (length(setdiff(full_seq, AR_naive)) == 0) {
        AR = AR_naive
      } else {
        AR = AR_naive[AR_naive > max(setdiff(full_seq, AR_naive))]
      }
      AR = min(AR):length(nn)
    }
    return(AR)
  }
  
  
  ###AIC model selection
  
  AIC <- function(n, yeff) {
    aic = rep(0, length(n[n != 0]))
    for (l in 1:length(n[n != 0])) {
      if (l == 1) {
        jstar = length(n[n != 0])
        sm = sum(yeff[1:jstar])
        nm = sum(n[1:jstar])
        ql = sm / nm
        qfit = pava(y = ql, w = nm)
        ql_hat = qfit[1]
        likhood = (ql_hat ^ sm) * (1 - ql_hat) ^ (nm - sm)
        aic[l] = 2 * length(qfit) - 2 * log(likhood)
      } else{
        jstar = length(n[n != 0])
        sk = (yeff[which(n != 0)])[1:(l - 1)]
        nk = (n[n != 0])[1:(l - 1)]
        qk = sk / nk
        sm = sum(yeff[l:jstar])
        nm = sum(n[l:jstar])
        ql = sm / nm
        qfit = pava(y = c(qk, ql), w = c(nk, nm))
        qk_hat = qfit[1:(l - 1)]
        ql_hat = qfit[l]
        likhood = prod((qk_hat ^ sk) * (1 - qk_hat) ^ (nk - sk)) * (ql_hat ^
                                                                      sm) * (1 - ql_hat) ^ (nm - sm)
        aic[l] = 2 * length(qfit) - 2 * log(likhood)
      }
    }
    return(aic)
  }
  
  dselect = rep(0, ntrial)
  ndose = nrow(pmat_T)
  N = matrix(rep(0, ndose * ntrial), ncol = ndose)
  YTOX = matrix(rep(0, ndose * ntrial), ncol = ndose)
  YEFF = matrix(rep(0, ndose * ntrial), ncol = ndose)
  
  
  for (trial in 1:ntrial) {
    ytox = rep(0, ndose)
    yeff = rep(0, ndose)
    n = rep(0, ndose)
    d = 1 ##dose start from level 1
    
    ### dose-finding procedure
    for (i in 1:ncohort) {
      ytox0 = rowSums(rmultinom(cohortsize,1,prob = pmat_T[d,]))
      ytox[d] = ytox[d] + sum(ytox0*wt_Tox)/max(wt_Tox)
      yeff0 = rowSums(rmultinom(cohortsize,1,prob = pmat_E[d,]))
      yeff[d] = yeff[d] + sum(yeff0*wt_Eff)/max(wt_Eff)
      
      n[d] = n[d] + cohortsize
      AT = adm_tox(n = n, ytox = ytox)
      AR = adm_eff(n = n, yeff = yeff)
      try = length(n[n != 0])
      if ((try %in% AT) & (try < ndose)) {
        d = d + 1
      } else {
        A = intersect(AT, AR)
        if (length(A) == 0) {
          d = 0
          break
        } else {
          OBD = A[which.min(AIC(n = n, yeff = yeff)[A])]
          if (OBD > d) {
            d = d + 1
          } else if (OBD < d) {
            d = d - 1
          } else {
            d = d
          }
        }
      }
    }
    if (d == 0) {
      dselect[trial] = 0
      N[trial, ] = n
      YTOX[trial, ] = ytox
      YEFF[trial, ] = yeff
    } else{
      AT = adm_tox(n = n, ytox = ytox)
      AR = adm_eff(n = n, yeff = yeff)
      A = intersect(AT, AR)
      if (length(A) == 0){
        dselect[trial] = 0
        N[trial, ] = n
        YTOX[trial, ] = ytox
        YEFF[trial, ] = yeff
      }
      else{
        pT = rep(phiT,length(n[n!=0]))
        pT = (ytox[n!=0])/(n[n!=0])
        probT = pava(pT,n[n!=0])
        min_dif = min(abs(probT[A]-phiT))
        mtd = max(which(abs(probT-phiT)==min_dif))
        A = intersect(A,c(1:mtd))
        dselect[trial] = A[which.min(AIC(n = n, yeff = yeff)[A])]
        N[trial, ] = n
        YTOX[trial, ] = ytox
        YEFF[trial, ] = yeff
      }
    }
  }
  
  
  selpercent = rep(0, ndose + 1)
  patpercent = matrix(rep(0, ntrial * ndose), ncol = ntrial, nrow = ndose)
  efficacy = rep(0, ntrial)
  toxicity = rep(0, ntrial)
  f <- function(x) {
    x[i] / sum(x)
  }
  ## Summarize results
  
  print("True ETS")
  cat(formatC(c(ttox %*% wt_tox)/max(wt_tox), digits = 2, format = "f"),
      sep = " ", "\n")
  
  print("True EES")
  cat(formatC(c(teff %*% wt_eff)/max(wt_eff), digits = 2, format = "f"),
      sep = " ", "\n")
  
  for (i in 0:ndose) {
    selpercent[(i + 1)] = sum(dselect == i) / ntrial * 100
  }
  print("selection probablity")
  cat(formatC(selpercent, digits = 1, format = "f"), sep = " ", "\n")
  for (i in 1:ndose) {
    patpercent[i, ] = apply(N, 1, f)
  }
  print("average percent of patients")
  cat(formatC(
    apply(patpercent, 1, mean) * 100,
    digits = 1,
    format = "f"
  ),
  sep = " ", "\n")
  print("average number of patients")
  cat(formatC(c(apply(N, 2, mean), sum(apply(
    N, 2, mean
  ))), digits = 1, format = "f"),
  sep = " ", "\n")
  print("average number of patients response to efficacy")
  cat(formatC(c(apply(YEFF, 2, mean), sum(
    apply(YEFF, 2, mean)
  )), digits = 1, format = "f"),
  sep = " ", "\n")
  for (i in 1:ntrial) {
    efficacy[i] = sum(YEFF[i, ]) / sum(N[i, ])
  }
  
  print("average number of patients response to toxicity")
  cat(formatC(c(apply(YTOX, 2, mean), sum(
    apply(YTOX, 2, mean)
  )), digits = 1, format = "f"),
  sep = " ", "\n")
  for (i in 1:ntrial) {
    toxicity[i] = sum(YTOX[i, ]) / sum(N[i, ])
  }
  
}

######### continous endpoints ########
UNITED.oc = function(ttox,
                    teff,
                    wt_tox,
                    sev_wt,
                    cohortsize = 3,
                    ncohort = 20,
                    ntrial = 5000,
                    phiT = 3.04,
                    muT = 0.95,
                    phiE = 0.3,
                    muE = 0.95) {
  
  library("Iso")
  
  ###find admissible set of toxicity
  adm_tox <- function(n, ytox) {
    nn = n[n != 0]
    
    mu_hat = sapply(1:length(nn), function(i) mean(ytox[[i]]))
    
    s_square = sapply(1:length(nn), function(i) 
      sum((ytox[[i]] - mu_hat[i])^2) / (nn[i]*(nn[i]-1)))
    
    Tox_prob = 1 - pt((phiT-mu_hat)/sqrt(s_square), nn-1)
    AT_naive = which(Tox_prob < muT)
    if (length(AT_naive)==0){
      AT=AT_naive
    }
    else{
      full_seq = seq(min(AT_naive),max(AT_naive),1)
      if (length(setdiff(full_seq, AT_naive)) == 0) {
        AT = AT_naive
      }
      else{
        AT = AT_naive[AT_naive < min(setdiff(full_seq, AT_naive))]
      }
      AT = 1:max(AT)
    }
    return(AT)
  }
  
  ###find admissible set of efficacy
  adm_eff <- function(n, yeff) {
    nn = n[n != 0]
    
    mu_hat = sapply(1:length(nn), function(i) mean(yeff[[i]]))
    
    s_square = sapply(1:length(nn), function(i) 
      sum((yeff[[i]] - mu_hat[i])^2) / (nn[i]*(nn[i]-1)))
    
    Eff_prob = pt((phiE-mu_hat)/sqrt(s_square), nn-1)
    AR_naive = which(Eff_prob < muE)
    if (length(AR_naive)==0){
      AR=AR_naive
    }
    else{
      full_seq = seq(min(AR_naive),max(AR_naive),1)
      if (length(setdiff(full_seq, AR_naive)) == 0) {
        AR = AR_naive
      } else {
        AR = AR_naive[AR_naive > max(setdiff(full_seq, AR_naive))]
      }
      AR = min(AR):length(nn)
    }
    return(AR)
  }
  
  
  AIC <- function( n,y, tol = 1e-3, max_iter = 1000) {
    results <- list()
    ni<-n[n!=0]
    y_bar = mean(unlist(y))
    y_bar = rep(y_bar, length(ni))
    variance<-sapply(1:length(ni), function(i) sum((y[[i]] - y_bar[i])^2) / ni[i])
    wt <- ni / variance
    iso_fit <- pava(y_bar, wt)
    residuals <- unlist(lapply(1:length(ni), function(i) y[[i]] - iso_fit[i]))
    log_likelihood <-  - 0.5 * sum(ni * log(variance)) - 
      0.5 * sum(sapply(1:length(ni), function(i) sum((y[[i]] - iso_fit[i])^2) / variance[i]))
    k <- 1 + length(variance)
    AIC <- 2 * k - 2 * log_likelihood
    results[[1]] <- list(AIC = AIC, 
                         iso_fit = iso_fit, variance = variance, 
                         iterations = 1, inflection_point = 1)
    
    
    for (inflection_point in 1:(length(ni) - 1)) {
      
      y_bar <- sapply(1:inflection_point, function(i) mean(y[[i]]))
      if (inflection_point < length(ni)) {
        y_bar <- c(y_bar, 
                   rep(mean(unlist(y[(inflection_point + 1):length(ni)])),length(ni) - inflection_point))
        
      }
      
      variance<-sapply(1:length(ni), function(i) sum((y[[i]] - y_bar[i])^2) / ni[i])
      
      wt <- ni / variance
      
      for (iter in 1:max_iter) {
        iso_fit <- pava(y_bar, wt)
        
        new_variance <- sapply(1:length(ni), function(i) sum((y[[i]] - iso_fit[i])^2) / ni[i])
        
        if (max(abs(iso_fit - y_bar)) <= tol) {
          break
        }
        
        variance <- new_variance
        wt <- ni / variance
        y_bar <- iso_fit
      }
      
      residuals <- unlist(lapply(1:length(ni), function(i) y[[i]] - iso_fit[i]))
      log_likelihood <-  - 0.5 * sum(ni * log(variance)) - 
        0.5 * sum(sapply(1:length(ni), function(i) sum((y[[i]] - iso_fit[i])^2) / variance[i]))
      k <- (inflection_point+1) + length(variance)
      AIC <- 2 * k - 2 * log_likelihood
      
      results[[inflection_point+1]] <- list(AIC = AIC, iso_fit = iso_fit, variance = variance, iterations = iter, inflection_point = inflection_point+1)
    }
    
    aic_result <- sapply(results, function(res) res$AIC)
    
    aic_result
  }
  
  get_prob <- function(matrix, col_num) {
    
    # Extract the specified column
    column_data <- matrix[, col_num]
    
    # Define the lengths for each component
    #lengths <- c(2, 2, 2, 3, 2, 2)
    lengths <- c(4, 2, 3, 2, 2)
    # Initialize an empty list
    result_list <- vector("list", length(lengths))
    
    # Populate the list with the specified lengths
    start_index <- 1
    for (i in 1:length(lengths)) {
      end_index <- start_index + lengths[i] - 1
      result_list[[i]] <- column_data[start_index:end_index]
      start_index <- end_index + 1
    }
    
    return(result_list)
  }
  
  YT_multinomial <- function(probability,cohortsize) {
    
    # Splitting the probability vector into separate vectors for each category
    prob_list <- lapply(probability, function(p) {
      c(1 - sum(p),p)
    })
    # Generating the multinomial random variables
    result <- lapply(prob_list, function(p) {
      rmultinom(cohortsize, size = 1, prob = p)
    })
    
    return(result)
  }
  
  calculate_TTB <- function(generated_data, sev_weig) {
    weighted_sums <- mapply(function(data, weights) {
      apply(data, 2, function(col) sum(col * weights))
    }, generated_data, sev_weig, SIMPLIFY = FALSE)
    
    # Sum the values for each corresponding item in the list
    summed_columns <- Reduce("+", weighted_sums)
    
    return(summed_columns)
  }
  
  
  dselect = rep(0, ntrial)
  ndose = length(mu_eff)
  N = matrix(rep(0, ndose * ntrial), ncol = ndose)
  
  for (trial in 1:ntrial) {
    ytox = vector("list", length = ndose)
    yeff = vector("list", length = ndose)
    n = rep(0, ndose)
    d = 1 ##dose start from level 1
    
    ### dose-finding procedure
    for (i in 1:ncohort) {
      yeff0=rnorm(cohortsize,mean=mu_eff[d],sd=0.2*mu_eff[d])
      yeff[[d]]=c(yeff[[d]],yeff0)
      ytox0 = calculate_TTB(YT_multinomial(get_prob(ttox,d),cohortsize),
                            sev_wt)
      ytox[[d]]=c(ytox[[d]],ytox0)
      
      n[d] = n[d] + cohortsize
      
      AT = adm_tox(n = n, ytox = ytox)
      AR = adm_eff(n = n, yeff = yeff)
      try = length(n[n != 0])
      
      if ((try %in% AT) & (try < ndose)) {
        d = d + 1
      } else {
        A = intersect(AT, AR)
        if (length(A) == 0) {
          d = 0
          break
        } else {
          OBD = A[which.min(AIC(n, yeff)[A])]
          if (OBD > d) {
            d = d + 1
          } else if (OBD < d) {
            d = d - 1
          } else {
            d = d
          }
        }
      }
      
    }
    
    if (d == 0) {
      dselect[trial] = 0
      N[trial, ] = n
    } else{
      AT = adm_tox(n = n, ytox = ytox)
      AR = adm_eff(n = n, yeff = yeff)
      A = intersect(AT, AR)
      if (length(A) == 0){
        dselect[trial] = 0
        N[trial, ] = n
      }
      else{
        pT = rep(phiT,length(n[n!=0]))
        pT = sapply(1:length(n[n!=0]), function(i) mean(ytox[[i]]))
        probT = pava(pT,n[n!=0])
        min_dif = min(abs(probT[A]-phiT))
        mtd = max(which(abs(probT-phiT)==min_dif))
        A = intersect(A,c(1:mtd))
        dselect[trial] = A[which.min(AIC(n, yeff)[A])]
        N[trial, ] = n
      }
    }
  }
  
  selpercent = rep(0, ndose + 1)
  patpercent = matrix(rep(0, ntrial * ndose), ncol = ntrial, nrow = ndose)
  f <- function(x) {
    x[i] / sum(x)
  }
  ## Summarize results
  
  print("True TTB")
  cat(formatC(c(wt_tox%*%ttox), digits = 2, format = "f"),
      sep = " ", "\n")
  
  print("True efficacy")
  cat(formatC(c(teff), digits = 2, format = "f"),
      sep = " ", "\n")
  
  for (i in 0:ndose) {
    selpercent[(i + 1)] = sum(dselect == i) / ntrial * 100
  }
  print("selection probablity")
  cat(formatC(selpercent, digits = 1, format = "f"), sep = " ", "\n")
  for (i in 1:ndose) {
    patpercent[i, ] = apply(N, 1, f)
  }
  print("average percent of patients")
  cat(formatC(
    apply(patpercent, 1, mean) * 100,
    digits = 1,
    format = "f"
  ),
  sep = " ", "\n")
  print("average number of patients")
  cat(formatC(c(apply(N, 2, mean), sum(apply(
    N, 2, mean
  ))), digits = 1, format = "f"),
  sep = " ", "\n")
  
}
