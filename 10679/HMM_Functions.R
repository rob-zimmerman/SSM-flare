# compute logged forward probabilities
logalpha_t <- function(t, init, tpm, data, fx, all = T) {
  Fx <- function(x) {
    diag(fx(x))
  }
  
  alpha <- matrix(0L, nrow = nrow(tpm), ncol = t)
  
  at_tilde <- init %*% Fx(data[, 1])
  norm_at_tilde <- sum(at_tilde)
  at_hat <- at_tilde/norm_at_tilde
  c <- log(norm_at_tilde)
  
  alpha[, 1] <- c + log(at_hat)
  
  for (s in 2:t) {
    at_tilde <- at_hat %*% tpm %*% Fx(data[, s])
    norm_at_tilde <- sum(at_tilde)
    at_hat <- at_tilde/norm_at_tilde
    c <- c + log(norm_at_tilde)
    
    alpha[, s] <- c + log(at_hat)
  }
  
  return(alpha)
  
}

# compute logged backwards probabilities
logbeta_t <- function(t, tpm, data, fx, all = T) {
  Fx <- function(x) {
    diag(fx(x))
  }
  
  TT <- ncol(data)
  
  beta <- matrix(0, nrow = nrow(tpm), ncol = TT)
  
  bt_tilde <- tpm %*% Fx(data[, TT]) %*% rep(1, times = nrow(tpm))
  norm_bt_tilde <- sum(bt_tilde)
  bt_hat <- bt_tilde/norm_bt_tilde
  c <- log(norm_bt_tilde)
  
  beta[, TT - 1] <- c + log(bt_hat)
  
  for (s in (TT - 1):1) {
    bt_tilde <- tpm %*% Fx(data[, s]) %*% bt_hat
    norm_bt_tilde <- sum(bt_tilde)
    bt_hat <- bt_tilde/norm_bt_tilde
    c <- c + log(norm_bt_tilde)
    
    beta[, s - 1] <- c + log(bt_hat)
  }
  
  return(beta)
}

# local decoding via forward-backward algorithm (adapted from Zucchini et al (2016))
local_dec <- function(init, tpm, data, Fx) {
  TT <- ncol(data)
  K <- ncol(tpm)
  la <- logalpha_t(TT, init, tpm, data, fx = Fx, all = T)
  lb <- logbeta_t(TT, tpm, data, fx = Fx, all = T)
  c <- max(la[, TT])
  llk <- c + log(sum(exp(la[, TT] - c)))
  stateprobs <- matrix(NA, ncol = TT, nrow = K)
  for (i in 1:TT) stateprobs[, i] <- exp(la[, i] + lb[, i] - llk)
  ild <- rep(NA, TT)
  for (i in 1:TT) ild[i] <- which.max(stateprobs[, i])
  return(list(probs = stateprobs, preds = ild))
}

# convert signed Infs to max-precision values (from 'copula' package)asFinite <- function(x) {  if (any(nifi <- !is.finite(x)))    x[nifi] <- sign(x[nifi]) * .Machine$double.xmax  x}