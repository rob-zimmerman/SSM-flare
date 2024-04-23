# local decoding via forward-backward algorithm (adapted from Zucchini et al (2016))
local_dec <- function(init, tpm, data, Fx) {
  TT <- ncol(data)
  K <- ncol(tpm)
  
  # precompute F(data)
  Fdata <- matrix(0L, nrow=K, ncol=TT)
  for (t in 1:TT) {
    Fdata[,t] <- Fx(data[, t])
  }
  
  # compute logged forward probabilities
  alpha <- matrix(0L, nrow = nrow(tpm), ncol = TT)
  
  at_tilde <- init * Fdata[,1]
  norm_at_tilde <- sum(at_tilde)
  at_hat <- at_tilde/norm_at_tilde
  c <- log(norm_at_tilde)
  
  alpha[, 1] <- c + log(at_hat)
  
  for (s in 2:TT) {
    print(s)
    at_tilde <- at_hat %*% (t(t(tpm) * Fdata[,s]))
    norm_at_tilde <- sum(at_tilde)
    at_hat <- at_tilde/norm_at_tilde
    c <- c + log(norm_at_tilde)
    
    alpha[, s] <- c + log(at_hat)
  }
  
  
  
  # compute logged backward probabilities
  beta <- matrix(0, nrow = nrow(tpm), ncol = TT)
  
  bt_tilde <- t(t(tpm) * Fdata[,TT]) %*% rep(1, times = nrow(tpm))
  norm_bt_tilde <- sum(bt_tilde)
  bt_hat <- bt_tilde/norm_bt_tilde
  c <- log(norm_bt_tilde)
  
  beta[, TT-1] <- c + log(bt_hat)
  
  for (s in (TT-1):1) {
    print(s)
    bt_tilde <- t(t(tpm) * Fdata[, s]) %*% bt_hat
    norm_bt_tilde <- sum(bt_tilde)
    bt_hat <- bt_tilde/norm_bt_tilde
    c <- c + log(norm_bt_tilde)
    
    beta[, s-1] <- c + log(bt_hat)
  }
  
  c <- max(alpha[, TT])
  llk <- c + log(sum(exp(alpha[, TT] - c)))
  stateprobs <- matrix(0, ncol = TT, nrow = K)
  for (i in 1:TT) {
    stateprobs[, i] <- exp(alpha[, i] + beta[, i] - llk)
  }
  preds <- rep(0, TT)
  for (i in 1:TT) {
    preds[i] <- which.max(stateprobs[, i])
  }
  return(list(probs = stateprobs, preds = preds))
}

# convert signed Infs to max-precision values (from 'copula' package)
asFinite <- function(x) {
  if (any(nifi <- !is.finite(x)))
    x[nifi] <- sign(x[nifi]) * .Machine$double.xmax
  x
}