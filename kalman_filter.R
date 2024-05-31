#in this script we will write the equations for the multidimensions Kalman Filter

#Step 1: Calculating the predicted state (vector same size as beta)
pred_state <- function(A, prev_state, Sigma_prev){
  return(A %*% prev_state)
}

#Step 2: Predicted process covariance matrix
pred_cov_mat <- function(A, Sigma_prev, delta){
  P <- A %*% Sigma_prev %*% t(A)
  return(P + ((1 - delta)/delta) * P)
}

#Step 3: Kalman Gain
kalman_gain <- function(H, Sigma_pred, Omega){
  return((Sigma_pred %*% H) %*% solve((H %*% Sigma_pred %*% t(H) + Omega)))
}

#Step 4: Computing the current state
current_state <- function(H, A, pred_state, kalman_gain, measurement){
  return(pred_state + kalman_gain %*% (measurement - H %*% A %*% pred_state))
}

#Step 5: Update process covariance matrix
new_Sigma <- function(kalman_gain, H, Sigma_prev){
  return((diag(1, nrow(Sigma_prev)) - kalman_gain %*% H) %*% Sigma_prev)
}

################################################################################
################################################################################
################################################################################

kalman_filter <- function(time_pts,
                          k,
                          cov,
                          beta,
                          omega,
                          delta,
                          A = NULL,
                          H = NULL,
                          Sigma_init = NULL,
                          beta_0 = NULL){
  
  ################################################################################
  ################################################################################
  ################################################################################
  #This function computes the kalman filter for a set of estimates, I am going to
  #be using it for a REM, but it should work for the estimates obtained from any
  #model.
    #time_pts: the total number of time points 
    #k: the lag in the time series (For instance, if the data is separater weekly, then k = 7. Because Monday should depend on Monday, Tuesday on Tuesday and so on..)
    #cov: number of covariates in the model
    #beta: MLE estimates from the model
    #omega: MLE's standard error from the model
    #delta: is the discount factor used to increase the uncertainty from t-1 to t
    #A: the A matrix
    #H: the H matrix
    #Sigma_init: the initial covariance matrix of the process (just a guess)
    #beta_0: the initial state (just a guess)
  ################################################################################
  ################################################################################
  ################################################################################
  if(is.null(A)){
    A <- diag(1, cov)
  }
  if(is.null(H)){
    H <- diag(1, cov)
  }
  if(is.null(Sigma_init)){
    Sigma_init <- diag(1000, cov) 
  }
  if(is.null(beta_0)){
    beta_0 <- rep(0, cov)
  }
  
  #Storing the function output
  current_state_mat <- matrix(NA, nrow = cov, ncol = total)
  kalman_gain_mat <- array(NA, dim = c(cov, cov, total))
  process_cov_mat <- array(NA, dim = c(cov, cov, total))
  
  for(i in 1:total){
    if(i %in% 1:k){
      predicted_state <- pred_state(A, beta_0, Sigma_init)
      predicted_cov_mat <- pred_cov_mat(A, Sigma_init, delta)
      kalman_gain_mat[,,i] <- kalman_gain(H, predicted_cov_mat, omega[,,i])
      current_state_mat[,i] <- current_state(H, A, predicted_state, kalman_gain_mat[,,i], beta[,i])
      process_cov_mat[,,i] <- new_Sigma(kalman_gain_mat[,,i], H, predicted_cov_mat)
    } else {
      predicted_state <- pred_state(A, current_state_mat[,i-k], process_cov_mat[,,i-k])
      predicted_cov_mat <- pred_cov_mat(A, process_cov_mat[,,i-k], delta)
      kalman_gain_mat[,,i] <- kalman_gain(H, predicted_cov_mat, omega[,,i])
      current_state_mat[,i] <- current_state(H, A, predicted_state, kalman_gain_mat[,,i], beta[,i])
      process_cov_mat[,,i] <- new_Sigma(kalman_gain_mat[,,i], H, predicted_cov_mat)
    }
  }
  
  return(list(beta = current_state_mat,
              kalman_gain = kalman_gain_mat,
              sigma = process_cov_mat)) 
  
}
