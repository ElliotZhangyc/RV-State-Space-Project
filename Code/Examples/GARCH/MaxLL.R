# Let's see how we do minimizing our log likelihood.  We also use this
# script to load the functions defined below.

## IMPORTANT ASSUMPTIONS ##

# We assume theta = (alpha, beta, omega).

# We assume we have the return data and the conditional variance data
# stored in r.data and cv.data.

## END ASSUMPTIONS ##

# Let's use tau.current for the current value of tau.
# tau.current = 1.0;

# The function -2 * log likelihood.
negloglikelihood <- function(theta){
  n = length(r.data);
  # Let h_t = \sigma_t^2.  Recall that \sigma_t^2 is recursively
  # defined by an initial value \sigma_1^2 and the vector of
  # parameters \theta = (\alpha, \beta, \omega).
  h = rep(0, n);
  # Use our known (and in this case, it is truly known) initial
  # conditional variance.
  h[1] = cv.data[1];
  # The first term of the sum.
  llhood = log(h[1]) + tau.current * r.data[1]^2/h[1];
  for(i in 2:n){
    # The next term in the sequence of conditional variances.
    h[i] = theta[1]*r.data[i-1]^2 + theta[2]*h[i-1] + theta[3];
    # Add this to our cummulative sum.
    llhood = llhood + log(h[i]) + tau.current * r.data[i]^2/h[i];
  }
  llhood;
}

# The gradient of -2 * log likelihood.
D.negloglikelihood <- function(theta){
  n = length(r.data);
  # Let h_t = \sigma_t^2.  Recall that \sigma_t^2 is recursively
  # defined by an initial value \sigma_1^2 and the vector of
  # parameters \theta = (\alpha, \beta, \omega).
  h = rep(0, n);
  # We need to keep track of the three components of our gradient.
  dh.domega = rep(0,n);
  dh.dalpha = rep(0,n);
  dh.dbeta = rep(0,n);
  # We know the first element of our sequence of conditional variances.
  h[1] = cv.data[1];
  # We know that the likelihood is zero to start since dh[1]/dwhatever = 0.
  dllhood.domega = 0;
  dllhood.dalpha = 0;
  dllhood.dbeta = 0;
  # Now recursively calculate things.
  for(i in 2:n){
    # The conditional variance.
    h[i] = theta[1]*r.data[i-1]^2 + theta[2]*h[i-1] + theta[3];
    # The derivatives.
    dh.domega[i] = theta[2] * dh.domega[i-1] + 1;
    dh.dalpha[i] = theta[2] * dh.dalpha[i-1] + r.data[i-1]^2;
    dh.dbeta[i] = h[i-1] + theta[3] * dh.dbeta[i-1];
    # This term is present in each component of the gradient.
    placeholder = (1/h[i] - tau.current * r.data[i]^2/h[i]^2);
    # Keep summing it all up.
    dllhood.domega = dllhood.domega + placeholder * dh.domega[i];
    dllhood.dalpha = dllhood.dalpha + placeholder * dh.dalpha[i];
    dllhood.dbeta = dllhood.dbeta + placeholder * dh.dbeta[i];
  }
  # Return the gradient.
  c(dllhood.dalpha, dllhood.dbeta, dllhood.domega);
}

# The likelihood.
# Remember that h_t = \sigma_t^2.
likelihoodkernel <- function(theta){
  h = rep(0, T);
  h[1] = cv.data[1];
  llhood = log(h[1]) + tau.current * r.data[1]^2/h[1];
  for(i in 2:T){
    h[i] = theta[1]*r.data[i-1]^2 + theta[2]*h[i-1] + theta[3];
    llhood = llhood + log(h[i]) + tau.current * r.data[i]^2/h[i];
  }
  exp(-0.5 * llhood);
}

# The initial values to start our search.
iv = c(0.5, 0.5, 0.5);

# Let's see what happens.  I have included several different examples.

# The NULL parameter here tells optim to calculate the gradient using
# finite differences.  If you can calulate the gradient analytically,
# you may include that function in place of NULL.

#output1 = optim(iv, negloglikelihood, NULL, method="L-BFGS-B",
#  lower=c(0.01, 0.01, 0.01), upper=c(1, 1, 0.3),
#  control = c(factr=1e3), hessian=TRUE);

# Now we use our gradient of the loglikelihood function.
output = optim(iv, negloglikelihood, D.negloglikelihood, method="L-BFGS-B",
  lower=c(0.01, 0.01, 0.01), upper=c(1, 1, 0.3),
  hessian=TRUE);

#output = optim(iv, negloglikelihood, D.negloglikelihood, method="CG",
#  hessian=TRUE);

# In included this as a second optimiztion.  The idea is that you do a
# first run optimization to get close to the solution and then
# optimization again to do some small scale refinement.  I think this
# is desirable given how optim decides when it has found an estremum.

output3 = optim(output$par, negloglikelihood, D.negloglikelihood,
  method="L-BFGS-B",
  lower=c(0.01, 0.01, 0.01), upper=c(1, 1, 0.3),
  control = c(factr=1e3), hessian=TRUE);
