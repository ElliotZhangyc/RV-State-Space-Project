# This script calibrates our GARCH(1,1) model, when our innovations
# are normal.

# We want to use a multivariate normal and multivariate T distribution
# so we need to load the mvtnorm library.
# library(mvtnorm);

## ASSUMPTIONS ##

# We assume that r.data holds our return data and that cv.data holds
# our conditional variance data.  We only use the first term of the
# conditional variance data, so we are cheating only a little.  We
# assume that we have normal innovations.

# We also may assume that we have already run MaxLL.R, but this may
# not be the case.

## END ASSUMPTIONS ##

# The length of our time series.
T = length(r.data);

# The number of steps we use in our Gibbs sampler.
num.samples = 1000;

## DATA STRUCTURES ##

# Set up our data structures.
tau = rep(1.0, num.samples);
theta = matrix(0, num.samples, 3);

## PRIOR PARAMETERS ##

# We assume that the prior for tau is Ga(tau.a/2, tau.b/2).
tau.a = 3.0;
tau.b = 2.0;

# We need to let everyone know that tau.current is global.
tau.current = 1.0;

# For our accept/reject.
tune = 1.0;

# Seed the starting values.
tau[1] = 1;
theta[1,] = c(0.2, 0.5, 0.05);

# SigSq
SigSq = rep(0, T);
SigSq[1] = cv.data[1];
for(i in 2:T){
  SigSq[i] = theta[1] * r.data[i-1]^2 + theta[2] * SigSq[i-1] + theta[3];
}

## HELPER FUNCTIONS ##

tau.current = 1;

## CONDITIONAL POSTERIOR DISTRIBUTIONS FOR GIBBS SAMPLING ##

tau.cond.post <- function(r.data, h.data, tau.a, tau.b){
  n = length(r.data);
  tau.a.post = a + n;
  tau.b.post = tau.b + sum( r.data^2 / h.data);
  tau.draw = rgamma(1, tau.a.post/2, rate=tau.b.post/2);
}

theta.cond.post <- function(r.data, tau, theta.previous){
  pass.test = FALSE;
  while(!pass.test){
    print(pass.test)
    cand.draw = rmvnorm(1, m.cand, C.cand);
    pass.test = prod(0 <= cand.draw) * prod(cand.draw <= 1);
  }
  # print(cand.draw);
  ratio = likelihoodkernel(cand.draw)/likelihoodkernel(theta.previous)*
    dmvnorm(theta.previous, m.cand, C.cand)/dmvnorm(cand.draw, m.cand, C.cand);
  a = min(1, ratio);
  print(a);
  u = runif(1);
  draw = cand.draw * (u<a) + theta.previous * (u>a);
}

## HELPER MISCELLANY ##

# If we assume that tau is known, then we only need to infer \theta.
# In that case we only need to find our mode and Hessian one time.
tau.current = tau[1];
# Find the optimium.
source('MaxLL.R');
# Candidate mode.
m.cand = output$par;
gradient = matrix(D.negloglikelihood(output$par), 3, 1);
hess = 0.5 * exp(-0.5 * output$value) *
  (0.5 * gradient %*% t(gradient) - output$hessian);
# Candidate covariance.
C.cand = -1 * tune * hess;

# Find the optimium.
# source('MaxLL.R');
# Candidate mode.
#m.cand = output$par;
#gradient = matrix(D.negloglikelihood(output$par), 3, 1);
#hess = exp(output$value) *
#  (gradient %*% t(gradient) + output$hessian);

## GIBBS SAMPLING ##
for(i in 2:num.samples){
  theta[i,] = theta.cond.post(r.data, tau[i-1], theta[i-1,]);
  # print(i);
  # SigSq[1] = cv.data[1];
  # for(i in 2:T){
  #   SigSq[i] = theta[1] * r.data[i-1]^2 + theta[2] * SigSq[i-1] + theta[3];
  # }
  # tau[i] = tau.cond.post(r.data, SigSq, tau.a, tau.b);
}
