# We are considering a dynamic linear model governed by
#  y_t = z_t + \nu_t, \nu_t \sim N(m.nu ,V)
#  z_t = \mu + \phi (z_{t-1} - \mu) + \omega_t, \omega_t \sim N(0,W).

# We want to generate the posterior distribution through
# simulation for our x values and our parameter values.  We will use a
# Gibbs sampler.  Thus we break our draw from
#   p(phi, W, mu, x[0:N] | y[1:N])
# into steps
#   p(phi | everything else)
#   p(W | everything else)
#   p(mu | everything else)
#   p(x[0:n] | everything else).

# We chose conjugate priors in most cases, However, we assume that z_t
# has reached equilibrium and hence z_0 has disribution N(mu,
# W/(1-\phi^2)).  This is identical to assuming that x_t = z_t -
# \mu has distribution N(0, W/(1-\phi^2))

## ASSUMPTIONS ##

# We assume our data is stored in y.data

# The length of our time series.
T = length(y.data);

## DATA STRUCTURES ##

# Now set up our data structures.  We set the values to the ``true
# values'' so that we may check out our Gibbs steps.
mu.gibbs = rep(0, mcmc$samples);
phi.gibbs = rep(0, mcmc$samples);
W.gibbs = rep(0, mcmc$samples);
# R is column centric!  But we didn't see a speed up when we changed this.
z.gibbs = matrix(0, T+1, mcmc$samples);

## SEED VALUES ##

# Let's seed our sampler with some values.

# mu[1] = mu.true;
mu.gibbs[1] = seed$mu;
# phi[1] = phi.true;
phi.gibbs[1] = seed$phi;
# W[1] = W.true;
W.gibbs[1] = seed$W;
# z[1,1] = 0.2;
z.gibbs[1,1] = seed$z.0;

# We need to generate the rest of the z values.
for(i in 2:(T+1)){
  z.gibbs[i,1] = mu.gibbs[1] +
    phi.gibbs[1]*(z.gibbs[i-1,1]-mu.gibbs[1]) +
    rnorm(1, 0, sqrt(W.gibbs[1]));
}

## TEMPORARY CALCULATIONS ##

# Useful for check Gibbs sampler.
#mu.gibbs = rep(true$mu, mcmc$samples);
#phi.gibbs = rep(true$phi, mcmc$samples);
#W.gibbs = rep(true$W, mcmc$samples);
#for(i in 1:mcmc$samples){
#  z.gibbs[,i] = z.data;
#}

## CONDITIONAL DENSITIES ##

# Let's set up the function, which will create the draws from our
# conditional posteriors.

# To draw from x | everything else.
z.cond.post <- function(y.data, mu, m.nu, V, phi, W, m.0, C.0){
  # We need to keep track of some information as we feed forward.  The
  # procedure we use here can be found in my notes.  We need an extra
  # spot for m.0 and C.0 in our m an C arrays.  Since R forces us to
  # start indexing at 1, we have that m[1] = m.0 and C[1] = C.0.  We
  # must adjust our indexing accordingly.
  n = length(y.data);
  R = rep(0, n);
  m = rep(0, n+1); m[1] = m.0;
  C = rep(0, n+1); C[1] = C.0;
  # Now let's feed forward.
  for(i in 1:n){
    R[i] = phi^2 * C[i] + W;
    Q = R[i] + V;
    A = R[i]/Q;
    # NOTE: We include m.mu here to accomodate nu with nonzero mean.
    m[i+1] = mu + phi*(m[i]-mu) + A*( y.data[i] - (mu + phi*(m[i]-mu) + m.nu) );
    C[i+1] = A*V;
  }
  # Now let's go backwards.
  # We need a data structure to hold our draw from z.
  z.draw = rep(0, n+1);
  # First draw z_n | y.data.
  z.draw[n+1] = rnorm(1, m[n+1], sqrt(C[n+1]) );
  # Now we draw backwards, conditionally on our previous values of z.
  for(i in n:1){
    A = phi*C[i]/R[i];
    m.z.draw = m[i] + A * (z.draw[i+1] - (mu + phi*(m[i]-mu)) );
    C.z.draw  = C[i]*W/R[i];
    # Now draw z_{i-1}, which is z[i] since we have an z_0 value.
    z.draw[i] = rnorm(1, m.z.draw, sqrt(C.z.draw) );
  }
  # Return z.draw.
  z.draw;
  # Return parameters.
  # matrix( c(m,C), T+1, 2 );
}

# To draw from phi | everything else.
# Remeber that z.data[1:(T+1)] = {z_0, \ldots, z_n}.
# x.data = z.data - mu.
phi.cond.post <- function(x.data, W, m.phi, C.phi, phi.prev){
  n = length(x.data) - 1;
  g.0 = x.data[1:n] %*% x.data[1:n];
  g.1 = x.data[2:(n+1)] %*% x.data[1:n];
  # Some useful quantities.  Refer to notes.
  # These values correspond to our ``intermediate'' density.
  Q = C.phi + W / g.0;
  m.inter = m.phi * (W/g.0) / Q + (g.1/g.0) * C.phi / Q;
  C.inter = C.phi * (W/g.0) / Q;
  # Now we do our accept reject algorithm
  # phi.cand = rnorm(1, m.inter, sqrt(C.inter));
  #numer = dnorm(phi.cand, m.inter, sqrt(C.inter)) *
  #        dnorm(x.data[1], 0, sqrt(W/(1-phi.cand^2)));
  #denom = dnorm(phi.prev, m.inter, sqrt(C.inter)) *
  #        dnorm(x.data[1], 0, sqrt(W/(1-phi.prev^2)));
  #ratio = numer/denom;
  #phi.draw = phi.prev;
  #if (phi.cand<0 || phi.cand>1) ratio = 0;
  #if (runif(1,0,1) < ratio) phi.draw = phi.cand;
  #phi.draw = phi.cand;
  #if (phi.cand<0 || phi.cand>1) phi.draw = phi.prev;
  phi.draw = rtnorm(1, m.inter, sqrt(C.inter), 0, 1);
  #phi.draw;
}

# To draw from mu | everything else.
mu.cond.post <- function(z.data, phi, W, m.mu, C.mu){
  # Refer to notes for specifics.
  n = length(z.data)-1;
  # We need the quantity a from the notes.
  # In particular, look at the notes for ``Model 2.''
  a = mean(z.data[2:(n+1)]- phi*z.data[1:n]);
  Q = (1-phi)^2 * n * C.mu + W;
  m.mu.draw = m.mu * W / Q + a * (1-phi) * n * C.mu / Q;
  C.mu.draw = C.mu * W / Q;
  mu.draw = rnorm(1, m.mu.draw, sqrt(C.mu.draw));
}

# To draw from W | everything else.
W.cond.post <- function(x.data, phi, W.a, W.b){
  # Again, review notes for details.
  n = length(x.data)-1;
  a.post = W.a + n + 1;
  diff = x.data[2:(n+1)] - phi*x.data[1:n];
  b.post = W.b + diff %*% diff + (1-phi^2) * x.data[1];
  W.recip.draw = rgamma(1, a.post/2, rate=b.post/2);
  W.draw = 1/W.recip.draw;
}

## THE FOLLOWING IS FOR DRAWING (phi, W) SIMULTANEOUSLY ##

# The proposal density for (phi, W).
phi.W.prop.dens <- function(phi, W, n.m, n.v, x.0, a, b){
  ig.term = W^(-0.5*a - 1) * exp(-0.5*b/W);
  n.term = dnorm(phi, n.m, n.v);
  ans = n.term * ig.term;
}

# The target density for (phi, W).
phi.W.target.dens <- function(phi, W, n.m, n.v, x.0, a, b){
  ig.term = W^(-0.5*a - 1) * exp(-0.5*b/W);
  n.term = dnorm(phi, n.m, n.v) * sqrt(1-phi^2);
  ans = n.term * ig.term;
  if (phi < 0 || phi > 1) ans = 0;
  ans;
}

# importance = target / proposal.  W is built into the proposal and
# target mean and variance.  We get cancellation of the inverse gamma
# distributions.
log.phi.W.importance <- function(phi, target.m, target.v, S){
  diff = 0.5 * log(S) - 0.5 * (S-1)/S * (phi - target.m)^2 / target.v +
    0.5 * log(1-phi^2);
  if (phi < 0 || phi > 1) diff = -Inf;
  diff;
}
  
# To draw (phi, W) | everything else.
phi.W.cond.post <- function(x.data, m.phi, C.phi, W.a, W.b, phi.prev, W.prev){
  # We can factor (phi, W) | e.e. as something that looks very similar
  # to a normal-inverse Gamma distribution.  See notes.
  n = length(x.data) - 1;
  # Some preliminary data.
  g.0 = x.data[2:n] %*% x.data[2:n];
  g.1 = x.data[2:(n+1)] %*% x.data[1:n];
  phi.1 = g.1 / g.0; # A point estimate of phi.
  # First, we have the IG portion.
  a.post = W.a + n;
  b.post = W.b + (g.0 + x.data[1]^2 + x.data[n+1]^2) - g.1^2/g.0;
  W.recip.draw = rgamma(1, a.post/2, rate=b.post/2);
  W.draw = 1/W.recip.draw;
  W = W.draw; # To keep things brief.
  # Second, we have the normal portion.
  # We calculate the precision, variance, and mean.
  # post.prec = 1/C.phi + g.0/W.draw;
  llh.var = W / g.0;
  post.var = C.phi * llh.var / (C.phi + llh.var);
  post.mean = C.phi * phi.1 / (C.phi + llh.var) +
              llh.var * m.phi / (C.phi + llh.var);
  # We will draw from this distribution with a little larger variance.
  # In West's notes he suggests drawing from a distribution with
  # slighly larger variance when there is only one dimension.
  S = 2; # Maybe I should try 2.4.
  prop.mean = post.mean;
  prop.sd = sqrt(S*post.var);
  phi.draw = rnorm(1, prop.mean, prop.sd);
  # In the end, we will only accept this if it is in (0,1).
  # Set the draw.
  draw = c(phi.prev, W.prev);
  # Only do this if we can accept the draw.
  if(phi.draw > 0 && phi.draw <1){
  # Now we need to calculate our acceptance ratio.  See p. 393 of R&C.
  # First, we need to calculate some values so that we can evaluate
  # the importance function for the the previous phi and W.
  llh.var.prev = W.prev / g.0;
  post.var.prev = C.phi * llh.var / (C.phi + llh.var);
  post.mean.prev = C.phi * phi.1 / (C.phi + llh.var) +
              llh.var * m.phi / (C.phi + llh.var);
  # Now calculate the acceptance probability.
  log.ratio = log.phi.W.importance(phi.draw, post.mean, post.var, S) -
    log.phi.W.importance(phi.prev, post.mean.prev, post.var.prev, S);
  ratio = exp(log.ratio);
  accept = min(1, ratio);
  #print(c(post.mean, post.var, post.mean.prev, post.var.prev, ratio));
  # Now determine if we will accept the draw.
  u = runif(1);
  if (u < accept) draw = c(phi.draw, W.draw);
  }
  draw;
}

## GIBBS SAMPLING ##

# To keep track of the time.
the.time = rep(0, floor(mcmc$samples/100)+1);
delta.time = rep(0, length(the.time)-1);
the.time[1] = proc.time()[[1]];

# Now we can go ahead and do our Gibbs sampling.  Again refer to the
# notes for a description.  Make sure to sample mu and z last.
for(i in 2:mcmc$samples){
  # A change of variables is useful.
  x.data = z.gibbs[,i-1] - mu.gibbs[i-1];
  # Sample phi | e.e.
#  phi.gibbs[i]=phi.cond.post(x.data, W.gibbs[i-1],
#             prior$m.phi, prior$v.phi, phi.gibbs[i-1]);
  # Sample W | e.e.
#  W.gibbs[i] = W.cond.post(x.data, phi.gibbs[i], prior$a.W, prior$b.W);
  # Sample (phi, W) | e.e.
  phi.and.W = phi.W.cond.post(x.data, prior$m.phi, prior$v.phi,
    prior$a.W, prior$b.W, phi.gibbs[i-1], W.gibbs[i-1]);
  phi.gibbs[i] = phi.and.W[1];
  W.gibbs[i] = phi.and.W[2];
  # Sample mu | e.e.
  mu.gibbs[i] = mu.cond.post(z.gibbs[,i-1], phi.gibbs[i],
            W.gibbs[i], prior$m.mu, prior$v.mu);
  # Sample z | e.e.
  z.gibbs[,i] = z.cond.post(y.data, mu.gibbs[i], nu$m, nu$v,
           phi.gibbs[i], W.gibbs[i], mu.gibbs[i],
           W.gibbs[i]/(1-phi.gibbs[i]^2));
  # Record the time if this is our 100th draw and output time estimates.
  if (i %% 100 == 0){
    idx = i/100 + 1; # We must increment by 1 because our array starts at 0.
    the.time[idx] = proc.time()[[1]];
    delta.time[idx-1] = the.time[idx] - the.time[idx-1];
    ave.time = mean(delta.time[1:(idx-1)]);
    min.to.complete = (mcmc$samples - i)/6000*ave.time;
    print(paste("Iter:", i, "|",
                round(delta.time[idx-1], 3), "s last 100", "|",
                round(ave.time, 3), "s per 100", "|",
                round(min.to.complete, 3), "min to done"));
  }
}

## RECORD AVERAGE TIME ##
mcmc$ave.time = ave.time;

## CHECKING OUTPUT ##

# For normal posterior distributions.  QoI = quantity of interest.
#range = 3000:num.samples;
#QoI = mu[range];
#m.QoI = mean(QoI);
#sd.QoI = sd(QoI);
#hist(QoI, breaks=40, prob=TRUE);
#grid = seq(min(QoI), max(QoI), 0.01);
#dens = dnorm(grid, m.QoI, sd.QoI);
#lines(grid, dens);
#print(c(m.QoI, sd.QoI));
