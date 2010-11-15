# We are considering a dynamic linear model governed by
#  y_t = z_t + \nu_t, \nu_t \sim Normal Mixture !!!
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

# Remember nu$m and $nu$v do not matter any more since we are using a
# mixture of normals.

## DATA STRUCTURES ##

# Now set up our data structures.  We set the values to the ``true
# values'' so that we may check out our Gibbs steps.
phi.gibbs = rep(0, mcmc$samples);
W.gibbs = rep(0, mcmc$samples);
mu.gibbs = rep(0, mcmc$samples);
z.gibbs = matrix(0, T+1, mcmc$samples);
q.gibbs = matrix(0, T, mcmc$samples);

## SEED THE ARRAYS ##

# For mu
mu.gibbs[1] = seed$mu;
# For phi
phi.gibbs[1] = seed$phi;
# For W
W.gibbs[1] = seed$W;

# F or z.
z.gibbs[1,1] = seed$z.0;
for(i in 2:(T+1)){
  z.gibbs[i,1] = mu.gibbs[1] + phi.gibbs[1]*(z.gibbs[i-1,1]-mu.gibbs[1]) +
    rnorm(1, 0, sqrt(W.gibbs[1]));
}

# For q.  We really do not need to do this given how we've set up our sampler.
q.gibbs[,1] = sample(1:nrow(nmix), T, replace=TRUE, prob=nmix$q);

## TEMPORARY CALCULATIONS ##

#mu = rep(true$mu, mcmc$samples);
#phi.gibbs = rep(true$phi, mcmc$samples);
#W = rep(true$phi, mcmc$samples);

# For q
#for(i in 1:mcmc$samples){
#  q.gibbs[,i] = q.data;
#}

#for(i in 1:mcmc$samples){
#  z.gibbs[,i] = z.data;
#}

## CONDITIONAL DENSITIES ##

# Let's set up the function, which will create the draws from our
# conditional posteriors.

# To draw from x | everything else.
z.cond.post <- function(y.data, mu, q, phi, W, m.0, C.0){
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
    Q = R[i] + nmix$v[q[i]];
    A = R[i]/Q;
    # NOTE: We include m.mu here to accomodate nu with nonzero mean.
    m[i+1] = mu + phi*(m[i]-mu) +
             A*( y.data[i] - (mu + phi*(m[i]-mu) + nmix$b[q[i]]) );
    C[i+1] = A*nmix$v[q[i]];
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
  # Draw from an approximate posterior distribution.
  phi.draw = rtnorm(1, m.inter, sqrt(C.inter), 0, 1);
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

# To generate the conditional posterior distribution of q.
q.cond.post <- function(y.data, z.data){
  n = length(y.data);
  num.from.mix = nrow(nmix);
  q.draw = 1:n;
  q.t.post = 1:num.from.mix;
  # To speed things up:
  log.n.over.sqrt.v = log(nmix$q) - 0.5 * log(nmix$v);
  # Crunch through things:
  for(t in 1:n){
    # To speed things up:
    y.from.z = y.data[t] - z.data[t+1];
    #for(i in 1:num.from.mix){
    #  q.t.post[i] = exp(-0.5*(y.data[t] - z.data[t+1] - nmix$b[i])^2/nmix$v[i]) *
    #                nmix$q[i]/sqrt(nmix$v[i]);
    #}
    # HUGE SPEED UP BY VECTORIZING.  SMALLER IMPROVEMENT BY REDUCING OPERATIONS.
    q.t.post = exp(-0.5 * (y.from.z - nmix$b)^2 / nmix$v + log.n.over.sqrt.v);
    #print(c(t,q.t.post));
    nrm = sum(q.t.post);
    q.t.post = q.t.post / nrm;
    q.draw[t] = sample(1:num.from.mix, 1, replace=TRUE, prob=q.t.post);
  }
  q.draw;
}

## GIBBS SAMPLING ##

# To keep track of the time.
the.time = rep(0, floor(mcmc$samples/100)+1);
delta.time = rep(0, length(the.time)-1);
the.time[1] = proc.time()[[1]];

# Now we can go ahead and do our Gibbs sampling.
# Again refer to the notes for a description.
for(i in 2:mcmc$samples){
  # A change of variables is useful.
  x.data = z.gibbs[,i-1] - mu.gibbs[i-1];
  # Sample q | e.e.
  q.gibbs[,i] = q.cond.post(y.data, z.gibbs[,i-1]);
  # Sample phi | e.e.
  phi.gibbs[i] = phi.cond.post(x.data, W.gibbs[i-1],
             prior$m.phi, prior$v.phi, phi.gibbs[i-1]);
  # Sample W | e.e.
  W.gibbs[i] = W.cond.post(x.data, phi.gibbs[i], prior$a.W, prior$b.W);
  # Sample mu | e.e.
  mu.gibbs[i] = mu.cond.post(z.gibbs[,i-1], phi.gibbs[i],
            W.gibbs[i], prior$m.mu, prior$v.mu);
  # Sample z | e.e.
  z.gibbs[,i] = z.cond.post(y.data, mu.gibbs[i], q.gibbs[,i],
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
#range = 3000:mcmc$samples;
#QoI = mu[range];
#m.QoI = mean(QoI);
#sd.QoI = sd(QoI);
#hist(QoI, breaks=40, prob=TRUE);
#grid = seq(min(QoI), max(QoI), 0.01);
#dens = dnorm(grid, m.QoI, sd.QoI);
#lines(grid, dens);
#print(c(m.QoI, sd.QoI));
