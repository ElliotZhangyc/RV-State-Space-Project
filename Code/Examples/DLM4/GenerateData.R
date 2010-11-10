# We are considering a dynamics linear model governed by
#  y_t = x_t + \nu_t, \nu_t \sim Normal Mixture !!!
#  z_t = \mu + \phi (z_{t-1}-\mu) + \omega_t, \omega_t \sim N(0,W).

# The length of our observed time series, (y_t).
T = 10000;

# The values of our parameters.
mu.true = 0.5;
phi.true = 0.8;
W.true = 0.5;

# Our data structures, y and x.  For plotting purposes we keep track
# of x.  y is indexed from 1 to n.  x is indexed from 0 to n.
y.true = rep(0,T);
z.true = rep(0,T+1);
q.true = rep(0,T);

# Set x_0 to something.  Since x_0 mean reverts to zero we chose
# something that is not too far away from 0.
z.true[1] = 0.1;

# For our Normal Mixture
q.prior = c(0.0073, 0.1056, 0.2575, 0.34, 0.2456, 0.0440);
b.mix = c(-5.7002, -2.6216, -1.1793, -0.3255, 0.2624, 0.7537);
v.mix = c(1.4490, 0.6534, 0.3157, 0.16, 0.0851, 0.0418);

# Now generate our data.
for (i in 1:T){
  z.true[i+1] = mu.true + phi.true*(z.true[i]-mu.true) +
                rnorm(1, 0, sqrt(W.true));
  q.true[i] = sample(1:6, 1, replace=TRUE, prob=q.prior);
  y.true[i] = z.true[i+1] + rnorm(1, b.mix[q.true[i]], sqrt(v.mix[q.true[i]]));
}

# For purposes of calibrating.  Our CalibrateData.R script assumes
# there is a variable y.data with the observed y-values.
z.data = z.true;
y.data = y.true;
q.data = q.true;
