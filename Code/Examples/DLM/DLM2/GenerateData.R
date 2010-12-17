# We are considering a dynamics linear model governed by
#  y_t = x_t + \nu_t, \nu_t \sim N(m.nu ,V)
#  z_t = \mu + \phi (z_{t-1}-\mu) + \omega_t, \omega_t \sim N(0,W).

# The length of our observed time series, (y_t).
T = 10000;

# The values of our parameters.
mu.true = 0.5;
m.nu.true = 0.3;
V.true = 0.2;
phi.true = 0.8;
W.true = 0.5;

# Our data structures, y and x.  For plotting purposes we keep track
# of x.  y is indexed from 1 to n.  x is indexed from 0 to n.
y.true = rep(0,T);
z.true = rep(0,T+1);

# Set x_0 to something.  Since x_0 mean reverts to zero we chose
# something that is not too far away from 0.
z.true[1] = 0.1;

# Now generate our data.
for (i in 1:T){
  z.true[i+1] = mu.true + phi.true*(z.true[i]-mu.true) +
    rnorm(1, 0, sqrt(W.true));
  y.true[i] = z.true[i+1] + rnorm(1, m.nu.true, sqrt(V.true));
}

# For purposes of calibrating.  Our CalibrateData.R script assumes
# there is a variable y.data with the observed y-values.
z.data = z.true;
y.data = y.true;
