# We are considering a dynamics linear model governed by
#  y_t = \mu + x_t + \nu_t, \nu_t \sim N(m.nu ,V)
#  x_t = \phi x_{t-1} + \omega_t, \omega_t \sim N(0,W).

# The length of our observed time series, (y_t).
T = 100;

# The values of our parameters.
mu.true = 0.5;
m.nu.true = 0.3;
V.true = 0.2;
phi.true = 0.3;
W.true = 0.5;

# Our data structures, y and x.  For plotting purposes we keep track
# of x.  y is indexed from 1 to n.  x is indexed from 0 to n.
y = rep(0,T);
x = rep(0,T+1);

# Set x_0 to something.  Since x_0 mean reverts to zero we chose
# something that is not too far away from 0.
x[1] = 0.1;

# Now generate our data.
for (i in 1:T){
  x[i+1] = phi.true*x[i] + rnorm(1, 0, sqrt(W.true));
  y[i] = mu.true + x[i+1] + rnorm(1, m.nu.true, sqrt(V.true));
}

# For purposes of calibrating.  Our CalibrateData.R script assumes
# there is a variable y.data with the observed y-values.
x.data = x;
y.data = y;
