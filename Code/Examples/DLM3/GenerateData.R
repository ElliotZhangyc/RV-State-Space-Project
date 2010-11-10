# We are considering a dynamics linear model governed by
#  y_t = x_t + \nu_t, \nu_t \sim N(m.nu ,V)
#  z_t = \mu + \phi (z_{t-1}-\mu) + \omega_t, \omega_t \sim N(0,W).

# The length of our observed time series, (y_t).
T = 1000;

# The true values of our parameters.
true = data.frame(
  mu = -3.80,
  m.nu = -0.62,
  v.nu = 1.16,
  phi = 0.8,
  W = 0.5,
);
  
# Our data structures, y and x.  For plotting purposes we keep track
# of x.  y is indexed from 1 to n.  x is indexed from 0 to n.
y.true = rep(0,T);
z.true = rep(0,T+1);

# Set z_0 to something.  Since z_0 mean reverts to zero we chose
# something that is not too far away from 0.5.  It turns out
# 0.1 puts up pretty far away so it takes a while to get back to 0.5.
# We actually know here that z_0 should come from N(mu, W/(1-phi^2)),
# but we refrain from drawing from that distribution.
z.true[1] = 0.1;

# Now generate our data.
for (i in 1:T){
  z.true[i+1] = true$mu + true$phi*(z.true[i]-true$mu) +
                rnorm(1, 0, sqrt(true$W));
  y.true[i] = z.true[i+1] + rnorm(1, true$m.nu, sqrt(true$v.nu));
}

# For purposes of calibrating.  Our CalibrateData.R script assumes
# there is a variable y.data with the observed y-values.
z.data = z.true;
y.data = y.true;
