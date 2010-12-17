# We are considering a dynamics linear model governed by
#  y_t = z_t + \nu_t, \nu_t \sim Normal Mixture !!!
#  z_t = \mu + \phi (z_{t-1}-\mu) + \omega_t, \omega_t \sim N(0,W).

# Our data structures for y, z, and q.  For plotting purposes we keep track
# of z.  y is indexed from 1 to T.  z is indexed from ``0 to T''.
y.true = rep(0,T);
z.true = rep(0,T+1);
q.true = rep(0,T);

# Set z_0 to something.  Since z_0 mean reverts to zero we chose
# something that is not too far away from the mean mu.
z.true[1] = true$z.0;

# Now generate our data.
for (i in 1:T){
  z.true[i+1] = true$mu + true$phi*(z.true[i]-true$mu) +
                rnorm(1, 0, sqrt(true$W));
  q.true[i] = sample(1:6, 1, replace=TRUE, prob=nmix$q);
  y.true[i] = z.true[i+1] + rnorm(1, nmix$b[q.true[i]], sqrt(nmix$v[q.true[i]]));
}

# For purposes of calibrating.  Our CalibrateData.R script assumes
# there is a variable y.data with the observed y-values.
z.data = z.true;
y.data = y.true;
q.data = q.true;
