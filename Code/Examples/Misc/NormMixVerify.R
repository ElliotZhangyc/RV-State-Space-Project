# These quantities describe a normal mixture that approximates a 0.5 *
# log(chi-square(1)) distribution.  In particular the distribution is
# governed by p(x | i) \sim N( b[i], sqrt(w[i]) ), p(i) = q(i).
q = c(0.0073, 0.1056, 0.2575, 0.34, 0.2456, 0.0440);
b = c(-5.7002, -2.6216, -1.1793, -0.3255, 0.2624, 0.7537);
w = c(1.4490, 0.6534, 0.3157, 0.16, 0.0851, 0.0418);

n = length(q);

N = 10000;
data = 1:N;

#draws = 1000;
#data = sample(1:n, draws, replace=TRUE, prob=q);

# When plotting a histogram things do not look like if you try to turn
# the histogram into a probability density function using prob=TRUE.

for(i in 1:N){
  idx = sample(1:n, 1, replace=TRUE, prob=q);
  data[i] = rnorm(1, b[idx], sqrt(w[idx]));
}

# A chi-square(k) distribution is Ga(k/2, rate=1/2).
chisquare = rgamma(N, 1/2, rate=1/2);
logchisq = 0.5 * log(chisquare);

# To see that our normal mixture approximates a chisquare...
actual = hist(logchisq, breaks=40, plot=FALSE);
approx = hist(data, breaks=40, prob=TRUE);
lines(actual$mids, actual$density);
