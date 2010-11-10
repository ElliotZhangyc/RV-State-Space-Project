#phi.v2 = phi;
#W.v2 = W;
#mu.v2 = mu;
#z.v2 = z;

#phi.v2.hist = hist(phi.v2[3000:num.samples], breaks=40, plot=FALSE);
#W.v2.hist = hist(W.v2[3000:num.samples], breaks=40, plot=FALSE);
#mu.v2.hist = hist(mu.v2[3000:num.samples], breaks=40, plot=FALSE);
#z.v2.hist = hist(z.v2[3000:num.samples,50], breaks=40, plot=FALSE);

hist(phi[3000:num.samples], breaks=40, prob=TRUE)
lines(phi.v2.hist$mids, phi.v2.hist$density)
#hist(mu[3000:num.samples], breaks=40, prob=TRUE)
#lines(mu.v2.hist$mids, mu.v2.hist$density)
#hist(W[3000:num.samples], breaks=40, prob=TRUE)
#lines(W.v2.hist$mids, W.v2.hist$density)
#hist(z[3000:num.samples,50], breaks=40, prob=TRUE)
#lines(z.v2.hist$mids, z.v2.hist$density)
