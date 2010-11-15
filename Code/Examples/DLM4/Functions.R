# A collection of useful functions.

# A two dimensional histogram: The data we are counting is x and y.
# The histogram bins things according to xbins and ybins.  xbins is a
# sequence that defines how we define the intervals characterizing the
# bins for the x values.  For instance, if xbins=c(0, 0.5, 1), then
# ths bins are (basically) (0, 0.5] and (0.5, 1].  We handle values
# outside the bins by putting them in the first or last most bin.
hist.2d = function(x, y, xbins, ybins){
  # Set up the data structure for the histogram.
  n.xbins = length(xbins)-1;
  n.ybins = length(ybins)-1;
  the.hist = matrix(0, n.xbins, n.ybins);
  # Get the number of samples.
  N = length(x);
  # Calculate the lengths of the intervals.
  len.x = xbins[n.xbins]-xbins[1];
  len.y = ybins[n.ybins]-ybins[1];
  # Now fill the histogram.
  for(i in 1:N){
    ratio.x = (x[i]-xbins[1])/len.x;
    idx.x = ceiling(n.xbins*ratio.x);
    ratio.y = (y[i]-ybins[1])/len.y;
    idx.y = ceiling(n.ybins*ratio.y);
    idx.x = min(n.xbins, max(1, idx.x));
    idx.y = min(n.ybins, max(1, idx.y));
    the.hist[idx.x, idx.y] = the.hist[idx.x, idx.y] + 1;
  }
  the.hist;
}

#hist.2d = function(x, y, xlim, ylim, xbin, ybin){
#  # Set up the data structure for the histogram.
#  the.hist = matrix(0, xbin, ybin);
#  # Get the number of samples.
#  N = length(x);
#  # Calculate the lengths of the intervals.
#  len.x = xlim[2]-xlim[1];
#  len.y = ylim[2]-ylim[1];
#  # Now fill the histogram.
#  for(i in 1:N){
#    ratio.x = x[i]/len.x;
#    idx.x = ceiling(xbin*ratio.x);
#    ratio.y = y[i]/len.y;
#    idx.y = ceiling(ybin*ratio.y);
#    idx.x = min(xbin, max(1, idx.x));
#    idx.y = min(ybin, max(1, idx.y));
#    the.hist[idx.x, idx.y] = the.hist[idx.x, idx.y] + 1;
#  }
#  the.hist;
#}
