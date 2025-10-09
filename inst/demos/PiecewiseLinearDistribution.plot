%!demo
%! ## Create a Piecewise Linear distribution and plot its PDF.
%! load patients
%! [f, x] = ecdf (Weight);
%! pd = PiecewiseLinearDistribution (x(1:5:end), f(1:5:end));
%! plot (pd)
%! title ("Piecewise Linear distribution from Weight data")

%!demo
%! ## Create a Piecewise Linear distribution from data and plot its CDF.
%! load patients
%! [f, x] = ecdf (Weight);
%! pd = PiecewiseLinearDistribution (x(1:5:end), f(1:5:end));
%! plot (pd, "PlotType", "cdf")
%! title ("CDF of Piecewise Linear distribution from Weight data")
%!
%! ## Use this to visualize the CDF of the Piecewise Linear distribution,
%! ## useful for comparing to empirical CDF.
