%!demo
%! ## Plot various CDFs from the Piecewise Linear distribution
%! load patients
%! [f1, x1] = ecdf (Weight);
%! [f2, x2] = ecdf (Height);
%! pd1 = PiecewiseLinearDistribution (x1(1:10:end), f1(1:10:end));
%! pd2 = PiecewiseLinearDistribution (x2(1:10:end), f2(1:10:end));
%! vals = 50:0.1:250;
%! p1 = cdf (pd1, vals);
%! p2 = cdf (pd2, vals);
%! plot (vals, p1, "-b", vals, p2, "-r")
%! grid on
%! legend ({"Weight", "Height"}, "location", "southeast")
%! title ("Piecewise Linear CDF")
%! xlabel ("values in x")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different Piecewise Linear distributions, showing how probability
%! ## accumulates across the defined segments.
