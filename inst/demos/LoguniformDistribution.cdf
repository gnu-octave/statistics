%!demo
%! ## Plot various CDFs from the Log-uniform distribution
%! x = 0:0.01:10;
%! pd1 = LoguniformDistribution (1, 4);
%! pd2 = LoguniformDistribution (1, 6);
%! pd3 = LoguniformDistribution (1, 8);
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"Lower=1, Upper=4", "Lower=1, Upper=6", "Lower=1, Upper=8"}, ...
%!         "location", "southeast")
%! title ("Log-uniform CDF")
%! xlabel ("values in x (Lower <= x <= Upper)")
%! ylabel ("Cumulative probability")
%! 
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different Log-uniform distributions, showing how probability
%! ## accumulates over the range, useful in uncertainty modeling.
