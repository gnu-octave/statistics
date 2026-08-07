%!demo
%! ## Plot various CDFs from the Generalized Extreme Value distribution
%! x = -2:0.01:5;
%! pd1 = makedist ("GeneralizedExtremeValue", "k", 0, "sigma", 1, "mu", 0);
%! pd2 = makedist ("GeneralizedExtremeValue", "k", 0.2, "sigma", 1, "mu", 0);
%! pd3 = makedist ("GeneralizedExtremeValue", "k", -0.2, "sigma", 1, "mu", 0);
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"k = 0, sigma = 1, mu = 0", "k = 0.2, sigma = 1, mu = 0", "k = -0.2, sigma = 1, mu = 0"}, ...
%!         "location", "southeast")
%! title ("Generalized Extreme Value CDF")
%! xlabel ("Value")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function for
%! ## different GEV distributions, showing how probability accumulates over extreme
%! ## values with varying shape parameters.
