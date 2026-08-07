%!demo
%! ## Plot various CDFs from the Extreme Value distribution
%! x = -5:0.01:10;
%! pd1 = makedist ("ExtremeValue", "mu", 0, "sigma", 0.5);
%! pd2 = makedist ("ExtremeValue", "mu", 0, "sigma", 1);
%! pd3 = makedist ("ExtremeValue", "mu", 0, "sigma", 1.5);
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"mu = 0, sigma = 0.5", "mu = 0, sigma = 1", "mu = 0, sigma = 1.5"}, ...
%!         "location", "southeast")
%! title ("Extreme Value CDF")
%! xlabel ("Values")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different Extreme Value distributions, showing how probability
%! ## accumulates over extreme values.
