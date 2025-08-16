%!demo
%! ## Plot various CDFs from the Exponential distribution
%! x = 0:0.01:10;
%! pd1 = makedist ("Exponential", "mu", 1);
%! pd2 = makedist ("Exponential", "mu", 2);
%! pd3 = makedist ("Exponential", "mu", 3);
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"mu = 1", "mu = 2", "mu = 3"}, "location", "southeast")
%! title ("Exponential CDF")
%! xlabel ("values in x")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different Exponential distributions, showing how probability
%! ## accumulates over values.
