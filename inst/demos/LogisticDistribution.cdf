%!demo
%! ## Plot various CDFs from the Logistic distribution
%! x = -5:0.01:5;
%! pd1 = makedist ("LogisticDistribution", "mu", 0, "sigma", 0.5);
%! pd2 = makedist ("LogisticDistribution", "mu", 0, "sigma", 1);
%! pd3 = makedist ("LogisticDistribution", "mu", 0, "sigma", 1.5);
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"mu = 0, sigma = 0.5", "mu = 0, sigma = 1", "mu = 0, sigma = 1.5"}, ...
%!         "location", "southeast")
%! title ("Logistic CDF")
%! xlabel ("Value")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different Logistic distributions, showing how probability accumulates
%! ## over values, useful in regression or classification tasks.
