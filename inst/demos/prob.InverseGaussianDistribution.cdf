%!demo
%! ## Plot various CDFs from the Inverse Gaussian distribution
%! x = 0:0.01:5;
%! pd1 = makedist ("InverseGaussian", "mu", 1, "lambda", 1);
%! pd2 = makedist ("InverseGaussian", "mu", 1, "lambda", 2);
%! pd3 = makedist ("InverseGaussian", "mu", 1, "lambda", 3);
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"mu = 1, lambda = 1", "mu = 1, lambda = 2", "mu = 1, lambda = 3"}, ...
%!         "location", "southeast")
%! title ("Inverse Gaussian CDF")
%! xlabel ("values")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different Inverse Gaussian distributions, showing how probability
%! ## accumulates over values.
