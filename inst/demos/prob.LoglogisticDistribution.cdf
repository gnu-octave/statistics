%!demo
%! ## Plot various CDFs from the Log-logistic distribution
%! x = 0:0.01:10;
%! data1 = loglrnd (0, 0.5, 10000, 1);
%! data2 = loglrnd (0, 1, 10000, 1);
%! data3 = loglrnd (0, 2, 10000, 1);
%! pd1 = fitdist (data1, "Loglogistic");
%! pd2 = fitdist (data2, "Loglogistic");
%! pd3 = fitdist (data3, "Loglogistic");
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"mu = 0, sigma = 0.5", "mu = 0, sigma = 1", "mu = 0, sigma = 2"}, ...
%!         "location", "southeast")
%! title ("Log-logistic CDF")
%! xlabel ("values in x (x > 0)")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different Log-logistic distributions, showing how probability
%! ## accumulates for positive values, useful in survival analysis or risk modeling.
