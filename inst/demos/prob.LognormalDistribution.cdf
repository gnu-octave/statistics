%!demo
%! ## Plot various CDFs from the Lognormal distribution
%! x = 0:0.01:10;
%! data1 = lognrnd (0, 0.5, 10000, 1);
%! data2 = lognrnd (0, 1.0, 10000, 1);
%! data3 = lognrnd (0, 1.5, 10000, 1);
%! pd1 = fitdist (data1, "Lognormal");
%! pd2 = fitdist (data2, "Lognormal");
%! pd3 = fitdist (data3, "Lognormal");
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"mu = 0, sigma = 0.5", "mu = 0, sigma = 1", "mu = 0, sigma = 1.5"}, ...
%!         "location", "southeast")
%! title ("Lognormal CDF")
%! xlabel ("values in x (x > 0)")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different Lognormal distributions, showing how probability accumulates
%! ## for positive skewed data, useful in finance or biology modeling.
