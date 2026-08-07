%!demo
%! ## Plot various CDFs from the Half-normal distribution
%! x = -1:0.01:5;
%! data1 = 0.5 * abs (randn (10000, 1));
%! data2 = 1.0 * abs (randn (10000, 1));
%! data3 = 2.0 * abs (randn (10000, 1));
%! pd1 = fitdist (data1, "HalfNormal");
%! pd2 = fitdist (data2, "HalfNormal");
%! pd3 = fitdist (data3, "HalfNormal");
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"mu = 0, sigma = 0.5", "mu = 0, sigma = 1", "mu = 0, sigma = 2"}, ...
%!         "location", "southeast")
%! title ("Half-normal CDF")
%! xlabel ("values in x (x >= mu)")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different Half-normal distributions, showing how probability
%! ## accumulates for positive values, useful in reliability or error modeling.
