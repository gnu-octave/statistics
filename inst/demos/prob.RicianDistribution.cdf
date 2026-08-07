%!demo
%! ## Plot various CDFs from the Rician distribution
%! x = 0:0.01:5;
%! data1 = ricernd (1, 0.5, [10000, 1]);
%! data2 = ricernd (1, 1, [10000, 1]);
%! data3 = ricernd (1, 2, [10000, 1]);
%! pd1 = fitdist (data1, "Rician");
%! pd2 = fitdist (data2, "Rician");
%! pd3 = fitdist (data3, "Rician");
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"s = 1, sigma = 0.5", "s = 1, sigma = 1", "s = 1, sigma = 2"}, ...
%!         "location", "southeast")
%! title ("Rician CDF")
%! xlabel ("values in x (x >= 0)")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different Rician distributions, showing how probability accumulates
%! ## for non-negative signal magnitudes, useful in signal processing.
