%!demo
%! ## Plot various CDFs from the t Location-Scale distribution
%! x = -5:0.01:5;
%! data1 = tlsrnd (0, 0.5, 5, 10000, 1);
%! data2 = tlsrnd (0, 1, 5, 10000, 1);
%! data3 = tlsrnd (0, 2, 5, 10000, 1);
%! pd1 = fitdist (data1, "tLocationScale");
%! pd2 = fitdist (data2, "tLocationScale");
%! pd3 = fitdist (data3, "tLocationScale");
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"mu = 0, sigma = 0.5, nu = 5", "mu = 0, sigma = 1, nu = 5", ...
%!         "mu = 0, sigma = 2, nu = 5"}, "location", "southeast")
%! title ("t Location-Scale CDF")
%! xlabel ("Values in x")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different t Location-Scale distributions, showing how probability
%! ## accumulates, useful in risk analysis or hypothesis testing.
