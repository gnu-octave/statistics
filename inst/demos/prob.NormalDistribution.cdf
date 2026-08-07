%!demo
%! ## Plot various CDFs from the Normal distribution
%! x = -5:0.01:5;
%! data1 = 0 + 0.5 * randn (10000, 1);
%! data2 = 0 + 1.0 * randn (10000, 1);
%! data3 = 0 + 2.0 * randn (10000, 1);
%! pd1 = fitdist (data1, "Normal");
%! pd2 = fitdist (data2, "Normal");
%! pd3 = fitdist (data3, "Normal");
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"mu = 0, sigma = 0.5", "mu = 0, sigma = 1", "mu = 0, sigma = 2"}, ...
%!         "location", "southeast")
%! title ("Normal CDF")
%! xlabel ("values in x")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different Normal distributions, showing how probability accumulates,
%! ## essential in hypothesis testing or confidence intervals.
