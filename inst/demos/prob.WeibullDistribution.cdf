%!demo
%! ## Plot various CDFs from the Weibull distribution
%! x = 0:0.01:5;
%! data1 = wblrnd (1, 0.5, 10000, 1);
%! data2 = wblrnd (1, 1, 10000, 1);
%! data3 = wblrnd (1, 2, 10000, 1);
%! pd1 = fitdist (data1, "Weibull");
%! pd2 = fitdist (data2, "Weibull");
%! pd3 = fitdist (data3, "Weibull");
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"lambda = 1, k = 0.5", "lambda = 1, k = 1", "lambda = 1, k = 2"}, ...
%!         "location", "southeast")
%! title ("Weibull CDF")
%! xlabel ("values in x (x >= 0)")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different Weibull distributions, showing how probability accumulates
%! ## for positive values, useful in survival analysis or time-to-event modeling.
