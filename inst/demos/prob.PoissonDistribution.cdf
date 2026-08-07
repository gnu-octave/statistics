%!demo
%! ## Plot various CDFs from the Poisson distribution
%! x = 0:20;
%! data1 = poissrnd (2, 10000, 1);
%! data2 = poissrnd (5, 10000, 1);
%! data3 = poissrnd (10, 10000, 1);
%! pd1 = fitdist (data1, "Poisson");
%! pd2 = fitdist (data2, "Poisson");
%! pd3 = fitdist (data3, "Poisson");
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"lambda = 2", "lambda = 5", "lambda = 10"}, "location", "southeast")
%! title ("Poisson CDF")
%! xlabel ("values in x (non-negative integers)")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different Poisson distributions, showing the probability of observing
%! ## at most k events, useful in risk assessment or queueing theory.
