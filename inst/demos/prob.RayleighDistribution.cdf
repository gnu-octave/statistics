%!demo
%! ## Plot various CDFs from the Rayleigh distribution
%! x = -1:0.01:5;
%! data1 = raylrnd (0.5, 10000, 1);
%! data2 = raylrnd (1.0, 10000, 1);
%! data3 = raylrnd (2.0, 10000, 1);
%! pd1 = fitdist (data1, "Rayleigh");
%! pd2 = fitdist (data2, "Rayleigh");
%! pd3 = fitdist (data3, "Rayleigh");
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"sigma = 0.5", "sigma = 1", "sigma = 2"}, ...
%!         "location", "southeast")
%! title ("Rayleigh CDF")
%! xlabel ("values in x (x >= 0)")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different Rayleigh distributions, showing how probability
%! ## accumulates for nonnegative values, useful in signal processing or physics.
