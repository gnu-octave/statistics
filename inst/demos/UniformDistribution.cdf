%!demo
%! ## Plot various CDFs from the Uniform distribution
%! x = -1:0.01:12;
%! data1 = unifrnd (0, 5, 10000, 1);
%! data2 = unifrnd (2, 8, 10000, 1);
%! data3 = unifrnd (4, 10, 10000, 1);
%! pd1 = fitdist (data1, "Uniform");
%! pd2 = fitdist (data2, "Uniform");
%! pd3 = fitdist (data3, "Uniform");
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"Lower = 0, Upper = 5", "Lower = 2, Upper = 8", "Lower = 4, Upper = 10"}, ...
%!         "location", "southeast")
%! title ("Uniform CDF")
%! xlabel ("Values in x")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for Uniform distributions, showing how probability accumulates over the
%! ## defined interval, useful in probability assessments or simulations.
