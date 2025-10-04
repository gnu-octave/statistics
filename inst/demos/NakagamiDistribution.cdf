%!demo
%! ## Plot various CDFs from the Nakagami distribution
%! x = -1:0.01:5;
%! data1 = nakarnd (0.5, 1, 10000, 1);
%! data2 = nakarnd (1, 1, 10000, 1);
%! data3 = nakarnd (2, 1, 10000, 1);
%! pd1 = fitdist (data1, "Nakagami");
%! pd2 = fitdist (data2, "Nakagami");
%! pd3 = fitdist (data3, "Nakagami");
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"mu = 0.5, omega = 1", "mu = 1, omega = 1", "mu = 2, omega = 1"}, ...
%!         "location", "southeast")
%! title ("Nakagami CDF")
%! xlabel ("values in x (x >= 0)")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different Nakagami distributions, showing how probability accumulates
%! ## for signal amplitudes, useful in wireless communications analysis.
