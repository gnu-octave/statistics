%!demo
%! ## Plot various CDFs from the Birnbaum-Saunders distribution
%! x = 0:0.1:5;
%! pd1 = makedist ("BirnbaumSaunders", "beta", 1, "gamma", 0.2);
%! pd2 = makedist ("BirnbaumSaunders", "beta", 1, "gamma", 0.5);
%! pd3 = makedist ("BirnbaumSaunders", "beta", 1, "gamma", 0.8);
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r", "LineWidth", 2)
%! grid on
%! legend ({"beta = 1, gamma = 0.2", "beta = 1, gamma = 0.5", "beta = 1, gamma = 0.8"}, ...
%!         "location", "southeast")
%! title ("Birnbaum-Saunders CDF")
%! xlabel ("Time to failure")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different Birnbaum-Saunders distributions, showing how probability
%! ## accumulates over time-to-failure.
