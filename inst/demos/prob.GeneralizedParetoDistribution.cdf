%!demo
%! ## Plot various CDFs from the Generalized Pareto distribution
%! x = 0:0.01:5;
%! pd1 = makedist ("GeneralizedPareto", "k", 0.2, "sigma", 1, "theta", 0);
%! pd2 = makedist ("GeneralizedPareto", "k", 0.5, "sigma", 1, "theta", 0);
%! pd3 = makedist ("GeneralizedPareto", "k", 0.8, "sigma", 1, "theta", 0);
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"k = 0.2, sigma = 1, theta = 0", "k = 0.5, sigma = 1, theta = 0", "k = 0.8, sigma = 1, theta = 0"}, ...
%!         "location", "southeast")
%! title ("Generalized Pareto CDF")
%! xlabel ("Exceedance value")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different Generalized Pareto distributions, showing how probability
%! ## accumulates over exceedance values in extreme value modeling.
