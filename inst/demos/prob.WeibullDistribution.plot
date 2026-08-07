%!demo
%! ## Create a Weibull distribution with fixed parameters lambda = 1 and
%! ## k = 1 and plot its PDF.
%!
%! data = wblrnd (1, 1, 10000, 1);
%! pd = fitdist (data, "Weibull");
%! plot (pd)
%! title ("Fixed Weibull distribution with lambda = 1 and k = 1")

%!demo
%! ## Generate a data set of 100 random samples from a Weibull distribution
%! ## with parameters lambda = 1 and k = 1. Fit a Weibull distribution to
%! ## this data and plot its CDF superimposed over an empirical CDF.
%!
%! rand ("seed", 21);
%! data = wblrnd (1, 1, 100, 1);
%! pd_fitted = fitdist (data, "Weibull");
%! plot (pd_fitted, "PlotType", "cdf")
%! txt = "Fitted Weibull distribution with lambda = %0.2f and k = %0.2f";
%! title (sprintf (txt, pd_fitted.lambda, pd_fitted.k))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## Use this to visualize the fitted CDF compared to the empirical CDF of the
%! ## data, useful for assessing model fit in time-to-failure scenarios.

%!demo
%! ## Generate a data set of 200 random samples from a Weibull distribution
%! ## with parameters lambda = 1 and k = 1. Display a probability plot for
%! ## the Weibull distribution fit to the data.
%!
%! rand ("seed", 21);
%! data = wblrnd (1, 1, 200, 1);
%! pd_fitted = fitdist (data, "Weibull");
%! plot (pd_fitted, "PlotType", "probability")
%! txt = strcat ("Probability plot of fitted Weibull distribution", ...
%!               " with lambda = %0.2f and k = %0.2f");
%! title (sprintf (txt, pd_fitted.lambda, pd_fitted.k))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## This creates a probability plot to compare the fitted distribution to the
%! ## data, useful for validating the Weibull assumption in reliability studies.
