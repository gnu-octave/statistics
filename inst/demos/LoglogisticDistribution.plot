%!demo
%! ## Create a Log-logistic distribution with fixed parameters mu = 0 and
%! ## sigma = 1 and plot its PDF.
%!
%! data = loglrnd (0, 1, 10000, 1);
%! pd = fitdist (data, "Loglogistic");
%! plot (pd)
%! title ("Fixed Log-logistic distribution with mu = 0 and sigma = 1")

%!demo
%! ## Generate a data set of 100 random samples from a Log-logistic
%! ## distribution with parameters mu = 0 and sigma = 1. Fit a Log-logistic
%! ## distribution to this data and plot its CDF superimposed over an empirical
%! ## CDF.
%!
%! rand ("seed", 21);
%! data = loglrnd (0, 1, 100, 1);
%! pd_fitted = fitdist (data, "Loglogistic");
%! plot (pd_fitted, "PlotType", "cdf")
%! txt = "Fitted Log-logistic distribution with mu = %0.2f and sigma = %0.2f";
%! title (sprintf (txt, pd_fitted.mu, pd_fitted.sigma))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## Use this to visualize the fitted CDF compared to the empirical CDF of the
%! ## data, useful for assessing model fit in skewed positive data.

%!demo
%! ## Generate a data set of 200 random samples from a Log-logistic
%! ## distribution with parameters mu = 0 and sigma = 1. Display a probability
%! ## plot for the Log-logistic distribution fit to the data.
%!
%! rand ("seed", 21);
%! data = loglrnd (0, 1, 200, 1);
%! pd_fitted = fitdist (data, "Loglogistic");
%! plot (pd_fitted, "PlotType", "probability")
%! txt = strcat ("Probability plot of fitted Log-logistic", ...
%!               " distribution with mu = %0.2f and sigma = %0.2f");
%! title (sprintf (txt, pd_fitted.mu, pd_fitted.sigma))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## This creates a probability plot to compare the fitted distribution to the
%! ## data, useful for checking if the Log-logistic model captures the tail behavior.
