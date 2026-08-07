%!demo
%! ## Create an Extreme Value distribution with fixed parameters mu = 0 and
%! ## sigma = 1 and plot its PDF.
%!
%! pd = makedist ("ExtremeValue", "mu", 0, "sigma", 1)
%! plot (pd)
%! title ("Fixed Extreme Value distribution with mu = 0 and sigma = 1")

%!demo
%! ## Generate a data set of 100 random samples from an Extreme Value
%! ## distribution with parameters mu = 0 and sigma = 1. Fit an Extreme Value
%! ## distribution to this data and plot its CDF superimposed over an empirical
%! ## CDF.
%!
%! pd_fixed = makedist ("ExtremeValue", "mu", 0, "sigma", 1)
%! rand ("seed", 21);
%! data = random (pd_fixed, 100, 1);
%! pd_fitted = fitdist (data, "ExtremeValue")
%! plot (pd_fitted, "PlotType", "cdf")
%! txt = "Fitted Extreme Value distribution with mu = %0.2f and sigma = %0.2f";
%! title (sprintf (txt, pd_fitted.mu, pd_fitted.sigma))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## Use this to visualize the fitted CDF compared to the empirical CDF of the
%! ## data, useful for assessing model fit.

%!demo
%! ## Generate a data set of 200 random samples from an Extreme Value
%! ## distribution with parameters mu = 0 and sigma = 1. Display a probability
%! ## plot for the Extreme Value distribution fit to the data.
%!
%! pd_fixed = makedist ("ExtremeValue", "mu", 0, "sigma", 1)
%! rand ("seed", 21);
%! data = random (pd_fixed, 200, 1);
%! pd_fitted = fitdist (data, "ExtremeValue")
%! plot (pd_fitted, "PlotType", "probability")
%! txt = strcat ("Probability plot of fitted Extreme Value", ...
%!               " distribution with mu = %0.2f and sigma = %0.2f");
%! title (sprintf (txt, pd_fitted.mu, pd_fitted.sigma))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## This creates a probability plot to compare the fitted distribution to the
%! ## data, useful for checking if the Extreme Value model is appropriate.
