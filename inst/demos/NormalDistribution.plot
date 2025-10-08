%!demo
%! ## Create a Normal distribution with fixed parameters mu = 0 and
%! ## sigma = 1 and plot its PDF.
%!
%! data = randn (10000, 1);
%! pd = fitdist (data, "Normal");
%! plot (pd)
%! title ("Fixed Normal distribution with mu = 0 and sigma = 1")

%!demo
%! ## Generate a data set of 100 random samples from a Normal
%! ## distribution with parameters mu = 0 and sigma = 1. Fit a Normal
%! ## distribution to this data and plot its CDF superimposed over an empirical
%! ## CDF.
%!
%! rand ("seed", 21);
%! data = randn (100, 1);
%! pd_fitted = fitdist (data, "Normal");
%! plot (pd_fitted, "PlotType", "cdf")
%! txt = "Fitted Normal distribution with mu = %0.2f and sigma = %0.2f";
%! title (sprintf (txt, pd_fitted.mu, pd_fitted.sigma))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## Use this to visualize the fitted CDF compared to the empirical CDF of the
%! ## data, useful for assessing model fit in symmetric distributions.

%!demo
%! ## Generate a data set of 200 random samples from a Normal
%! ## distribution with parameters mu = 0 and sigma = 1. Display a probability
%! ## plot for the Normal distribution fit to the data.
%!
%! rand ("seed", 21);
%! data = randn (200, 1);
%! pd_fitted = fitdist (data, "Normal");
%! plot (pd_fitted, "PlotType", "probability")
%! txt = strcat ("Probability plot of fitted Normal", ...
%!               " distribution with mu = %0.2f and sigma = %0.2f");
%! title (sprintf (txt, pd_fitted.mu, pd_fitted.sigma))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## This creates a probability plot to compare the fitted distribution to the
%! ## data, useful for checking normality assumptions.
