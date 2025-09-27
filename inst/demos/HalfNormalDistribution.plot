%!demo
%! ## Create a Half-normal distribution with fixed parameters mu = 0 and
%! ## sigma = 1 and plot its PDF.
%!
%! data = abs (randn (10000, 1));
%! pd = fitdist (data, "HalfNormal");
%! plot (pd)
%! title ("Fixed Half-normal distribution with mu = 0 and sigma = 1")

%!demo
%! ## Generate a data set of 100 random samples from a Half-normal
%! ## distribution with parameters mu = 0 and sigma = 1. Fit a Half-normal
%! ## distribution to this data and plot its CDF superimposed over an empirical
%! ## CDF.
%!
%! rand ("seed", 21);
%! data = abs (randn (100, 1));
%! pd_fitted = fitdist (data, "HalfNormal");
%! plot (pd_fitted, "PlotType", "cdf")
%! txt = "Fitted Half-normal distribution with mu = %0.2f and sigma = %0.2f";
%! title (sprintf (txt, pd_fitted.mu, pd_fitted.sigma))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## Use this to visualize the fitted CDF compared to the empirical CDF of the
%! ## data, useful for assessing model fit.

%!demo
%! ## Generate a data set of 200 random samples from a Half-normal
%! ## distribution with parameters mu = 0 and sigma = 1. Display a probability
%! ## plot for the Half-normal distribution fit to the data.
%!
%! rand ("seed", 21);
%! data = abs (randn (200, 1));
%! pd_fitted = fitdist (data, "HalfNormal");
%! plot (pd_fitted, "PlotType", "probability")
%! txt = strcat ("Probability plot of fitted Half-normal", ...
%!               " distribution with mu = %0.2f and sigma = %0.2f");
%! title (sprintf (txt, pd_fitted.mu, pd_fitted.sigma))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## This creates a probability plot to compare the fitted distribution to the
%! ## data, useful for checking if the Half-normal model is appropriate.
