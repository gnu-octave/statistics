%!demo
%! ## Create a Generalized Extreme Value distribution with fixed parameters k = 0,
%! ## sigma = 1, mu = 0 and plot its PDF.
%! pd = makedist ("GeneralizedExtremeValue", "k", 0, "sigma", 1, "mu", 0)
%! plot (pd)
%! title ("Fixed Generalized Extreme Value distribution with k = 0, sigma = 1, mu = 0")
%!
%! ## Use this to visualize the PDF of a GEV distribution with specified parameters.

%!demo
%! ## Generate a data set of 100 random samples from a Generalized Extreme Value
%! ## distribution with parameters k = 0, sigma = 1, mu = 0. Fit a GEV
%! ## distribution to this data and plot its CDF superimposed over an empirical CDF.
%! pd_fixed = makedist ("GeneralizedExtremeValue", "k", 0, "sigma", 1, "mu", 0)
%! rand ("seed", 21);
%! data = random (pd_fixed, 100, 1);
%! pd_fitted = fitdist (data, "GeneralizedExtremeValue")
%! plot (pd_fitted, "PlotType", "cdf")
%! txt = "Fitted Generalized Extreme Value distribution with k = %0.2f, sigma = %0.2f, mu = %0.2f";
%! title (sprintf (txt, pd_fitted.k, pd_fitted.sigma, pd_fitted.mu))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## Use this to visualize the fitted CDF compared to the empirical CDF of the
%! ## data, useful for assessing model fit.

%!demo
%! ## Generate a data set of 200 random samples from a Generalized Extreme Value
%! ## distribution with parameters k = 0, sigma = 1, mu = 0. Display a probability
%! ## plot for the GEV distribution fit to the data.
%! pd_fixed = makedist ("GeneralizedExtremeValue", "k", 0, "sigma", 1, "mu", 0)
%! rand ("seed", 21);
%! data = random (pd_fixed, 200, 1);
%! pd_fitted = fitdist (data, "GeneralizedExtremeValue")
%! plot (pd_fitted, "PlotType", "probability")
%! txt = strcat ("Probability plot of fitted Generalized Extreme Value", ...
%!               " distribution with k = %0.2f, sigma = %0.2f, mu = %0.2f");
%! title (sprintf (txt, pd_fitted.k, pd_fitted.sigma, pd_fitted.mu))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## This creates a probability plot to compare the fitted GEV distribution to
%! ## the data, useful for checking if the GEV model is appropriate.
