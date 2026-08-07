%!demo
%! ## Create a Rician distribution with fixed parameters s = 1 and sigma = 1
%! ## and plot its PDF.
%! pd = RicianDistribution (1, 1);
%! plot (pd)
%! title ("Rician distribution with s = 1 and sigma = 1")
%!
%! ## Use this to visualize the PDF of a Rician distribution with fixed parameters,
%! ## useful for understanding the shape of the distribution.

%!demo
%! ## Generate a data set of 100 random samples from a Rician distribution
%! ## with parameters s = 1 and sigma = 1. Fit a Rician distribution to this
%! ## data and plot its CDF superimposed over an empirical CDF.
%! rand ("seed", 21);
%! data = ricernd (1, 1, [100, 1]);
%! pd_fitted = fitdist (data, "Rician");
%! plot (pd_fitted, "PlotType", "cdf")
%! txt = "Fitted Rician distribution with s = %0.2f and sigma = %0.2f";
%! title (sprintf (txt, pd_fitted.s, pd_fitted.sigma))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## Use this to visualize the fitted CDF compared to the empirical CDF of the
%! ## data, useful for assessing model fit.

%!demo
%! ## Generate a data set of 200 random samples from a Rician distribution
%! ## with parameters s = 1 and sigma = 1. Display a probability plot for the
%! ## Rician distribution fit to the data.
%! rand ("seed", 21);
%! data = ricernd (1, 1, [200, 1]);
%! pd_fitted = fitdist (data, "Rician");
%! plot (pd_fitted, "PlotType", "probability")
%! txt = strcat ("Probability plot of fitted Rician distribution", ...
%!               " with s = %0.2f and sigma = %0.2f");
%! title (sprintf (txt, pd_fitted.s, pd_fitted.sigma))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## This creates a probability plot to compare the fitted distribution to the
%! ## data, useful for checking if the Rician model is appropriate.
