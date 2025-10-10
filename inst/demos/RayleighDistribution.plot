%!demo
%! ## Create a Rayleigh distribution with fixed parameter sigma = 1 and plot its PDF.
%!
%! data = raylrnd (1, 10000, 1);
%! pd = fitdist (data, "Rayleigh");
%! plot (pd)
%! title ("Fixed Rayleigh distribution with sigma = 1")

%!demo
%! ## Generate a data set of 100 random samples from a Rayleigh
%! ## distribution with parameter sigma = 1. Fit a Rayleigh
%! ## distribution to this data and plot its CDF superimposed over an empirical
%! ## CDF.
%!
%! rand ("seed", 21);
%! data = raylrnd (1, 100, 1);
%! pd_fitted = fitdist (data, "Rayleigh");
%! plot (pd_fitted, "PlotType", "cdf")
%! txt = "Fitted Rayleigh distribution with sigma = %0.2f";
%! title (sprintf (txt, pd_fitted.sigma))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## Use this to visualize the fitted CDF compared to the empirical CDF of the
%! ## data, useful for assessing model fit.

%!demo
%! ## Generate a data set of 200 random samples from a Rayleigh
%! ## distribution with parameter sigma = 1. Display a probability
%! ## plot for the Rayleigh distribution fit to the data.
%!
%! rand ("seed", 21);
%! data = raylrnd (1, 200, 1);
%! pd_fitted = fitdist (data, "Rayleigh");
%! plot (pd_fitted, "PlotType", "probability")
%! txt = strcat ("Probability plot of fitted Rayleigh", ...
%!               " distribution with sigma = %0.2f");
%! title (sprintf (txt, pd_fitted.sigma))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## This creates a probability plot to compare the fitted distribution to the
%! ## data, useful for checking if the Rayleigh model is appropriate.
