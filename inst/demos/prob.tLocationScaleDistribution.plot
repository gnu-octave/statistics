%!demo
%! ## Create a t Location-Scale distribution with fixed parameters and plot its PDF
%! pd = tLocationScaleDistribution (0, 1, 5);
%! plot (pd)
%! title ("t Location-Scale distribution with mu = 0, sigma = 1, nu = 5")
%!
%! ## Use this to visualize the PDF of a t Location-Scale distribution with
%! ## fixed parameters, helpful for theoretical exploration.

%!demo
%! ## Generate a data set and plot the CDF of a fitted t Location-Scale distribution
%! rand ("seed", 21);
%! data = tlsrnd (0, 1, 5, 100, 1);
%! pd_fitted = fitdist (data, "tLocationScale");
%! plot (pd_fitted, "PlotType", "cdf")
%! txt = "Fitted t Location-Scale distribution with mu = %0.2f, sigma = %0.2f, nu = %0.2f";
%! title (sprintf (txt, pd_fitted.mu, pd_fitted.sigma, pd_fitted.nu))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## Use this to visualize the fitted CDF compared to the empirical CDF of the
%! ## data, useful for assessing model fit.

%!demo
%! ## Generate a data set and display a probability plot for a fitted t Location-Scale distribution
%! rand ("seed", 21);
%! data = tlsrnd (0, 1, 5, 200, 1);
%! pd_fitted = fitdist (data, "tLocationScale");
%! plot (pd_fitted, "PlotType", "probability")
%! txt = strcat ("Probability plot of fitted t Location-Scale", ...
%!               " distribution with mu = %0.2f, sigma = %0.2f, nu = %0.2f");
%! title (sprintf (txt, pd_fitted.mu, pd_fitted.sigma, pd_fitted.nu))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## This creates a probability plot to compare the fitted distribution to the
%! ## data, useful for checking if the t Location-Scale model is appropriate.
