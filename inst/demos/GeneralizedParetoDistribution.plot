%!demo
%! ## Create a Generalized Pareto distribution with fixed parameters k = 0.5, sigma = 1, theta = 0
%! ## and plot its PDF.
%!
%! pd = makedist ("GeneralizedPareto", "k", 0.3, "sigma", 1, "theta", 0)
%! plot (pd)
%! title ("Fixed Generalized Pareto distribution with k = 0.3, sigma = 1, theta = 0")

%!demo
%! ## Generate a data set of 100 random samples from a Generalized Pareto
%! ## distribution with parameters k = 0.5, sigma = 1, theta = 0. Fit a Generalized Pareto
%! ## distribution to this data (fixing theta=0) and plot its CDF superimposed over an empirical
%! ## CDF.
%!
%! pd_fixed = makedist ("GeneralizedPareto", "k", 0.5, "sigma", 1, "theta", 0)
%! rand ("seed", 21);
%! data = random (pd_fixed, 100, 1);
%! pd_fitted = GeneralizedParetoDistribution.fit (data, 0)
%! plot (pd_fitted, "PlotType", "cdf")
%! txt = "Fitted Generalized Pareto distribution with k = %0.2f, sigma = %0.2f, theta = 0";
%! title (sprintf (txt, pd_fitted.k, pd_fitted.sigma))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## Use this to visualize the fitted CDF compared to the empirical CDF of the
%! ## data, useful for assessing model fit in extreme value analysis.

%!demo
%! ## Generate a data set of 200 random samples from a Generalized Pareto
%! ## distribution with parameters k = 0.5, sigma = 1, theta = 0. Display a probability
%! ## plot for the Generalized Pareto distribution fit to the data (fixing theta=0).
%!
%! pd_fixed = makedist ("GeneralizedPareto", "k", 0.5, "sigma", 1, "theta", 0)
%! rand ("seed", 21);
%! data = random (pd_fixed, 200, 1);
%! pd_fitted = GeneralizedParetoDistribution.fit (data, 0)
%! plot (pd_fitted, "PlotType", "probability")
%! txt = strcat ("Probability plot of fitted Generalized Pareto", ...
%!               " distribution with k = %0.2f, sigma = %0.2f, theta = 0");
%! title (sprintf (txt, pd_fitted.k, pd_fitted.sigma))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## This creates a probability plot to compare the fitted distribution to the
%! ## data, useful for checking if the Generalized Pareto model is appropriate for tails.
