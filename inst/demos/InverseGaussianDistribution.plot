%!demo
%! ## Create an Inverse Gaussian distribution with fixed parameters μ = 1 and
%! ## λ = 2 and plot its PDF.
%!
%! pd = makedist ("InverseGaussian", "mu", 1, "lambda", 2)
%! plot (pd)
%! title ("Fixed Inverse Gaussian distribution with mu = 1 and lambda = 2")

%!demo
%! ## Generate a data set of 100 random samples from an Inverse Gaussian
%! ## distribution with parameters μ = 1 and λ = 2. Fit an Inverse Gaussian
%! ## distribution to this data and plot its CDF superimposed over an empirical
%! ## CDF.
%!
%! pd_fixed = makedist ("InverseGaussian", "mu", 1, "lambda", 2)
%! rand ("seed", 21);
%! data = random (pd_fixed, 100, 1);
%! pd_fitted = fitdist (data, "InverseGaussian")
%! plot (pd_fitted, "PlotType", "cdf")
%! txt = "Fitted Inverse Gaussian distribution with μ = %0.2f and λ = %0.2f";
%! title (sprintf (txt, pd_fitted.mu, pd_fitted.lambda))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## Use this to visualize the fitted CDF compared to the empirical CDF of the
%! ## data, useful for assessing model fit.

%!demo
%! ## Generate a data set of 200 random samples from an Inverse Gaussian
%! ## distribution with parameters μ = 1 and λ = 2. Display a probability
%! ## plot for the Inverse Gaussian distribution fit to the data.
%!
%! pd_fixed = makedist ("InverseGaussian", "mu", 1, "lambda", 2)
%! rand ("seed", 21);
%! data = random (pd_fixed, 200, 1);
%! pd_fitted = fitdist (data, "InverseGaussian")
%! plot (pd_fitted, "PlotType", "probability")
%! txt = strcat ("Probability plot of fitted Inverse Gaussian", ...
%!               " distribution with μ = %0.2f and λ = %0.2f");
%! title (sprintf (txt, pd_fitted.mu, pd_fitted.lambda))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## This creates a probability plot to compare the fitted distribution to the
%! ## data, useful for checking if the Inverse Gaussian model is appropriate.
