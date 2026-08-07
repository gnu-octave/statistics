%!demo
%! ## Create a Birnbaum-Saunders distribution with fixed parameters β = 1 and
%! ## γ = 0.5 and plot its PDF.
%!
%! pd = makedist ("BirnbaumSaunders", "beta", 1, "gamma", 0.5)
%! plot (pd)
%! title ("Fixed Birnbaum-Saunders distribution with beta = 1 and gamma = 0.5")

%!demo
%! ## Generate a data set of 100 random samples from a Birnbaum-Saunders
%! ## distribution with parameters β = 1 and γ = 0.5. Fit a Birnbaum-Saunders
%! ## distribution to this data and plot its CDF superimposed over an empirical
%! ## CDF.
%!
%! pd_fixed = makedist ("BirnbaumSaunders", "beta", 1, "gamma", 0.5)
%! randg ("seed", 21);
%! data = random (pd_fixed, 100, 1);
%! pd_fitted = fitdist (data, "BirnbaumSaunders")
%! plot (pd_fitted, "PlotType", "cdf")
%! txt = "Fitted Birnbaum-Saunders distribution with β = %0.2f and γ = %0.2f";
%! title (sprintf (txt, pd_fitted.beta, pd_fitted.gamma))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## Use this to visualize the fitted CDF compared to the empirical CDF of the
%! ## data, useful for assessing model fit.

%!demo
%! ## Generate a data set of 200 random samples from a Birnbaum-Saunders
%! ## distribution with parameters β = 1 and γ = 0.5. Display a probability
%! ## plot for the Birnbaum-Saunders distribution fit to the data.
%!
%! pd_fixed = makedist ("BirnbaumSaunders", "beta", 1, "gamma", 0.5)
%! randg ("seed", 21);
%! data = random (pd_fixed, 200, 1);
%! pd_fitted = fitdist (data, "BirnbaumSaunders")
%! plot (pd_fitted, "PlotType", "probability")
%! txt = strcat ("Probability plot of fitted Birnbaum-Saunders", ...
%!               " distribution with β = %0.2f and γ = %0.2f");
%! title (sprintf (txt, pd_fitted.beta, pd_fitted.gamma))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## This creates a probability plot to compare the fitted distribution to the
%! ## data, useful for checking if the Birnbaum-Saunders model is appropriate.
