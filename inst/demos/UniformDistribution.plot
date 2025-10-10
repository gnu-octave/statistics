%!demo
%! ## Create a Uniform distribution with fixed parameters Lower = 0 and Upper = 5
%! ## and plot its PDF.
%! pd = makedist ("Uniform", "Lower", 0, "Upper", 5)
%! plot (pd)
%! title ("Uniform distribution with Lower = 0 and Upper = 5")
%!
%! ## Use this to visualize the PDF of a Uniform distribution with fixed bounds,
%! ## useful for understanding the uniform probability density.

%!demo
%! ## Generate a data set of 100 random samples from a Uniform distribution
%! ## with parameters Lower = 0 and Upper = 5. Fit a Uniform distribution to this
%! ## data and plot its CDF superimposed over an empirical CDF.
%! rand ("seed", 21);
%! pd_fixed = makedist ("Uniform", "Lower", 0, "Upper", 5);
%! data = random (pd_fixed, 100, 1);
%! pd_fitted = fitdist (data, "Uniform");
%! plot (pd_fitted, "PlotType", "cdf")
%! txt = "Fitted Uniform distribution with Lower = %0.2f and Upper = %0.2f";
%! title (sprintf (txt, pd_fitted.Lower, pd_fitted.Upper))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## Use this to visualize the fitted CDF compared to the empirical CDF of the
%! ## data, useful for assessing model fit.
