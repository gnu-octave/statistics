%!demo
%! ## Create a t Location-Scale distribution by fitting to data
%! data = tlsrnd (0, 1, 5, 10000, 1);  % Generate data with mu=0, sigma=1, nu=5
%! pd = fitdist (data, "tLocationScale");
%!
%! ## Query parameter 'mu' (location parameter)
%! pd.mu
%!
%! ## Set parameter 'mu'
%! pd.mu = 1
%!
%! ## Use this to initialize or modify the location parameter of a t Location-Scale
%! ## distribution. The location parameter (mu) is a real scalar that shifts the
%! ## distribution, useful for modeling data centered around a specific value.

%!demo
%! ## Create a t Location-Scale distribution object by calling its constructor
%! pd = tLocationScaleDistribution (2, 1, 5);
%!
%! ## Query parameter 'mu'
%! pd.mu
%!
%! ## This demonstrates direct construction with a specific location parameter,
%! ## ideal for modeling data with a known center, such as test scores or residuals.
