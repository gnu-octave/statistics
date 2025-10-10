%!demo
%! ## Create a t Location-Scale distribution with fitted parameters
%! data = tlsrnd (0, 1, 5, 10000, 1);  % Generate data with mu=0, sigma=1, nu=5
%! pd = fitdist (data, "tLocationScale");
%!
%! ## Query parameter 'nu' (degrees of freedom)
%! pd.nu
%!
%! ## Set parameter 'nu'
%! pd.nu = 10
%!
%! ## Use this to initialize or modify the degrees of freedom, which controls the
%! ## tail heaviness of the t Location-Scale distribution. Nu must be a positive
%! ## real scalar, useful for modeling heavy-tailed data like stock returns.

%!demo
%! ## Create a t Location-Scale distribution object by calling its constructor
%! pd = tLocationScaleDistribution (0, 1, 3);
%!
%! ## Query parameter 'nu'
%! pd.nu
%!
%! ## This demonstrates setting the degrees of freedom directly via the constructor,
%! ## ideal for modeling data with specific tail behavior, such as outlier-prone datasets.
