%!demo
%! ## Create a Rician distribution by fitting to data
%! data = ricernd (1, 1, [10000, 1]);  % Generate data with s=1, sigma=1
%! pd = fitdist (data, "Rician");
%!
%! ## Query parameter 's' (noncentrality parameter)
%! pd.s
%!
%! ## Set parameter 's'
%! pd.s = 1.5
%!
%! ## Use this to initialize or modify the noncentrality parameter of a Rician
%! ## distribution. The noncentrality parameter 's' must be a non-negative real
%! ## scalar, representing the magnitude of the signal in the presence of noise.

%!demo
%! ## Create a Rician distribution object by calling its constructor
%! pd = RicianDistribution (2, 1)
%!
%! ## Query parameter 's'
%! pd.s
%!
%! ## This demonstrates direct construction with a specific noncentrality
%! ## parameter, useful for modeling data with a known signal strength.
