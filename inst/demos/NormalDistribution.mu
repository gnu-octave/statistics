%!demo
%! ## Create a Normal distribution with default parameters
%! data = randn (10000, 1);
%! pd = fitdist (data, "Normal");
%!
%! ## Query parameter 'mu' (location parameter, mean)
%! pd.mu
%!
%! ## Set parameter 'mu'
%! pd.mu = 1
%!
%! ## Use this to initialize or modify the mean of a Normal distribution.
%! ## The mean parameter must be a real scalar, representing the center of symmetry,
%! ## useful for shifting the distribution, such as modeling centered data like errors.

%!demo
%! ## Create a Normal distribution object by calling its constructor
%! pd = NormalDistribution (1.5, 2)
%!
%! ## Query parameter 'mu'
%! pd.mu
%!
%! ## This demonstrates direct construction with a specific mean,
%! ## suitable for modeling data with a known central tendency, like IQ scores or heights.
