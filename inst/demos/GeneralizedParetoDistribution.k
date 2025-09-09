%!demo
%! ## Create a Generalized Pareto distribution with default parameters
%! pd = makedist ("GeneralizedPareto")
%!
%! ## Query parameter 'k' (shape parameter)
%! pd.k
%!
%! ## Set parameter 'k'
%! pd.k = 0.5
%!
%! ## Use this to initialize or modify the shape parameter of a Generalized Pareto
%! ## distribution. The shape parameter can be any real scalar and determines the
%! ## tail behavior (heavy-tailed if k>0).

%!demo
%! ## Create a Generalized Pareto distribution object by calling its constructor
%! pd = GeneralizedParetoDistribution (0.5, 1, 0)
%!
%! ## Query parameter 'k'
%! pd.k
%!
%! ## This demonstrates direct construction with a specific shape parameter,
%! ## useful for modeling extreme values with known tail index.
