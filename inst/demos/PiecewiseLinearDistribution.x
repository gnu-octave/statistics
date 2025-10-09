%!demo
%! ## Create a Piecewise Linear distribution object by calling its constructor
%! pd = PiecewiseLinearDistribution ([0; 1; 2; 3], [0; 0.3; 0.7; 1])
%!
%! ## Query parameter 'x'
%! pd.x
%!
%! ## This demonstrates direct construction with specific x values, useful for
%! ## modeling custom empirical distributions or interpolating between known points.
