%!demo
%! ## Create a Piecewise Linear distribution object by calling its constructor
%! pd = PiecewiseLinearDistribution ([0; 1; 2; 3], [0; 0.3; 0.7; 1])
%!
%! ## Query parameter 'Fx'
%! pd.Fx
%!
%! ## This shows how to set CDF values directly via the constructor, ideal for
%! ## defining piecewise linear approximations of arbitrary distributions.
