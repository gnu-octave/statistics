%!demo
%! ## Create a Nakagami distribution with default parameters
%! data = nakarnd (1, 1, 10000, 1);
%! pd = fitdist (data, "Nakagami");
%!
%! ## Query parameter 'mu' (shape parameter)
%! pd.mu
%!
%! ## Set parameter 'mu'
%! pd.mu = 2
%!
%! ## Use this to initialize or modify the shape parameter of a Nakagami
%! ## distribution. The shape parameter must be a real scalar >= 0.5, controlling
%! ## the fading severity in signal modeling; higher mu indicates less fading.

%!demo
%! ## Create a Nakagami distribution object by calling its constructor
%! pd = NakagamiDistribution (1.5, 2)
%!
%! ## Query parameter 'mu'
%! pd.mu
%!
%! ## This demonstrates direct construction with a specific shape parameter,
%! ## useful for modeling wireless channel fading with known characteristics.
