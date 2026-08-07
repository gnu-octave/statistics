%!demo
%! ## Create a Negative Binomial distribution with default parameters
%! data = nbinrnd(5, 0.5, 10000, 1);
%! pd = fitdist (data, "NegativeBinomial");
%!
%! ## Query parameter 'R' (number of successes)
%! pd.R
%!
%! ## Set parameter 'R'
%! pd.R = 10
%!
%! ## Use this to initialize or modify the number of successes parameter in a
%! ## Negative Binomial distribution. R must be a positive scalar, controlling
%! ## the shape and often fixed based on the problem context, like successes in trials.

%!demo
%! ## Create a Negative Binomial distribution object by calling its constructor
%! pd = NegativeBinomialDistribution(10, 0.3)
%!
%! ## Query parameter 'R'
%! pd.R
%!
%! ## This demonstrates direct construction with a specific number of successes,
%! ## useful for modeling scenarios with a known fixed number of events, such as
%! ## in reliability testing or count processes.
