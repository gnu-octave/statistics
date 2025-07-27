%!demo
%! ## Create a binomial distribution with default parameters
%! pd = makedist ("Binomial")
%!
%! ## Query parameter 'p' (probability of success)
%! pd.p
%!
%! ## Set parameter 'p'
%! pd.p = 0.4
%!
%! ## Use this to initialize or modify the probability of success in each
%! ## trial. The probability must be between 0 and 1.

%!demo
%! ## Create a binomial distribution object by calling its constructor
%! pd = BinomialDistribution (10, 0.3)
%!
%! ## Query parameter 'p'
%! pd.p
%!
%! ## This shows how to set the success probability directly via the
%! ## constructor, ideal for modeling specific success rates.
