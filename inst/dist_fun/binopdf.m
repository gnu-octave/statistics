## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2021 Nicholas R. Jankowski
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{y} =} binopdf (@var{x}, @var{n}, @var{ps})
##
## Binomial probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the binomial distribution with parameters @var{n} and @var{ps}, where
## @var{n} is the number of trials and @var{ps} is the probability of success.
## The size of @var{y} is the common size of @var{x}, @var{n}, and @var{ps}.  A
## scalar input functions as a constant matrix of the same size as the other
## inputs.
##
## Matlab incompatibility: Octave's @code{binopdf} does not allow complex
## input values.  Matlab 2021b returns values for complex inputs despite the
## documentation indicates integer and real value inputs are required.
##
## Further information about the binomial distribution can be found at
## @url{https://en.wikipedia.org/wiki/Binomial_distribution}
##
## @seealso{binocdf, binoinv, binornd, binofit, binolike, binostat, binotest}
## @end deftypefn

function y = binopdf (x, n, ps)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("binopdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, N, and PS
  if (! isscalar (x) || ! isscalar (n) || ! isscalar (ps))
    [retval, x, n, ps] = common_size (x, n, ps);
    if (retval > 0)
      error ("binopdf: X, N, and PS must be of common size or scalars.");
    endif
  endif

  ## Check for X, N, and PS being reals
  if (iscomplex (x) || iscomplex (n) || iscomplex (ps))
    error ("binopdf: X, N, and PS must not be complex.");
  endif

  sz_x = size (x); # save original size for reshape later
  x = x(:); n = n(:); ps = ps(:); # columns for easier vectorization

  ## Initialize output, preserve class of output if any are singles
  if (isa (x, "single") || isa (n, "single") || isa (ps, "single"));
    y = zeros (numel (x), 1, "single");
  else
    y = zeros (numel (x), 1);
  endif

  ## k - index of array locations needing calculation
  k = (x == fix (x)) & (n == fix (n)) & (n >= 0) & (ps >= 0) & (ps <= 1) ...
       & (x >= 0) & (x <= n);

  nx = n - x;
  q = 1 - ps;


  ## Catch special cases ahead of calculations:

  ## Matlab incompatibility: Matlab 2021b returns values for complex inputs
  ## despite documentation indicating integer and real value inputs required.
  ## Octave chooses to return an NaN instead.
  catch_special = (iscomplex (x) | iscomplex (n) | iscomplex (ps));
  k(catch_special) = false;
  y(catch_special) = NaN;

  ## x = 0 and x = n cases where ps != 0 or 1, respectively
  ## remove them from k, use alternate calculation to avoid /0
  catch_special = (x == 0)& (! catch_special);
  k(catch_special) = false;
  y(catch_special) = exp (n(catch_special) .* log (q(catch_special)));

  catch_special = (nx == 0) & (! catch_special);
  k(catch_special) = false;
  y(catch_special) = exp (n(catch_special) .* log (ps(catch_special)));


  ## Perform Loader pdf calculation on non-trivial elements
  if (any (k))
    y(k) = loader_expansion (x(k), n(k), ps(k), nx(k), q(k));
  endif

  ## Trivial case special outputs:
  ksp = ((ps == 0) & (x == 0)) | (ps == 1) & (x == n);
  y(ksp) = 1;

  ## Input NaN, n not pos int, or ps outside [0,1],
  ## set output to NaN (overrides 0 or 1)
  ksp = (n ~= fix (n)) | (n < 0) | (ps < 0) | (ps > 1) | isnan (x) ...
           | isnan (n) | isnan (ps);
  y(ksp) = NaN;

  y = reshape (y, sz_x); ## restore output to input shape
endfunction

function y = loader_expansion (x, n, ps, nx, q)
    ## Precalculated constants, d_n from n = 0 to 30
    ## extended from Loader using octave symbolic vpa
    ## out to n = 30
    d_n = [
      0.08106146679532725821967026359438236013860,
      0.04134069595540929409382208140711750802535,
      0.02767792568499833914878929274624466659538,
      0.02079067210376509311152277176784865633309,
      0.01664469118982119216319486537359339114739,
      0.01387612882307074799874572702376290856175,
      0.01189670994589177009505572411765943862013,
      0.01041126526197209649747856713253462919952,
      0.00925546218271273291772863663310013611743,
      0.00833056343336287125646931865962855220929,
      0.00757367548795184079497202421159508389293,
      0.00694284010720952986566415266347536265992,
      0.00640899418800420706843963108297831257520,
      0.00595137011275884773562441604646945832642,
      0.00555473355196280137103868995979228464907,
      0.00520765591960964044071799685790189865099,
      0.00490139594843473786071681819096755442865,
      0.00462915374933402859242721316419232323878,
      0.00438556024923232426828773634861946570116,
      0.00416631969199692245746292338221831613633,
      0.00396795421864085961728763680734281467287,
      0.00378761806844443457786667706893349200129,
      0.00362296022468309470738119836390285473489,
      0.00347202138297876696294511542270952959204,
      0.00333315563672809287580701911737271025035,
      0.00320497022805503801118415655381541759643,
      0.00308627868260877706325624133564397946129,
      0.00297606398355040882602116255686080370692,
      0.00287344936235246638755235148906672207372,
      0.00277767492975269360359490376220667282839
      ];
    stored_dn = numel(d_n);

    ## Indices for precalculated vs to-be-calculated values
    n_precalc = (n > 0) & (n < stored_dn);
    x_precalc =  (x > 0) & (x < stored_dn);
    nx_precalc = (nx > 0) & (nx < stored_dn);

    [delta_n, delta_x, delta_nx] = deal (zeros (size (x)));
    ## Fetch precalculated values
    delta_n(n_precalc) = d_n(n(n_precalc));
    delta_x(x_precalc) = d_n(x(x_precalc));
    delta_nx(nx_precalc) = d_n(nx(nx_precalc));

    ## Calculate any other d(n) values
    delta_n(!n_precalc) = delta_fn (n(!n_precalc));
    delta_x(!x_precalc) = delta_fn (x(!x_precalc));
    delta_nx(!nx_precalc) = delta_fn (nx(!nx_precalc));

    ## Calculate exp(log(pdf));
    y = exp ((delta_n - delta_x - delta_nx - ...
                       deviance (x, n .* ps) - ...
                        deviance (nx, n .* q)) - ...
                         0.5 * (log(2*pi) + log (x) + log (1-x./n)));
endfunction

function y = delta_fn (n)
  ## Stirling formula error term approximations based on Loader paper.
  ## exact expression, n^n overflows to Inf for n > ~145:
  ## = log (n!*exp(n)/(sqrt(2pi*n)*n^n));
  ##
  ## Rewritten to avoid overflow out to n> 1e305.  accurate to ~10^-12
  ## = n + gammaln (n+1) - (n+0.5) * log(n) - log(2*pi)/2;
  ##
  ## Approximated as:
  ## accurate to ~10^-16 for n=30. underflow to 0 at n~10^309
  ## = 1/(12n)-1/(360n^3)+1/(1260n^5)- 1/(1680n.^7)+1/(1188n^9) + O(n^-11);
  ##
  ## Factored to reduced operation count. Used by Loader and in R:
  ## 25% faster than unfactored form.
  ##   =(1/12-(1/360-(1/1260-(1/1680-(1/1188)/n^2)/n^2)/n^2)/n^2)/n;

   nn = n.^2;
   y = (0.08333333333333333333333333333333333333333 - ...
       (0.00277777777777777777777777777777777777778 - ...
       (0.00079365079365079365079365079365079365079 - ...
       (0.00059523809523809523809523809523809523810 - ...
       (0.00084175084175084175084175084175084175084)./nn)./nn)./nn)./nn)./n;
endfunction

function D = deviance (x, np)
  ## requires equal length column inputs
  epsilon = x ./ np;
  v = (epsilon - 1) ./ (epsilon + 1);
  vtest = abs (v) < 0.1;

  if (any (vtest))
    ## For abs(v) < 0.1, do taylor expansion for higher precision.  Expansion
    ## term: v^(2j+1)/(2j+1).  For abs(v)< 0.1, term drops slowest for max
    ## abs(v) = 0.1. (n+1)th term is <<eps (1st_term) by n=10 for doubles.
    ## jmax should be >= 10.
    jmax = 12;

    two_jpone = 2 * [1:jmax] + 1; # sum term 2*j+1 (row vector expansion)

    D = zeros (numel (epsilon), 1);

    ## D = (x-np)*v + 2*x*sum_over_j(v^2j+1 / 2j+1)
    D(vtest) = (x(vtest) - np(vtest)) .* v(vtest) + 2 .* x(vtest) .* ...
             sum (v(vtest).^(two_jpone) ./ two_jpone, 2);

    D(! vtest) = x(! vtest) .* (log (epsilon(! vtest)) - 1) + np(! vtest);
  else
    D = x.* (log (epsilon) - 1) + np;
  endif
endfunction

%!demo
%! ## Plot various PDFs from the binomial distribution
%! x = 0:40;
%! y1 = binopdf (x, 20, 0.5);
%! y2 = binopdf (x, 20, 0.7);
%! y3 = binopdf (x, 40, 0.5);
%! plot (x, y1, "*b", x, y2, "*g", x, y3, "*r")
%! grid on
%! ylim ([0, 0.25])
%! legend ({"n = 20, ps = 0.5", "n = 20, ps = 0.7", ...
%!          "n = 40, ps = 0.5"}, "location", "northeast")
%! title ("Binomial PDF")
%! xlabel ("values in x (number of successes)")
%! ylabel ("density")

## Test output
%!shared x, y
%! x = [-1 0 1 2 3];
%! y = [0 1/4 1/2 1/4 0];
%!assert (binopdf (x, 2 * ones (1, 5), 0.5 * ones (1, 5)), y, eps)
%!assert (binopdf (x, 2, 0.5 * ones (1, 5)), y, eps)
%!assert (binopdf (x, 2 * ones (1, 5), 0.5), y, eps)
%!assert (binopdf (x, 2 * [0 -1 NaN 1.1 1], 0.5), [0 NaN NaN NaN 0])
%!assert (binopdf (x, 2, 0.5 * [0 -1 NaN 3 1]), [0 NaN NaN NaN 0])
%!assert (binopdf ([x, NaN], 2, 0.5), [y, NaN], eps)
%!assert (binopdf (cat (3, x, x), 2, 0.5), cat (3, y, y), eps)

## Test Special input values
%!assert (binopdf (1, 1, 1), 1)
%!assert (binopdf (0, 3, 0), 1)
%!assert (binopdf (2, 2, 1), 1)
%!assert (binopdf (1, 2, 1), 0)
%!assert (binopdf (0, 1.1, 0), NaN)
%!assert (binopdf (1, 2, -1), NaN)
%!assert (binopdf (1, 2, 1.5), NaN)

## Test empty inputs
%!assert (binopdf ([], 1, 1), [])
%!assert (binopdf (1, [], 1), [])
%!assert (binopdf (1, 1, []), [])
%!assert (binopdf (ones (1, 0), 2, .5), ones(1, 0))
%!assert (binopdf (ones (0, 1), 2, .5), ones(0, 1))
%!assert (binopdf (ones (0, 1, 2), 2, .5), ones(0, 1, 2))
%!assert (binopdf (1, ones (0, 1, 2), .5), ones(0, 1, 2))
%!assert (binopdf (1, 2, ones (0, 1, 2)), ones(0, 1, 2))
%!assert (binopdf (ones (1, 0, 2), 2, .5), ones(1, 0, 2))
%!assert (binopdf (ones (1, 2, 0), 2, .5), ones(1, 2, 0))
%!assert (binopdf (ones (0, 1, 2), NaN, .5), ones(0, 1, 2))
%!assert (binopdf (ones (0, 1, 2), 2, NaN), ones(0, 1, 2))

## Test class of input preserved
%!assert (binopdf (single ([x, NaN]), 2, 0.5), single ([y, NaN]))
%!assert (binopdf ([x, NaN], single (2), 0.5), single ([y, NaN]))
%!assert (binopdf ([x, NaN], 2, single (0.5)), single ([y, NaN]))

## Test input validation
%!error<binopdf: function called with too few input arguments.> binopdf ()
%!error<binopdf: function called with too few input arguments.> binopdf (1)
%!error<binopdf: function called with too few input arguments.> binopdf (1, 2)
%!error<binopdf: function called with too many inputs> binopdf (1, 2, 3, 4)
%!error<binopdf: X, N, and PS must be of common size or scalars.> ...
%! binopdf (ones (3), ones (2), ones (2))
%!error<binopdf: X, N, and PS must be of common size or scalars.> ...
%! binopdf (ones (2), ones (3), ones (2))
%!error<binopdf: X, N, and PS must be of common size or scalars.> ...
%! binopdf (ones (2), ones (2), ones (3))
%!error<binopdf: X, N, and PS must not be complex.> binopdf (i, 2, 2)
%!error<binopdf: X, N, and PS must not be complex.> binopdf (2, i, 2)
%!error<binopdf: X, N, and PS must not be complex.> binopdf (2, 2, i)
