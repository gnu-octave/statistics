## Copyright (C) 2021 Nicholas R. Jankowski <jankowskin@asme.org>
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{muhat} =} expfit (@var{x})
## @deftypefnx {statistics} {[@var{muhat}, @var{muci}] =} expfit (@var{x})
## @deftypefnx {statistics} {[@var{muhat}, @var{muci}] =} expfit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@dots{}] =} expfit (@var{x}, @var{alpha}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} expfit (@var{x}, @var{alpha}, @var{censor}, @var{freq})
##
## Estimate mean and confidence intervals for the exponential distribution.
##
## @code{@var{muhat} = expfit (@var{x})} returns the maximum likelihood estimate
## of the mean parameter, @var{muhat}, of the exponential distribution given the
## data in @var{x}. @var{x} is expected to be a non-negative vector.  If @var{x}
## is an array, the mean will be computed for each column of @var{x}.  If any
## elements of @var{x} are NaN, that vector's mean will be returned as NaN.
##
## @code{[@var{muhat}, @var{muci}] = expfit (@var{x})} returns the 95%
## confidence intervals for the parameter estimate.  If @var{x} is a vector,
## @var{muci} is a two element column vector.  If @var{x} is an array, each
## column of data will have a confidence interval returned as a two-row array.
##
## @code{[@dots{}] = evfit (@var{x}, @var{alpha})} also returns the
## @qcode{100 * (1 - @var{alpha})} percent confidence intervals for the
## parameter estimates.  By default, the optional argument @var{alpha} is
## 0.05 corresponding to 95% confidence intervals.  Pass in @qcode{[]} for
## @var{alpha} to use the default values.  Any invalid values for @var{alpha}
## will return NaN for both CI bounds.
##
## @code{[@dots{}] = expfit (@var{x}, @var{alpha}, @var{censor})} accepts a
## logical or numeric array, @var{censor}, of the same size as @var{x} with
## @qcode{1}s for observations that are right-censored and @qcode{0}s for
## observations that are observed exactly.  Any non-zero elements are regarded
## as @qcode{1}s.  By default, or if left empty,
## @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @code{[@dots{}] = expfit (@var{x}, @var{alpha}, @var{censor}, @var{freq})}
## accepts a frequency array, @var{freq}, of the same size as @var{x}.
## @var{freq} typically contains integer frequencies for the corresponding
## elements in @var{x}, but it can contain any non-integer non-negative values.
## By default, or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## Matlab incompatibility: Matlab's @code{expfit} produces unpredictable results
## for some cases with higher dimensions (specifically 1 x m x n x ... arrays).
## Octave's implementation allows for @math{nxD} arrays, consistently performing
## calculations on individual column vectors.  Additionally, @var{censor} and
## @var{freq} can be used with arrays of any size, whereas Matlab only allows
## their use when @var{x} is a vector.
##
## A common alternative parameterization of the exponential distribution is to
## use the parameter @math{λ} defined as the mean number of events in an
## interval as opposed to the parameter @math{μ}, which is the mean wait time
## for an event to occur. @math{λ} and @math{μ} are reciprocals,
## i.e. @math{μ = 1 / λ}.
##
## Further information about the exponential distribution can be found at
## @url{https://en.wikipedia.org/wiki/Exponential_distribution}
##
## @seealso{expcdf, expinv, explpdf, exprnd, explike, expstat}
## @end deftypefn

function [muhat, muci] = expfit (x, alpha = 0.05, censor = [], freq = [])

  ## Check arguments
  if (nargin == 0 || nargin > 4 || nargout > 2)
    print_usage ();
  endif

  if (! (isnumeric (x) || islogical (x)))
    x = double (x);
  endif

  ## Guarantee working with column vectors
  if (isvector (x))
    x = x(:);
  endif

  if (any (x(:) < 0))
    error ("expfit: X cannot be negative.");
  endif

  sz_s = size (x);

  if (isempty (alpha))
    alpha = 0.05;
  elseif (! (isscalar (alpha)))
    error ("expfit: ALPHA must be a scalar quantity.");
  endif

  if (isempty (censor) && isempty (freq))
    ## Simple case without freq or censor, shortcut other validations
    muhat = mean (x, 1);

    if (nargout == 2)
      X = sum (x, 1);
      muci = [2*X./chi2inv(1 - alpha / 2, 2 * sz_s(1));...
              2*X./chi2inv(alpha / 2, 2 * sz_s(1))];
    endif
  else

    ## Input validation for censor and freq

    if (isempty (censor))
      ## Expand to full censor with values that don't affect results
      censor = zeros (sz_s);
    elseif (! (isnumeric (censor) || islogical (censor)))
      ## Check for incorrect freq type
      error ("expfit: CENSOR must be a numeric or logical array.")
    elseif (isvector (censor))
      ## Guarantee working with a column vector
      censor = censor(:);
    endif

    if (isempty (freq))
      ## Expand to full censor with values that don't affect results
      freq = ones (sz_s);
    elseif (! (isnumeric (freq) || islogical (freq)))
      ## Check for incorrect freq type
      error ("expfit: FREQ must be a numeric or logical array.")
    elseif (isvector (freq))
      ## Guarantee working with a column vector
      freq = freq(:);
    endif

    ## Check that size of censor and freq match x
    if (! (isequal (size (censor), sz_s)))
      error("expfit: X and CENSOR vectors mismatch.");
    elseif (! isequal (size (freq), sz_s))
      error("expfit: X and FREQ vectors mismatch.");
    endif

    ## Trivial case where censor and freq have no effect
    if (all (censor(:) == 0 & freq(:) == 1))

      muhat = mean (x, 1);

      if (nargout == 2)
        X = sum (x, 1);
        muci = [2*X./chi2inv(1 - alpha / 2, 2 * sz_s(1));...
                2*X./chi2inv(alpha / 2, 2 * sz_s(1))];
      endif

    ## No censoring, just adjust sample counts for freq
    elseif (all (censor(:) == 0))

      X = sum (x.*freq, 1);
      n = sum (freq, 1);
      muhat = X ./ n;

      if (nargout == 2)
        muci = [2*X./chi2inv(1 - alpha / 2, 2 * n);...
                2*X./chi2inv(alpha / 2, 2 * n)];
      endif

    ## Censoring, but no sample counts adjustment
    elseif (all (freq(:) == 1))

      censor = logical(censor); # convert any numeric censor'x to 0s and 1s
      X = sum (x, 1);
      r = sz_s(1) - sum (censor, 1);
      muhat = X ./ r;

      if (nargout == 2)
        muci = [2*X./chi2inv(1 - alpha / 2, 2 * r);...
                2*X./chi2inv(alpha / 2, 2 * r)];
      endif

    ## Both censoring and sample count adjustment
    else

      censor = logical (censor); # convert any numeric censor'x to 0s and 1s
      X = sum (x .* freq , 1);
      r = sum (freq .* (! censor), 1);
      muhat = X ./ r;

      if (nargout == 2)
        muci = [2*X./chi2inv(1 - alpha / 2, 2 * r);...
                2*X./chi2inv(alpha / 2, 2 * r)];
      endif
    endif

    ## compatibility check, NaN for columns where all censor's or freq's remove
    ## all samples
    null_columns = all (censor) | ! all (freq);
    muhat(null_columns) = NaN;

    if (nargout == 2)
      muci(:,null_columns) = NaN;
    endif
  endif

endfunction

%!demo
%! ## Sample 3 populations from 3 different exponential distributions
%! rande ("seed", 1);   # for reproducibility
%! r1 = exprnd (2, 4000, 1);
%! rande ("seed", 2);   # for reproducibility
%! r2 = exprnd (5, 4000, 1);
%! rande ("seed", 3);   # for reproducibility
%! r3 = exprnd (12, 4000, 1);
%! r = [r1, r2, r3];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, 48, 0.52);
%! h = findobj (gca, "Type", "patch");
%! set (h(1), "facecolor", "c");
%! set (h(2), "facecolor", "g");
%! set (h(3), "facecolor", "r");
%! hold on
%!
%! ## Estimate their mu parameter
%! muhat = expfit (r);
%!
%! ## Plot their estimated PDFs
%! x = [0:max(r(:))];
%! y = exppdf (x, muhat(1));
%! plot (x, y, "-pr");
%! y = exppdf (x, muhat(2));
%! plot (x, y, "-sg");
%! y = exppdf (x, muhat(3));
%! plot (x, y, "-^c");
%! ylim ([0, 0.6])
%! xlim ([0, 40])
%! legend ({"Normalized HIST of sample 1 with μ=2", ...
%!          "Normalized HIST of sample 2 with μ=5", ...
%!          "Normalized HIST of sample 3 with μ=12", ...
%!          sprintf("PDF for sample 1 with estimated μ=%0.2f", muhat(1)), ...
%!          sprintf("PDF for sample 2 with estimated μ=%0.2f", muhat(2)), ...
%!          sprintf("PDF for sample 3 with estimated μ=%0.2f", muhat(3))})
%! title ("Three population samples from different exponential distributions")
%! hold off

## Tests for mean
%!assert (expfit (1), 1)
%!assert (expfit (1:3), 2)
%!assert (expfit ([1:3]'), 2)
%!assert (expfit (1:3, []), 2)
%!assert (expfit (1:3, [], [], []), 2)
%!assert (expfit (magic (3)), [5 5 5])
%!assert (expfit (cat (3, magic (3), 2*magic (3))), cat (3,[5 5 5], [10 10 10]))
%!assert (expfit (1:3, 0.1, [0 0 0], [1 1 1]), 2)
%!assert (expfit ([1:3]', 0.1, [0 0 0]', [1 1 1]'), 2)
%!assert (expfit (1:3, 0.1, [0 0 0]', [1 1 1]'), 2)
%!assert (expfit (1:3, 0.1, [1 0 0], [1 1 1]), 3)
%!assert (expfit (1:3, 0.1, [0 0 0], [4 1 1]), 1.5)
%!assert (expfit (1:3, 0.1, [1 0 0], [4 1 1]), 4.5)
%!assert (expfit (1:3, 0.1, [1 0 1], [4 1 1]), 9)
%!assert (expfit (1:3, 0.1, [], [-1 1 1]), 4)
%!assert (expfit (1:3, 0.1, [], [0.5 1 1]), 2.2)
%!assert (expfit (1:3, 0.1, [1 1 1]), NaN)
%!assert (expfit (1:3, 0.1, [], [0 0 0]), NaN)
%!assert (expfit (reshape (1:9, [3 3])), [2 5 8])
%!assert (expfit (reshape (1:9, [3 3]), [], eye(3)), [3 7.5 12])
%!assert (expfit (reshape (1:9, [3 3]), [], 2*eye(3)), [3 7.5 12])
%!assert (expfit (reshape (1:9, [3 3]), [], [], [2 2 2; 1 1 1; 1 1 1]), ...
%! [1.75 4.75 7.75])
%!assert (expfit (reshape (1:9, [3 3]), [], [], [2 2 2; 1 1 1; 1 1 1]), ...
%! [1.75 4.75 7.75])
%!assert (expfit (reshape (1:9, [3 3]), [], eye(3), [2 2 2; 1 1 1; 1 1 1]), ...
%! [3.5 19/3 31/3])

## Tests for confidence intervals
%!assert ([~,muci] = expfit (1:3, 0), [0; Inf])
%!assert ([~,muci] = expfit (1:3, 2), [Inf; 0])
%!assert ([~,muci] = expfit (1:3, 0.1, [1 1 1]), [NaN; NaN])
%!assert ([~,muci] = expfit (1:3, 0.1, [], [0 0 0]), [NaN; NaN])
%!assert ([~,muci] = expfit (1:3, -1), [NaN; NaN])
%!assert ([~,muci] = expfit (1:3, 5), [NaN; NaN])
#!assert ([~,muci] = expfit ([1:3;1:3], -1), NaN(2, 3)]
#!assert ([~,muci] = expfit ([1:3;1:3], 5), NaN(2, 3)]
%!assert ([~,muci] = expfit (1:3), [0.830485728373393; 9.698190330474096], ...
%!             1000*eps)
%!assert ([~,muci] = expfit (1:3, 0.1), ...
%!                          [0.953017262058213; 7.337731146400207], 1000*eps)
%!assert ([~,muci] = expfit ([1:3;2:4]), ...
%!             [0.538440777613095, 0.897401296021825, 1.256361814430554; ...
%!             12.385982973214016, 20.643304955356694, 28.900626937499371], ...
%!             1000*eps)
%!assert ([~,muci] = expfit ([1:3;2:4], [], [1 1 1; 0 0 0]), ...
%!             100*[0.008132550920455, 0.013554251534091, 0.018975952147727; ...
%!             1.184936706156216, 1.974894510260360, 2.764852314364504], ...
%!             1000*eps)
%!assert ([~,muci] = expfit ([1:3;2:4], [], [], [3 3 3; 1 1 1]), ...
%!             [0.570302756652583, 1.026544961974649, 1.482787167296715; ...
%!             4.587722594914109, 8.257900670845396, 11.928078746776684], ...
%!             1000*eps)
%!assert ([~,muci] = expfit ([1:3;2:4], [], [0 0 0; 1 1 1], [3 3 3; 1 1 1]), ...
%!             [0.692071440311161, 1.245728592560089, 1.799385744809018; ...
%!             8.081825275395081, 14.547285495711145, 21.012745716027212], ...
%!             1000*eps)

%!test
%! x = reshape (1:8, [4 2]);
%! x(4) = NaN;
%! [muhat,muci] = expfit (x);
%! assert ({muhat, muci}, {[NaN, 6.5], ...
%!         [NaN, 2.965574334593430;NaN, 23.856157493553368]}, 1000*eps);

%!test
%! x = magic (3);
%! censor = [0 1 0; 0 1 0; 0 1 0];
%! freq = [1 1 0; 1 1 0; 1 1 0];
%! [muhat,muci] = expfit (x, [], censor, freq);
%! assert ({muhat, muci}, {[5 NaN NaN], ...
%!                 [[2.076214320933482; 24.245475826185242],NaN(2)]}, 1000*eps);

## Test input validation
%!error expfit ()
%!error expfit (1,2,3,4,5)
%!error [a b censor] = expfit (1)
%!error <ALPHA must be a scalar quantity> expfit (1, [1 2])
%!error <X cannot be negative> expfit ([-1 2 3 4 5])
%!error <CENSOR must be a numeric or logical array> expfit ([1:5], [], "test")
%!error <FREQ must be a numeric or logical array> expfit ([1:5], [], [], "test")
%!error <X and CENSOR vectors mismatch.> expfit ([1:5], [], [0 0 0 0])
%!error <X and FREQ vectors mismatch.> expfit ([1:5], [], [], [1 1 1 1])
