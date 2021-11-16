## Copyright (C) 2021 Nicholas R. Jankowski <jankowskin@asme.org>
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
## @deftypefn {Function File} {@var{mu} =} expfit (@var{s})
## @deftypefnx {Function File} {[@var{mu}, @var{ci}] =} expfit (@var{s})
## @deftypefnx {Function File} {[@var{mu}, @var{ci}] =} expfit (@var{s}, @var{alpha})
## @deftypefnx {Function File} {@dots{} =} expfit (@var{s}, @var{alpha}, @var{c})
## @deftypefnx {Function File} {@dots{} =} expfit (@var{s}, @var{alpha}, @var{c}, @var{f})
##
## Estimate the mean of the exponential probability distribution function from
## which sample data @var{s} has been taken.  @var{s} is expected to be a
## non-negative vector.  If @var{s} is an array, the mean will be computed for
## each column of @var{s}.  If any elements of @var{s} are NaN, that vector's
## mean will be returned as NaN.
##
## If the optional output variable @var{ci} is requested, @code{expfit} will
## also return the confidence interval bounds for the estimate as a two element
## column vector.  If @var{s} is an array, each column of data will have a
## confidence interval returned as a two row array.
##
## The optional scalar input @var{alpha} can be used to define the
## (1-@var{alpha}) confidence interval to be applied to all estimates as a
## value between 0 and 1.  The default is 0.05, resulting in a 0.95 or 95% CI.
## Any invalid values for alpha will return NaN for both CI bounds.
##
## The optional input @var{c} is a logical or numeric array of zeros and ones
## the same size as @var{s}, used to right-censor individual elements of
## @var{s}.  A value of 1 indicates the data should be censored from
## the mean estimation.  Any nonzero values in @var{c} are treated as a 1.
##
## The optional input @var{f} is a numeric array the same size as @var{s}, used
## to specify occurrence frequencies for the elements in @var{s}.  Values of
## @var{f} need not be integers.  Any NaN elements in the frequency array will
## produce a NaN output for @var{mu}.
##
## Options can be skipped by using [] to revert to the default.
##
## Matlab incompatibility: Matlab's @code{expfit} produces unpredictable results
## for some cases with higher dimensions (specifically 1 x m x n x ... arrays).
## Octave's implementation allows for n-D arrays, consistently performing
## calculations on individual column vectors.  Additionally, @var{c} and @var{f}
## can be used with arrays of any size, whereas Matlab only allows their use
## when @var{s} is a vector.
##
## @end deftypefn
##
## @seealso{expcdf, expinv, explike, exppdf, exprnd, expstat}
## @seealso{expstat, exprnd, expcdf, expinv}

function [m, v] = expfit (s, alpha = 0.05, c = [], f = [])

  ## Check arguments
  if (nargin ==0 || nargin > 4 || nargout > 2)
    print_usage ();
  endif

  if ! (isnumeric (s) || islogical (s))
    s = double(s);
  endif

  ## guarantee working with column vectors
  if isvector (s)

    s = s(:);
  endif

  if any (s(:) < 0)
    error("expfit: input data S cannot be negative");
  endif

  sz_s = size (s);

  if (isempty (alpha))
    alpha = 0.05;
  elseif !(isscalar (alpha))
    error ("expfit: ALPHA must be a scalar quantity");
  endif

  if (isempty (c) && isempty (f))
    ##simple case without f or c, shortcut other validations
    m = mean (s, 1);

    if (nargout == 2)
      S = sum (s, 1);
      v = [2*S ./ chi2inv(1 - alpha/2, 2*sz_s(1));...
           2*S ./ chi2inv(alpha/2, 2*sz_s(1))];
    endif
  else

    ## input validation for c and f

    if (isempty (c))
      ##expand to full c with values that don't affect results
      c = zeros (sz_s);
    elseif (! (isnumeric(c) || islogical (c)))
      #check for incorrect f type
      error ("expfit: C must be a numeric or logical array")
    elseif (isvector (c))
      ## guarantee working with a column vector
      c = c(:);
    endif

    if (isempty (f))
      ##expand to full c with values that don't affect results
      f = ones (sz_s);
    elseif (! (isnumeric(f) || islogical (f)))
      #check for incorrect f type
      error ("expfit: F must be a numeric or logical array")
    elseif (isvector (f))
      ## guarantee working with a column vector
      f = f(:);
    endif

    #check that size of c and f match s
    if !(isequal(size (c), sz_s))
      error("expfit: C must be the same size as S");
    elseif (! isequal(size (f), sz_s))
      error("expfit: F must be the same size as S");
    endif

    ## trivial case where c and f have no effect
    if (all (c(:) == 0 & f(:) == 1))

      m = mean (s, 1);

      if (nargout == 2)
        S = sum (s, 1);
        v = [2*S ./ chi2inv(1 - alpha/2, 2*sz_s(1));...
             2*S ./ chi2inv(alpha/2, 2*sz_s(1))];
      endif

    ## no censoring, just adjust sample counts for f
    elseif (all (c(:) == 0))

      S = sum (s.*f, 1);
      n = sum (f, 1);
      m = S ./ n;

      if (nargout == 2)
        v = [2*S ./ chi2inv(1 - alpha/2, 2*n);...
             2*S ./ chi2inv(alpha/2, 2*n)];
      endif

    ## censoring, but no sample counts adjustment
    elseif (all (f(:) == 1))

      c = logical(c); ##convert any numeric c's to 0s and 1s
      S = sum (s, 1);
      r = sz_s(1) - sum (c, 1);
      m = S ./ r;

      if (nargout == 2)
        v = [2*S ./ chi2inv(1 - alpha/2, 2*r);...
             2*S ./ chi2inv(alpha/2, 2*r)];
      endif

    ## both censoring and sample count adjustment
    else

      c = logical(c); ##convert any numeric c's to 0s and 1s
      S = sum (s.*f , 1);
      r = sum (f.*(!c), 1);
      m = S ./ r;

      if (nargout == 2)
        v = [2*S ./ chi2inv(1 - alpha/2, 2*r);...
             2*S ./ chi2inv(alpha/2, 2*r)];
      endif
    endif

    ## compatibility check, NaN for columns where all c's or f's remove all samples
    null_columns = all (c) | ! all (f);
    m(null_columns) = NaN;

    if (nargout == 2)
      v(:,null_columns) = NaN;
    endif
  endif

endfunction

##tests for mean
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
%!assert (expfit (reshape (1:9, [3 3]), [], [], [2 2 2; 1 1 1; 1 1 1]), [1.75 4.75 7.75])
%!assert (expfit (reshape (1:9, [3 3]), [], [], [2 2 2; 1 1 1; 1 1 1]), [1.75 4.75 7.75])
%!assert (expfit (reshape (1:9, [3 3]), [], eye(3), [2 2 2; 1 1 1; 1 1 1]), [3.5 19/3 31/3])

##tests for confidence intervals
%!assert ([~,v] = expfit (1:3, 0), [0; Inf])
%!assert ([~,v] = expfit (1:3, 2), [Inf; 0])
%!assert ([~,v] = expfit (1:3, 0.1, [1 1 1]), [NaN; NaN])
%!assert ([~,v] = expfit (1:3, 0.1, [], [0 0 0]), [NaN; NaN])
%!assert ([~,v] = expfit (1:3, -1), [NaN; NaN])
%!assert ([~,v] = expfit (1:3, 5), [NaN; NaN])
#!assert ([~,v] = expfit ([1:3;1:3], -1), NaN(2, 3)]
#!assert ([~,v] = expfit ([1:3;1:3], 5), NaN(2, 3)]
%!assert ([~,v] = expfit (1:3), [0.830485728373393; 9.698190330474096], 1000*eps)
%!assert ([~,v] = expfit (1:3, 0.1), [0.953017262058213; 7.337731146400207], 1000*eps)
%!assert ([~,v] = expfit ([1:3;2:4]), ...
%!             [0.538440777613095, 0.897401296021825, 1.256361814430554; ...
%!             12.385982973214016, 20.643304955356694, 28.900626937499371], 1000*eps)
%!assert ([~,v] = expfit ([1:3;2:4], [], [1 1 1; 0 0 0]), ...
%!             100*[0.008132550920455, 0.013554251534091, 0.018975952147727; ...
%!             1.184936706156216, 1.974894510260360, 2.764852314364504], 1000*eps)
%!assert ([~,v] = expfit ([1:3;2:4], [], [], [3 3 3; 1 1 1]), ...
%!             [0.570302756652583, 1.026544961974649, 1.482787167296715; ...
%!             4.587722594914109, 8.257900670845396, 11.928078746776684], 1000*eps)
%!assert ([~,v] = expfit ([1:3;2:4], [], [0 0 0; 1 1 1], [3 3 3; 1 1 1]), ...
%!             [0.692071440311161, 1.245728592560089, 1.799385744809018; ...
%!             8.081825275395081, 14.547285495711145, 21.012745716027212], 1000*eps)

%!test
%! s = reshape (1:8, [4 2]);
%! s(4) = NaN;
%! [m,v] = expfit (s);
%! assert ({m, v}, {[NaN, 6.5], [NaN, 2.965574334593430;NaN, 23.856157493553368]}, 1000*eps);

%!test
%! s = magic (3);
%! c = [0 1 0; 0 1 0; 0 1 0];
%! f = [1 1 0; 1 1 0; 1 1 0];
%! [m,v] = expfit (s, [], c, f);
%! assert ({m, v}, {[5 NaN NaN], [[2.076214320933482; 24.245475826185242],NaN(2)]}, 1000*eps);

## input validation
%!error expfit ()
%!error expfit (1,2,3,4,5)
%!error [a b c] = expfit (1)
%!error <ALPHA must be a scalar quantity> expfit (1, [1 2])
%!error <input data S cannot be negative> expfit ([-1 2 3 4 5])
%!error <C must be a numeric or logical array> expfit ([1:5], [], "test")
%!error <F must be a numeric or logical array> expfit ([1:5], [], [], "test")
%!error <C must be the same size as S> expfit ([1:5], [], [0 0 0 0])
%!error <F must be the same size as S> expfit ([1:5], [], [], [1 1 1 1])
