## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {[@var{r}, @var{tieadj}]} = tiedrank (@var{x})
## @deftypefnx {statistics} {[@var{r}, @var{tieadj}]} = tiedrank (@var{x}, @var{tieflag})
## @deftypefnx {statistics} {[@var{r}, @var{tieadj}]} = tiedrank (@var{x}, @var{tieflag}, @var{bidir})
##
## @code{[@var{r}, @var{tieadj}] = tiedrank (@var{x})} computes the ranks of the
## values in vector @var{x}.  If any values in @var{x} are tied, @code{tiedrank}
## computes their average rank.  The return value @var{tieadj} is an adjustment
## for ties required by the nonparametric tests @code{signrank} and
## @code{ranksum}, and for the computation of Spearman's rank correlation.
##
## @code{[@var{r}, @var{tieadj}] = tiedrank (@var{x}, 1)} computes the ranks of
## the values in the vector @var{x}. @var{tieadj} is a vector of three
## adjustments for ties required in the computation of Kendall's tau.
## @code{tiedrank (@var{x}, 0)} is the same as @code{tiedrank (@var{x})}.
##
## @code{[@var{r}, @var{tieadj}] = tiedrank (@var{x}, 0, 1)} computes the ranks
## from each end, so that the smallest and largest values get rank 1, the next
## smallest and largest get rank 2, etc.  These ranks are used in the
## Ansari-Bradley test.
##
## @end deftypefn

function [r, tieadj] = tiedrank (x, tieflag, bidir)
  ## Check input arguments and add defauls
  if (nargin < 1 || nargin > 3)
    print_usage ();
  endif
  if (! isvector (x))
    error ("tiedrank: X must be a vector.");
  endif
  if (nargin < 2)
    tieflag = false;
  elseif (! isscalar (tieflag) || ! (isnumeric (tieflag) || isbool (tieflag)))
    error ("tiedrank: TIEFLAG must be a numeric or boolean scalar.");
  endif
  if (nargin < 3)
    bidir = false;
  elseif (! isscalar (bidir) || ! (isnumeric (bidir) || isbool (bidir)))
    error ("tiedrank: BIDIR must be a numeric or boolean scalar.");
  endif

  ## Sort X and leave NaNs at the end of vector
  [sx, idx] = sort (x(:));
  NaNs = sum (isnan (x));
  xLen = length (x) - NaNs;

  ## Count ranks from low end
  if (! bidir)
    ranks = [1:xLen, NaN(1,NaNs)]';
  ## Count ranks from both ends
  else
    ## For even number of samples
    if (mod (xLen, 2) == 0)
      ranks = [(1:xLen/2), (xLen/2:-1:1), NaN(1,NaNs)]';
    ## For odd number of samples
    else
      ranks = [(1:(xLen+1)/2), ((xLen-1)/2:-1:1), NaN(1,NaNs)]';
    endif
  endif

  ## Define number of adjustments
  if (! tieflag)
    tieadj = 0;
  else
    tieadj = [0; 0; 0];
  endif

  ## Check precision of X
  if (isa (x, "single"))
    ranks = single (ranks);
    tieadj = single (tieadj);
  endif

  ## Adjust for ties
  ties = sx(1:xLen-1) >= sx(2:xLen);
  tieloc = [find(ties); xLen+2];
  maxTies = length (tieloc);
  tiecount = 1;
  while (tiecount < maxTies)
    tiestart = tieloc(tiecount);
    ntied = 2;
    while (tieloc(tiecount + 1) == tieloc(tiecount) + 1)
      tiecount = tiecount + 1;
      ntied = ntied + 1;
    endwhile
    if (! tieflag)
      tieadj = tieadj + ntied * (ntied - 1) * (ntied + 1) / 2;
    else
      n2minusn = ntied * (ntied - 1);
      tieadj = tieadj + [n2minusn/2; n2minusn*(ntied-2); n2minusn*(2*ntied+5)];
    endif
    ## Compute mean of tied ranks
    ranks(tiestart:tiestart+ntied-1) = ...
                  sum (ranks(tiestart:tiestart+ntied-1)) / ntied;
    tiecount = tiecount + 1;
  endwhile

  ## Reshape ranks including NaN where required.
  r(idx) = ranks;
  r = reshape (r, size (x));
endfunction

## testing against mileage data and results from Matlab
%!test
%! [r,tieadj] = tiedrank ([10, 20, 30, 40, 20]);
%! assert (r, [1, 2.5, 4, 5, 2.5]);
%! assert (tieadj, 3);
%!test
%! [r,tieadj] = tiedrank ([10; 20; 30; 40; 20]);
%! assert (r, [1; 2.5; 4; 5; 2.5]);
%! assert (tieadj, 3);
%!test
%! [r,tieadj] = tiedrank ([10, 20, 30, 40, 20], 1);
%! assert (r, [1, 2.5, 4, 5, 2.5]);
%! assert (tieadj, [1; 0; 18]);
%!test
%! [r,tieadj] = tiedrank ([10, 20, 30, 40, 20], 0, 1);
%! assert (r, [1, 2.5, 2, 1, 2.5]);
%! assert (tieadj, 3);
%!test
%! [r,tieadj] = tiedrank ([10, 20, 30, 40, 20], 1, 1);
%! assert (r, [1, 2.5, 2, 1, 2.5]);
%! assert (tieadj, [1; 0; 18]);

## Test input validation
%!error <tiedrank: X must be a vector.> tiedrank (ones (2))
%!error <tiedrank: TIEFLAG must be a numeric or boolean scalar.> ...
%! tiedrank ([1, 2, 3, 4, 5], [1, 1])
%!error <tiedrank: TIEFLAG must be a numeric or boolean scalar.> ...
%! tiedrank ([1, 2, 3, 4, 5], "A")
%!error <tiedrank: TIEFLAG must be a numeric or boolean scalar.> ...
%! tiedrank ([1, 2, 3, 4, 5], [true, true])
%!error <tiedrank: BIDIR must be a numeric or boolean scalar.> ...
%! tiedrank ([1, 2, 3, 4, 5], 0, [1, 1])
%!error <tiedrank: BIDIR must be a numeric or boolean scalar.> ...
%! tiedrank ([1, 2, 3, 4, 5], 0, "A")
%!error <tiedrank: BIDIR must be a numeric or boolean scalar.> ...
%! tiedrank ([1, 2, 3, 4, 5], 0, [true, true])
