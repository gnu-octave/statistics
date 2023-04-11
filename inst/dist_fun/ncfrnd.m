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
## @deftypefn  {statistics} {@var{r} =} ncfrnd (@var{df1}, @var{df2}, @var{delta})
## @deftypefnx {statistics} {@var{r} =} ncfrnd (@var{df1}, @var{df2}, @var{delta}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} ncfrnd (@var{df1}, @var{df2}, @var{delta}, [@var{sz}])
##
## Random arrays from the noncentral F distribution.
##
## @code{@var{x} = ncfrnd (@var{p}, @var{df1}, @var{df2}, @var{delta})} returns
## an array of random numbers chosen from the noncentral F distribution with
## parameters  @var{df1}, @var{df2}, @var{delta}).  The size of @var{r} is the
## common size of @var{df1}, @var{df2}, and @var{delta}.  A scalar input
## functions as a constant matrix of the same size as the other input.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## @code{ncfrnd} generates values using the definition of a noncentral F random
## variable, as the ratio of a noncentral chi-square and a (central) chi-square.
##
## @seealso{ncfcdf, ncfinv, ncfpdf, ncfstat}
## @end deftypefn

function r = ncfrnd (df1, df2, delta, varargin)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("ncfrnd: too few input arguments.");
  endif

  ## Check for appropriate class
  if (isa (df1, "single") || isa (df2, "single") || isa (delta, "single"));
    is_class = "single";
  else
    is_class = "double";
  endif

  ## Check for additional dimensions in varargin and get their size
  dim_vec = 1;
  if (nargin > 2)
    extra_varargin = numel (varargin(:));
    if (extra_varargin == 1)
      size_dim = varargin{1};
      ## Check for empty input argument
      if (isempty (size_dim))
        error (strcat (["ncfrnd: extra argument for size of output"], ...
                       [" array cannot be empty."]));
      endif
      dim_vec = zeros (size_dim, is_class);
    elseif (extra_varargin > 1)
      for i = 1:extra_varargin
        size_dim(i) = varargin{i};
      endfor
      dim_vec = zeros (size_dim, is_class);
    endif
  endif

  ## Check for common size of MU, SIGMA, and output based on given dimensions
  if (! isscalar (df1) || ! isscalar (df2) || ...
      ! isscalar (delta) || ! isscalar (dim_vec))
    [err, df1, df2, delta, dim_vec] = common_size (df1, df2, delta, dim_vec);
    if (err > 0)
      error (strcat (["ncfrnd: DF1, DF2, DELTA, and DIM vector must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  ## Get final dimensions of returning random array
  size_out = size (df1);

  ## Return NaNs for out of range values of DF and DELTA
  df1(df1 <= 0) = NaN;
  df2(df2 <= 0) = NaN;
  delta(delta <= 0) = NaN;

  r = (ncx2rnd (df1, delta, size_out) ./ df1) ./ ...
      (2 .* randg (df2 ./ 2, size_out) ./ df2);
endfunction

## Test input validation
%!error<ncfrnd: too few input arguments.> ncfrnd ()
%!error<ncfrnd: too few input arguments.> ncfrnd (1)
%!error<ncfrnd: too few input arguments.> ncfrnd (1, 2)
%!error<ncfrnd: extra argument for size of output array cannot be empty.> ...
%! ncfrnd (5, ones (3), ones (2), [])
%!error<ncfrnd: DF1, DF2, DELTA, and DIM vector must be of common size or scalars.> ...
%! ncfrnd (5, ones (3), ones (2))
%!error<ncfrnd: DF1, DF2, DELTA, and DIM vector must be of common size or scalars.> ...
%! ncfrnd (5, ones (2), ones (2), 3, 2)
%!error<ncfrnd: DF1, DF2, DELTA, and DIM vector must be of common size or scalars.> ...
%! ncfrnd (5, ones (2), ones (2), 1, 2)

## Output validation tests
%!assert (size (ncfrnd (5, 2, 3, 3, 5, 7)), [3, 5, 7])
%!assert (size (ncfrnd (5, 2, 3, [3, 5, 7])), [3, 5, 7])
%!assert (size (ncfrnd (5, ones (3, 5), 2 * ones (3, 5), [3, 5])), [3, 5])
%!assert (size (ncfrnd (2, 3, 5)), [1, 1])
