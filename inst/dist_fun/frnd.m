## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{r} =} frnd (@var{df1}, @var{df2})
## @deftypefnx {statistics} {@var{r} =} frnd (@var{df1}, @var{df2}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} frnd (@var{df1}, @var{df2}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} frnd (@var{df1}, @var{df2}, [@var{sz}])
##
## Random arrays from the F distribution.
##
## @code{@var{r} = frnd (@var{df1}, @var{df2})} returns an array of random
## numbers chosen from the F distribution with @var{df1} and @var{df2} degrees
## of freedom.  The size of @var{r} is the common size of @var{df1} and
## @var{df2}.  A scalar input functions as a constant matrix of the same size as
## the other inputs.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## @seealso{fcdf, finv, fpdf, fstat}
## @end deftypefn

function r = frnd (df1, df2, varargin)

  if (nargin < 2)
    print_usage ();
  endif

  if (! isscalar (df1) || ! isscalar (df2))
    [retval, df1, df2] = common_size (df1, df2);
    if (retval > 0)
      error ("frnd: DF1 and DF2 must be of common size or scalars.");
    endif
  endif

  if (iscomplex (df1) || iscomplex (df2))
    error ("frnd: DF1 and DF2 must not be complex.");
  endif

  if (nargin == 2)
    sz = size (df1);
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} >= 0)
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0))
      sz = varargin{1};
    else
      error (strcat (["frnd: dimension vector must be a row vector of"], ...
                     [" non-negative integers."]));
    endif
  elseif (nargin > 3)
    if (any (cellfun (@(x) (! isscalar (x) || x < 0), varargin)))
      error ("frnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  if (! isscalar (df1) && ! isequal (size (df1), sz))
    error ("frnd: DF1 and DF2 must be scalar or of size SZ.");
  endif

  if (isa (df1, "single") || isa (df2, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  if (isscalar (df1) && isscalar (df2))
    if ((df1 > 0) && (df1 < Inf) && (df2 > 0) && (df2 < Inf))
      r = df2/df1 * randg (df1/2, sz, cls) ./ randg (df2/2, sz, cls);
    else
      r = NaN (sz, cls);
    endif
  else
    r = NaN (sz, cls);

    k = (df1 > 0) & (df1 < Inf) & (df2 > 0) & (df2 < Inf);
    r(k) = df2(k) ./ df1(k) .* randg (df1(k)/2, cls) ./ randg (df2(k)/2, cls);
  endif

endfunction


%!assert (size (frnd (1,2)), [1, 1])
%!assert (size (frnd (ones (2,1), 2)), [2, 1])
%!assert (size (frnd (ones (2,2), 2)), [2, 2])
%!assert (size (frnd (1, 2*ones (2,1))), [2, 1])
%!assert (size (frnd (1, 2*ones (2,2))), [2, 2])
%!assert (size (frnd (1, 2, 3)), [3, 3])
%!assert (size (frnd (1, 2, [4 1])), [4, 1])
%!assert (size (frnd (1, 2, 4, 1)), [4, 1])

## Test class of input preserved
%!assert (class (frnd (1, 2)), "double")
%!assert (class (frnd (single (1), 2)), "single")
%!assert (class (frnd (single ([1 1]), 2)), "single")
%!assert (class (frnd (1, single (2))), "single")
%!assert (class (frnd (1, single ([2 2]))), "single")

## Test input validation
%!error frnd ()
%!error frnd (1)
%!error frnd (ones (3), ones (2))
%!error frnd (ones (2), ones (3))
%!error frnd (i, 2)
%!error frnd (2, i)
%!error frnd (1,2, -1)
%!error frnd (1,2, ones (2))
%!error frnd (1, 2, [2 -1 2])
%!error frnd (1,2, 1, ones (2))
%!error frnd (1,2, 1, -1)
%!error frnd (ones (2,2), 2, 3)
%!error frnd (ones (2,2), 2, [3, 2])
%!error frnd (ones (2,2), 2, 2, 3)
