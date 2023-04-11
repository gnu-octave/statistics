## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
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
## @deftypefn  {statistics} {@var{r} =} trnd (@var{df})
## @deftypefnx {statistics} {@var{r} =} trnd (@var{df}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} trnd (@var{df}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} trnd (@var{df}, [@var{sz}])
##
## Random arrays from the Student's T distribution.
##
## Return a matrix of random samples from the t (Student) distribution with
## @var{df} degrees of freedom.
##
## @code{@var{r} = trnd (@var{df})} returns an array of random numbers chosen
## from the Student's T distribution with @var{df} degrees of freedom.  The size
## of @var{r} is the size of @var{df}.  @var{df} must be a finite real number
## greater than 0, otherwise NaN is returned.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## @seealso{tcdf, tpdf, tpdf, tstat}
## @end deftypefn

function r = trnd (df, varargin)

  if (nargin < 1)
    print_usage ();
  endif

  if (nargin == 1)
    sz = size (df);
  elseif (nargin == 2)
    if (isscalar (varargin{1}) && varargin{1} >= 0)
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0))
      sz = varargin{1};
    else
      error (strcat (["trnd: dimension vector must be row vector of"], ...
                     [" non-negative integers."]));
    endif
  elseif (nargin > 2)
    if (any (cellfun (@(x) (! isscalar (x) || x < 0), varargin)))
      error ("trnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  if (! isscalar (df) && ! isequal (size (df), sz))
    error ("trnd: DF must be scalar or of size SZ.");
  endif

  if (iscomplex (df))
    error ("trnd: DF must not be complex.");
  endif

  if (isa (df, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  if (isscalar (df))
    if ((df > 0) && (df < Inf))
      r = randn (sz, cls) ./ sqrt (2*randg (df/2, sz, cls) / df);
    else
      r = NaN (sz, cls);
    endif
  else
    r = NaN (sz, cls);

    k = (df > 0) & (df < Inf);
    r(k) = randn (sum (k(:)), 1, cls) ...
             ./ sqrt (2*randg (df(k)/2, cls) ./ df(k))(:);
  endif

endfunction


%!assert (size (trnd (2)), [1, 1])
%!assert (size (trnd (ones (2,1))), [2, 1])
%!assert (size (trnd (ones (2,2))), [2, 2])
%!assert (size (trnd (1, 3)), [3, 3])
%!assert (size (trnd (1, [4 1])), [4, 1])
%!assert (size (trnd (1, 4, 1)), [4, 1])

## Test class of input preserved
%!assert (class (trnd (1)), "double")
%!assert (class (trnd (single (1))), "single")
%!assert (class (trnd (single ([1 1]))), "single")

## Test input validation
%!error trnd ()
%!error trnd (1, -1)
%!error trnd (1, ones (2))
%!error trnd (i)
%!error trnd (1, [2 -1 2])
%!error trnd (1, 2, ones (2))
%!error trnd (1, 2, -1)
%!error trnd (ones (2,2), 3)
%!error trnd (ones (2,2), [3, 2])
%!error trnd (ones (2,2), 2, 3)
