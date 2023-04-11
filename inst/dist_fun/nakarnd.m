## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 1995-2015 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{r} =} nakarnd (@var{m}, @var{w})
## @deftypefnx {statistics} {@var{r} =} nakarnd (@var{m}, @var{w}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} nakarnd (@var{m}, @var{w}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} nakarnd (@var{m}, @var{w}, [@var{sz}])
##
## Random arrays from the Nakagami distribution.
##
## @code{@var{r} = nakarnd (@var{m}, @var{w})} returns an array of random
## numbers chosen from the Nakagami distribution with parameters @var{m} and
## @var{w}.  The size of @var{r} is the common size of @var{m} and @var{w}.  A
## scalar input functions as a constant matrix of the same size as the other
## inputs.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## @seealso{nakacdf, nakainv, nakapdf}
## @end deftypefn

function r = nakarnd (m, w, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    print_usage ();
  endif

  ## Check for common size of M and W
  if (! isscalar (m) || ! isscalar (w))
    [retval, m, w] = common_size (m, w);
    if (retval > 0)
      error ("nakarnd: M and W must be of common size or scalars.");
    endif
  endif

  ## Check for M and W being reals
  if (iscomplex (m) || iscomplex (w))
    error ("nakarnd: M and W must not be complex.");
  endif

  ## Check for SIZE vector or DIMENSION input arguments
  if (nargin == 2)
    sz = size (m);
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} >= 0)
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0))
      sz = varargin{1};
    else
      error (strcat (["nakarnd: dimension vector must be row vector"], ...
                     [" of non-negative integers."]));
    endif
  elseif (nargin > 3)
    if (any (cellfun (@(x) (! isscalar (x) || x < 0), varargin)))
      error ("nakarnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (m) && ! isequal (size (w), sz))
    error ("nakagrnd: M and W must be scalar or of size SZ.");
  endif

  ## Check for appropriate class
  if (isa (m, "single") || isa (w, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Generate random sample from Nakagami distribution
  if (isscalar (m) && isscalar (w))
    if ((0 < m) && (m < Inf) && (0 < w) && (w < Inf))
      m_gamma = m;
      w_gamma = w / m;
      r = gamrnd (m_gamma, w_gamma, sz);
      r = sqrt (r);
    else
      r = NaN (sz, cls);
    endif
  else
    r = NaN (sz, cls);
    k = (0 < m) & (m < Inf) & (0 < w) & (w < Inf);
    m_gamma = m;
    w_gamma = w ./ m;
    r(k) = gamrnd (m_gamma(k), w_gamma(k));
    r(k) = sqrt (r(k));
  endif

endfunction


%!assert (size (nakarnd (1,1)), [1, 1])
%!assert (size (nakarnd (ones (2,1), 1)), [2, 1])
%!assert (size (nakarnd (ones (2,2), 1)), [2, 2])
%!assert (size (nakarnd (1, ones (2,1))), [2, 1])
%!assert (size (nakarnd (1, ones (2,2))), [2, 2])
%!assert (size (nakarnd (1,1, 3)), [3, 3])
%!assert (size (nakarnd (1,1, [4 1])), [4, 1])
%!assert (size (nakarnd (1,1, 4, 1)), [4, 1])

## Test class of input preserved
%!assert (class (nakarnd (1,1)), "double")
%!assert (class (nakarnd (single (1),1)), "single")
%!assert (class (nakarnd (single ([1 1]),1)), "single")
%!assert (class (nakarnd (1,single (1))), "single")
%!assert (class (nakarnd (1,single ([1 1]))), "single")

## Test input validation
%!error nakarnd ()
%!error nakarnd (1)
%!error nakarnd (zeros (3), ones (2))
%!error nakarnd (zeros (2), ones (3))
%!error nakarnd (i, 2)
%!error nakarnd (1, i)
%!error nakarnd (1,2, -1)
%!error nakarnd (1,2, ones (2))
%!error nakarnd (1, 2, [2 -1 2])
%!error nakarnd (1,2, 1, ones (2))
%!error nakarnd (1,2, 1, -1)
%!error nakarnd (ones (2,2), 2, 3)
%!error nakarnd (ones (2,2), 2, [3, 2])
%!error nakarnd (ones (2,2), 2, 2, 3)

