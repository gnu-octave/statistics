## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} @var{r} = nbinrnd (@var{n}, @var{ps})
## @deftypefnx {statistics} @var{r} = nbinrnd (@var{n}, @var{ps}, @var{rows})
## @deftypefnx {statistics} @var{r} = nbinrnd (@var{n}, @var{ps}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} @var{r} = nbinrnd (@var{n}, @var{ps}, [@var{sz}])
##
## Random arrays from the negative binomial distribution.
##
## @code{@var{r} = nbinrnd (@var{n}, @var{ps})} returns an array of random
## numbers chosen from the Laplace distribution with parameters @var{n} and
## @var{ps}.  The size of @var{r} is the common size of @var{n} and @var{ps}.
## A scalar input functions as a constant matrix of the same size as the other
## inputs.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## @seealso{nbininv, nbininv, nbinpdf, nbinstat}
## @end deftypefn

function rnd = nbinrnd (n, ps, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    print_usage ();
  endif

  ## Check for common size N and PS
  if (! isscalar (n) || ! isscalar (ps))
    [retval, n, ps] = common_size (n, ps);
    if (retval > 0)
      error ("nbinrnd: N and PS must be of common size or scalars.");
    endif
  endif

  ## Check for N and PS being reals
  if (iscomplex (n) || iscomplex (ps))
    error ("nbinrnd: N and PS must not be complex.");
  endif

  ## Check for SIZE vector or DIMENSION input arguments
  if (nargin == 2)
    sz = size (n);
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} >= 0)
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0))
      sz = varargin{1};
    else
      error (strcat (["nbinrnd: dimension vector must be row vector"], ...
                     [" of non-negative integers."]));
    endif
  elseif (nargin > 3)
    if (any (cellfun (@(x) (! isscalar (x) || x < 0), varargin)))
      error ("nbinrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (n) && ! isequal (size (n), sz))
    error ("nbinrnd: N and PS must be scalar or of size SZ.");
  endif

  ## Check for appropriate class
  if (isa (n, "single") || isa (ps, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Generate random sample from negative binomial distribution
  if (isscalar (n) && isscalar (ps))
    if ((n > 0) && (n < Inf) && (ps > 0) && (ps <= 1))
      rnd = randp ((1 - ps) ./ ps .* randg (n, sz, cls), cls);
    elseif ((n > 0) && (n < Inf) && (ps == 0))
      rnd = zeros (sz, cls);
    else
      rnd = NaN (sz, cls);
    endif
  else
    rnd = NaN (sz, cls);

    k = (n > 0) & (n < Inf) & (ps == 0);
    rnd(k) = 0;

    k = (n > 0) & (n < Inf) & (ps > 0) & (ps <= 1);
    rnd(k) = randp ((1 - ps(k)) ./ ps(k) .* randg (n(k), cls));
  endif

endfunction


%!assert (size (nbinrnd (2, 1/2)), [1, 1])
%!assert (size (nbinrnd (2*ones (2,1), 1/2)), [2, 1])
%!assert (size (nbinrnd (2*ones (2,2), 1/2)), [2, 2])
%!assert (size (nbinrnd (2, 1/2*ones (2,1))), [2, 1])
%!assert (size (nbinrnd (2, 1/2*ones (2,2))), [2, 2])
%!assert (size (nbinrnd (2, 1/2, 3)), [3, 3])
%!assert (size (nbinrnd (2, 1/2, [4 1])), [4, 1])
%!assert (size (nbinrnd (2, 1/2, 4, 1)), [4, 1])

## Test class of input preserved
%!assert (class (nbinrnd (2, 1/2)), "double")
%!assert (class (nbinrnd (single (2), 1/2)), "single")
%!assert (class (nbinrnd (single ([2 2]), 1/2)), "single")
%!assert (class (nbinrnd (2, single (1/2))), "single")
%!assert (class (nbinrnd (2, single ([1/2 1/2]))), "single")

## Test input validation
%!error nbinrnd ()
%!error nbinrnd (1)
%!error nbinrnd (ones (3), ones (2))
%!error nbinrnd (ones (2), ones (3))
%!error nbinrnd (i, 2)
%!error nbinrnd (2, i)
%!error nbinrnd (1,2, -1)
%!error nbinrnd (1,2, ones (2))
%!error nbinrnd (1, 2, [2 -1 2])
%!error nbinrnd (1,2, 1, ones (2))
%!error nbinrnd (1,2, 1, -1)
%!error nbinrnd (ones (2,2), 2, 3)
%!error nbinrnd (ones (2,2), 2, [3, 2])
%!error nbinrnd (ones (2,2), 2, 2, 3)
