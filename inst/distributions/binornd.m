## Copyright (C) 2015 Michael Leitner <michael.leitner@frm2.tum.de>
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
## @deftypefn  {statistics} {} binornd (@var{n}, @var{ps})
## @deftypefnx {statistics} {} binornd (@var{n}, @var{ps}, @var{r})
## @deftypefnx {statistics} {} binornd (@var{n}, @var{ps}, @var{r}, @var{c}, @dots{})
## @deftypefnx {statistics} {} binornd (@var{n}, @var{ps}, [@var{sz}])
##
## Random arrays from the Binomial distribution
##
## Return a matrix of random samples from the binomial distribution with
## parameters @var{n} and @var{ps}, where @var{n} is the number of trials
## and @var{ps} is the probability of success.  The size of @var{r} is the
## common size of @var{n} and @var{ps}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## If no size arguments are given then the result matrix is the common size of
## @var{n} and @var{ps}.
##
## @seealso{binocdf, binoinv, binopdf, binostat, binotest}
## @end deftypefn

function r = binornd (n, ps, varargin)

  if (nargin < 2)
    print_usage ();
  endif

  if (! isscalar (n) || ! isscalar (ps))
    [retval, n, ps] = common_size (n, ps);
    if (retval > 0)
      error ("binornd: N and PS must be of common size or scalars.");
    endif
  endif

  if (iscomplex (n) || iscomplex (ps))
    error ("binornd: N and PS must not be complex.");
  endif

  if (nargin == 2)
    sz = size (n);
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} >= 0)
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0))
      sz = varargin{1};
    else
      error (strcat (["binornd: dimension vector must be row vector of"], ...
                     [" non-negative integers."]));
    endif
  elseif (nargin > 3)
    if (any (cellfun (@(x) (! isscalar (x) || x < 0), varargin)))
      error ("binornd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  if (! isscalar (n) && ! isequal (size (n), sz))
    error ("binornd: N and PS must be scalar or of size SZ.");
  endif

  if (isa (n, "single") || isa (ps, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  if (isscalar (n) && isscalar (ps))
    if ((n > 0) && (n < Inf) && (n == fix (n)) && (ps >= 0) && (ps <= 1))
      nel = prod (sz);
      tmp = rand (n, nel);
      r = sum (tmp < ps, 1);
      r = reshape (r, sz);
      if (strcmp (cls, "single"))
        r = single (r);
      endif
    elseif ((n == 0) && (ps >= 0) && (ps <= 1))
      r = zeros (sz, cls);
    else
      r = NaN (sz, cls);
    endif
  else
    r = zeros (sz, cls);

    k = !(n >= 0) | !(n < Inf) | !(n == fix (n)) | !(ps >= 0) | !(ps <= 1);
    r(k) = NaN;

    k = (n > 0) & (n < Inf) & (n == fix (n)) & (ps >= 0) & (ps <= 1);
    if (any (k(:)))
      L = sum (k(:));
      ind = repelems ((1 : L), [(1 : L); n(k)(:)'])';
      p_ext = ps(k)(ind)(:);
      r(k) = accumarray (ind, rand (sum(n(k)(:)), 1) < p_ext);
    endif
  endif

endfunction


%!assert (binornd (0, 0, 1), 0)
%!assert (binornd ([0, 0], [0, 0], 1, 2), [0, 0])

%!assert (size (binornd (2, 1/2)), [1, 1])
%!assert (size (binornd (2*ones (2,1), 1/2)), [2, 1])
%!assert (size (binornd (2*ones (2,2), 1/2)), [2, 2])
%!assert (size (binornd (2, 1/2*ones (2,1))), [2, 1])
%!assert (size (binornd (2, 1/2*ones (2,2))), [2, 2])
%!assert (size (binornd (2, 1/2, 3)), [3, 3])
%!assert (size (binornd (2, 1/2, [4 1])), [4, 1])
%!assert (size (binornd (2, 1/2, 4, 1)), [4, 1])

## Test class of input preserved
%!assert (class (binornd (2, 0.5)), "double")
%!assert (class (binornd (single (2), 0.5)), "single")
%!assert (class (binornd (single ([2 2]), 0.5)), "single")
%!assert (class (binornd (2, single (0.5))), "single")
%!assert (class (binornd (2, single ([0.5 0.5]))), "single")

## Test input validation
%!error binornd ()
%!error binornd (1)
%!error binornd (ones (3), ones (2))
%!error binornd (ones (2), ones (3))
%!error binornd (i, 2)
%!error binornd (2, i)
%!error binornd (1,2, -1)
%!error binornd (1,2, ones (2))
%!error binornd (1,2, [2 -1 2])
%!error binornd (1,2, 1, ones (2))
%!error binornd (1,2, 1, -1)
%!error binornd (ones (2,2), 2, 3)
%!error binornd (ones (2,2), 2, [3, 2])
%!error binornd (ones (2,2), 2, 2, 3)
