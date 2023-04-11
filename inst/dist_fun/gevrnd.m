## Copyright (C) 2012 Nir Krakauer <nkrakauer@ccny.cuny.edu>
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
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
## @deftypefn  {statistics} {@var{r} =} gevrnd (@var{k}, @var{sigma}, @var{mu})
## @deftypefnx {statistics} {@var{r} =} gevrnd (@var{k}, @var{sigma}, @var{mu}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} gevrnd (@var{k}, @var{sigma}, @var{mu}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} gevrnd (@var{k}, @var{sigma}, @var{mu}, [@var{sz}])
##
## Random arrays from the generalized extreme value (GEV) distribution.
##
## @code{@var{r} = gevrnd (@var{k}, @var{sigma}, @var{mu}} returns an array of
## random numbers chosen from the GEV distribution with shape parameter @var{k},
## scale parameter @var{sigma}, and location parameter @var{mu}.  The size of
## @var{r} is the common size of the input arguments.  A scalar input functions
## as a constant matrix of the same size as the other inputs.
##
## When called with a single size argument, returns a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector @var{sz} of dimensions.
##
## @subheading References
##
## @enumerate
## @item
## Rolf-Dieter Reiss and Michael Thomas. @cite{Statistical Analysis of Extreme
## Values with Applications to Insurance, Finance, Hydrology and Other Fields}.
## Chapter 1, pages 16-17, Springer, 2007.
## @end enumerate
##
## @seealso{gevcdf, gevinv, gevpdf, gevfit, gevlike, gevstat}
## @end deftypefn

function r = gevrnd (k, sigma, mu, varargin)

  if (nargin < 3)
    print_usage ();
  endif

  if any (sigma <= 0)
    error ("gevrnd: sigma must be positive.");
  endif

  if (!isscalar (k) || !isscalar (sigma) || !isscalar (mu))
    [retval, k, sigma, mu] = common_size (k, sigma, mu);
    if (retval > 0)
      error ("gevrnd: k, sigma, mu must be of common size or scalars.");
    endif
  endif

  if (iscomplex (k) || iscomplex (sigma) || iscomplex (mu))
    error ("gevrnd: k, sigma, mu must not be complex");
  endif

  if (nargin == 3)
    sz = size (k);
  elseif (nargin == 4)
    if (isscalar (varargin{1}) && varargin{1} >= 0)
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0))
      sz = varargin{1};
    else
      error (strcat (["gevrnd: dimension vector must be row vector of"], ...
                     [" non-negative integers."]));
    endif
  elseif (nargin > 4)
    if (any (cellfun (@(x) (!isscalar (x) || x < 0), varargin)))
      error ("gevrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  if (!isscalar (k) && !isequal (size (k), sz))
    error ("gevrnd: k, sigma, mu must be scalar or of size SZ.");
  endif

  if (isa (k, "single") || isa (sigma, "single") || isa (mu, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  r = gevinv (rand(sz), k, sigma, mu);

  if (strcmp (cls, "single"))
    r = single (r);
  endif

endfunction


%!assert(size (gevrnd (1,2,1)), [1, 1]);
%!assert(size (gevrnd (ones(2,1), 2, 1)), [2, 1]);
%!assert(size (gevrnd (ones(2,2), 2, 1)), [2, 2]);
%!assert(size (gevrnd (1, 2*ones(2,1), 1)), [2, 1]);
%!assert(size (gevrnd (1, 2*ones(2,2), 1)), [2, 2]);
%!assert(size (gevrnd (1, 2, 1, 3)), [3, 3]);
%!assert(size (gevrnd (1, 2, 1, [4 1])), [4, 1]);
%!assert(size (gevrnd (1, 2, 1, 4, 1)), [4, 1]);

%% Test input validation
%!error gevrnd ()
%!error gevrnd (1, 2)
%!error gevrnd (ones(3),ones(2),1)
%!error gevrnd (ones(2),ones(3),1)
%!error gevrnd (i, 2, 1)
%!error gevrnd (2, i, 1)
%!error gevrnd (2, 0, 1)
%!error gevrnd (1,2, 1, -1)
%!error gevrnd (1,2, 1,  ones(2))
%!error gevrnd (1,2, 1,  [2 -1 2])
%!error gevrnd (1,2, 1,  1, ones(2))
%!error gevrnd (1,2, 1,  1, -1)
%!error gevrnd (ones(2,2), 2, 1, 3)
%!error gevrnd (ones(2,2), 2, 1, [3, 2])
%!error gevrnd (ones(2,2), 2, 1, 2, 3)

