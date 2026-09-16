## Copyright (C) 2024-2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {Private Function} {@var{B} =} __unique__ (@var{A})
## @deftypefnx {Private Function} {@var{B} =} __unique__ (@var{A}, @var{setOrder})
## @deftypefnx {Private Function} {@var{B} =} __unique__ (@var{A}, @var{occurrence})
## @deftypefnx {Private Function} {@var{B} =} __unique__ (@var{A}, @var{setOrder}, @var{occurrence})
## @deftypefnx {Private Function} {@var{B} =} __unique__ (@var{A}, @var{occurrence}, @var{setOrder})
## @deftypefnx {Private Function} {@var{B} =} __unique__ (@var{A}, @dots{}, @qcode{'rows'})
## @deftypefnx {Private Function} {@var{B} =} __unique__ (@var{A}, @qcode{'rows'}, @dots{})
## @deftypefnx {Private Function} {[@var{B}, @var{ixA}, @var{ixB}] =} __unique__ (@dots{})
##
## Return the unique elements of @var{x}.
##
## @end deftypefn

function [y, i, j] = __unique__ (x, varargin)

  if (nargin < 1)
    print_usage ();
  elseif (! (isnumeric (x) || islogical (x) || ischar (x) || iscellstr (x)))
    error ("unique: X must be an array or cell array of strings.");
  endif

  if (nargin > 1)
    ## Parse options
    if (! iscellstr (varargin))
      error ("unique: options must be strings.");
    endif

    optrows   = any (strcmp ("rows", varargin));
    optfirst  = any (strcmp ("first", varargin));
    optlast   = any (strcmp ("last", varargin));
    optsorted = any (strcmp ("sorted", varargin));
    optstable = any (strcmp ("stable", varargin));
    if (optrows && ndims (x) != 2)
      error ("unique: 'rows' applies only to 2-D matrices.");
    elseif (optfirst && optlast)
      error ("unique: cannot specify both 'first' and 'last'.");
    elseif (optsorted && optstable)
      error ("unique: cannot specify both 'sorted' and 'stable'.");
    elseif (optrows + optfirst + optlast + optsorted + optstable != nargin-1)
      error ("unique: invalid optional argument.");
    endif

    ## Set defaults if not set earlier.
    if (! optfirst && ! optlast)
      optfirst = true;
    endif
    if (! optsorted && ! optstable)
      optsorted = true;
    endif
  else
    optrows = false;
    optfirst = true;
    optsorted = true;
  endif

  if (optrows)
    n = rows (x);
    isrowvec = false;
  else
    n = numel (x);
    isrowvec = isrow (x);
  endif

  ## Special cases 0 and 1
  if (n == 0)
    y = x;
    if (! optrows && any (size (x)))
      if (iscellstr (x))
        y = cell (0, 1);
      else
        y = zeros (0, 1, class (x));
      endif
    endif
    i = j = [];
    return;
  elseif (n == 1)
    y = x;
    i = j = 1;
    return;
  endif

  ## Calculate y output
  if (optrows)
    if (nargout > 1 || ! optsorted)
      [y, j] = sortrows (x);
      j = j(:);
    else
      y = sortrows (x);
    endif
    if (iscellstr (x))
      match = all (cellfun (@isequal, y(1:n-1,:), y(2:n,:)), 2);
    else
      match = all (y(1:n-1,:) == y(2:n,:), 2);
    endif
    if (optsorted)
      y(match,:) = [];
    else
      y = x;
      if (optfirst)
        y(j([false; match]), :) = [];
      else
        y(j([match; false]), :) = [];
      endif
    endif
  else
    if (isvector (x))
      y = x;
    else
      y = x(:);
    endif
    if (nargout > 1 || ! optsorted)
      [y, j] = sort (y);
      j = j(:);
    else
      y = sort (y);
    endif
    if (iscellstr (y))
      match = strcmp (y(1:n-1), y(2:n));
    else
      match = (y(1:n-1) == y(2:n));
    endif
    if (optsorted)
      y(match) = [];
    else
      if (isvector (x))
        y = x;
      else
        y = x(:);
      endif
      if (optfirst)
        y(j([false; match(:)])) = [];
      else
        y(j([match(:); false])) = [];
      endif
    endif
  endif


  ## Calculate 2nd and 3rd outputs
  if (nargout > 1)
    if (optsorted)
      idx = find (match);

      if (optfirst)
        idx += 1;
      endif

      i = j;
      i(idx) = [];

      if (nargout > 2)
        j(j) = cumsum (! [false; match(:)]);
      endif

    else
      ## Get inverse of sort index j so that sort(x)(k) = x(j)(k) = x.
      k = j;  # cheap way to copy dimensions
      k(j) = 1:n;

      ## Generate logical index of sorted unique value locations.
      if (optfirst)
        uniquex = ! [false; match(:)];
      else
        uniquex = ! [match(:); false];
      endif

      ## Remap unique locations to unsorted x, such that y = x(i).
      i = find (uniquex(k));

      if (nargout > 2)

        ni = numel (i);

        u = find (uniquex);   # Linear index of unique elements of sort(x)
        p = j;                # cheap way to copy dimensions

        if (optfirst)
          l = u(cumsum (uniquex));  # Expand u for all elements in sort(x)
          p(i) = 1:ni;              # set p to contain the vector positions of i
        elseif (optrows)
          l = u(cumsum (uniquex, 'reverse'));
          p(i) = ni:-1:1;
        else
          l = u(cumsum (! [false; match(:)]));
          p(i) = 1:ni;
        endif

        j = p(j(l(k))); # Replace j with 3rd output mapping y->x.

      endif
    endif
  endif

endfunction
