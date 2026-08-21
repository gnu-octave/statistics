## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Private Function} parseResponseTransform (@var{rt}, @var{classname})
## Parse ResponseTransform for Regression objects.
## @end deftypefn

function [f, rt] = parseResponseTransform (ResponseTransform, classname)

  if (is_function_handle (ResponseTransform))
    ## nargin () raises on a handle to a built-in, so the handle is checked by
    ## what it does rather than by its declared arity.
    v = (1:5)';
    if (! isequal (size (v), size (ResponseTransform (v))))
      error (strcat ("%s: function handle for 'ResponseTransform' must", ...
                     " return the same size as its input."), classname);
    endif
    f = ResponseTransform;
    rt = 'custom function handle';
  elseif (ischar (ResponseTransform) && isrow (ResponseTransform))
    rt = tolower (ResponseTransform);
    switch (rt)
      case {'none', 'identity'}
        f = @(y) y;
        rt = 'none';
      case 'exp'
        f = @(y) exp (y);
      case 'log'
        f = @(y) log (y);
      otherwise
        error ("%s: unrecognized 'ResponseTransform' function.", classname);
    endswitch
  else
    error (strcat ("%s: 'ResponseTransform' must be a character vector", ...
                   " or a function handle."), classname);
  endif

endfunction
