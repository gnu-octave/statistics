## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Function File} [@var{err}, @var{varargout}] = fixvarsize (@var{varargin})
##
## Make all input arguments the same size.
##
## @code{fixvarsize} gets the size of each input argument and returns them
## according to the size of the first non-scalar input argument.  All scalar
## input arguments receive the size of the first non-scalar input argument, if
## all non-scalar input arguments share the same size.  In such case, the output
## argument @var{err} is 0.  In any other case, @var{err} equals 1, and all
## input arguments are returned in their original size.
##
## @code{fixvarsize} resembles the functionality of the MATLAB @code{distchck}
## function with somewhat different syntax, i.e. it does not require the number
## of input variables as the first inpuut argument.
##
## @end deftypefn

function [err, varargout] = fixvarsize (varargin)
  err = 0;
  varargout = varargin;
  ## Check for all input arguments being scalars
  is_scalar = (cellfun ("prodofsize", varargin) == 1);
  if (all (is_scalar))
    return;
  endif
  ## Get size of each input argument
  for j = 1:numel (varargin)
    sz{j} = size (varargin {j});
  endfor
  ## Find first non-scalar input argument
  t = sz(! is_scalar);
  size1 = t{1};
  ## Scalars receive the size of the first non-scalar input argument
  ## Other non-scalar input arguments must be of the same size.
  for j = 1:numel (varargin)
    sizej = sz{j};
    if (is_scalar (j))                  # scalars
      vj = varargin{j};
      if (isnumeric (vj))
        t = zeros (size1, class (vj));
      else
        t = zeros (size1);
      end
      t(:) = varargin{j};
      varargout{j} = t;
    elseif (! isequal (sizej, size1))   # arrays
      err = 1;
      varargout = varargin;
      return;
    end
  end
endfunction

## Output validation tests
%!test
%! [err, a, b, c] = fixvarsize (1, 2, [2, 3, 4]);
%! assert (err, 0);
%! assert (size (a), size (c));
%! assert (size (b), size (c));
%! assert (size (a), [1, 3]);
%!test
%! [err, a, b, c] = fixvarsize (1, [2; 2], [2, 3, 4]);
%! assert (err, 1);
%! assert (size (a), [1, 1]);
%! assert (size (b), [2, 1]);
%! assert (size (c), [1, 3]);
%!test
%! [err, a, b] = fixvarsize (1, ones (3,4,5));
%! assert (err, 0);
%! assert (size (a), [3, 4, 5]);
%! assert (size (b), [3, 4, 5]);
%!test
%! [err, a, b] = fixvarsize (zeros (2, 3), ones (3,4,5));
%! assert (err, 1);
%! assert (size (a), [2, 3]);
%! assert (size (b), [3, 4, 5]);
