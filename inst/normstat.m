## Copyright (C) 2006 Arno Onken
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{mn}, @var{vr}] =} normstat (@var{m}, @var{v})
## Returns mean and variance of the normal distribution, the given arguments
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{m} is the mean of the normal distribution
##
## @item
## @var{v} is the variance of the normal distribution.
## @var{v} must be positive
## @end itemize
## @var{m} and @var{v} must be of common size or one of them must be
## scalar
##
## @subheading Return values
##
## @itemize @bullet
## @item
## @var{mn} is the mean of the normal distribution
##
## @item
## @var{vr} is the variance of the normal distribution
## @end itemize
##
## @subheading Examples
##
## @example
## @group
## m = 1:6;
## v = 0:0.2:1;
## [mn, vr] = normstat (m, v)
## @end group
##
## @group
## [mn, vr] = normstat (0, v)
## @end group
## @end example
##
## @subheading References
##
## @enumerate
## @item
## Wendy L. Martinez and Angel R. Martinez. @cite{Computational Statistics
## Handbook with MATLAB}. Appendix E, pages 547-557, Chapman & Hall/CRC,
## 2001.
##
## @item
## Athanasios Papoulis. @cite{Probability, Random Variables, and Stochastic
## Processes}. McGraw-Hill, New York, second edition, 1984.
## @end enumerate
## @end deftypefn

## Author: Arno Onken <asnelt@asnelt.org>
## Description: Moments of the normal distribution

function [mn, vr] = normstat (m, v)

  # Check arguments
  if (nargin != 2)
    usage ("[mn, vr] = normstat (m, v)");
  endif

  if (! isempty (m) && ! ismatrix (m))
    error ("normstat: m must be a numeric matrix");
  endif
  if (! isempty (v) && ! ismatrix (v))
    error ("normstat: v must be a numeric matrix");
  endif

  if (! isscalar (m) || ! isscalar (v))
    [retval, m, v] = common_size (m, v);
    if (retval > 0)
      error ("normstat: m and v must be of common size or scalar");
    endif
  endif

  # Set moments
  mn = m;
  vr = v;

  # Continue argument check
  k = find (! (v > 0) | ! (v < Inf));
  if (any (k))
    mn (k) = NaN;
    vr (k) = NaN;
  endif

endfunction

%!test
%! m = 1:6;
%! v = 0.2:0.2:1.2;
%! [mn, vr] = normstat (m, v);
%! assert (mn, m);
%! assert (vr, v);

%!test
%! v = 0.2:0.2:1.2;
%! [mn, vr] = normstat (0, v);
%! expected_mn = [0, 0, 0, 0, 0, 0];
%! assert (mn, expected_mn);
%! assert (vr, v);
