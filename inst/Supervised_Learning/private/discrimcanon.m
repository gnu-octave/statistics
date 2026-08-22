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
## @deftypefn  {statistics} {@var{t} =} discrimcanon (@var{str})
## @deftypefnx {statistics} {[@var{t}, @var{fam}] =} discrimcanon (@var{str})
##
## Canonical spelling and family of a discriminant type.
##
## @var{t} is @var{str} respelled as the discriminant classes report it, which
## is case sensitive: @qcode{'diagLinear'}, not @qcode{'diaglinear'}.  It is
## empty when @var{str} names no discriminant type, which is how the callers
## detect an invalid value.
##
## @var{fam} is @qcode{'linear'} or @qcode{'quadratic'}, the family @var{t}
## belongs to, and empty alongside an empty @var{t}.  The family decides which
## covariances a fit has to keep: a linear one pools across classes and a
## quadratic one estimates a covariance per class, and no assignment ever moves
## a model between the two.
##
## @end deftypefn

function [t, fam] = discrimcanon (str)

  lin  = {'linear', 'diagLinear', 'pseudoLinear'};
  quad = {'quadratic', 'diagQuadratic', 'pseudoQuadratic'};

  t = '';
  fam = '';
  if (! (ischar (str) && isrow (str)))
    return;
  endif

  idx = find (strcmpi (str, lin));
  if (! isempty (idx))
    t = lin{idx};
    fam = 'linear';
    return;
  endif

  idx = find (strcmpi (str, quad));
  if (! isempty (idx))
    t = quad{idx};
    fam = 'quadratic';
  endif

endfunction
