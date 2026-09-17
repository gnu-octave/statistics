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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Private Function} {@var{ice} =} pdIce (@var{Mdl}, @var{F}, @var{predArgs})
##
## What a model answers for each observation on its own, one query point of the
## varied predictor at a time.
##
## These are the individual conditional expectation curves, one per row of the
## observations in @var{F}, and they are what
## @code{plotPartialDependence} draws in grey.  They come from @code{predict}
## for every model, a decision tree included: a curve belongs to one
## observation, so there is no distribution to walk.
##
## The line drawn over them is their mean, which for a tree, a bagged ensemble
## or a generalized additive model is **not** what @code{partialDependence}
## returns unless the observations are the ones the model was fitted on.
## MATLAB R2024a draws the mean of the curves here and the walked value where
## @qcode{'Conditional'} is @qcode{'none'}, measured 2026-09-17, and this
## follows it.
##
## @var{ice} is @math{nxnumX} for a regression model and
## @math{nxnumXxnum} for a classifier, the last dimension holding a page per
## class named.  Only one predictor may be varied, which is what MATLAB
## requires of a conditional plot.
##
## @end deftypefn

function ice = pdIce (Mdl, F, predArgs)

  if (nargin < 3)
    predArgs = {};
  endif

  qx = F.QP{1};
  nx = numel (qx);
  v1 = F.Vars(1);
  n = rows (F.X);
  if (F.IsClass)
    m = numel (F.LabelIdx);
  else
    m = 1;
  endif

  ice = zeros (n, nx, m);
  Z = F.X;
  for ii = 1:nx
    Z(:,v1) = qx(ii);
    if (F.IsClass)
      [~, s] = predict (Mdl, Z, predArgs{:});
      ice(:,ii,:) = reshape (s(:, F.LabelIdx), n, 1, m);
    else
      s = predict (Mdl, Z, predArgs{:});
      ice(:,ii) = s(:);
    endif
  endfor

endfunction
