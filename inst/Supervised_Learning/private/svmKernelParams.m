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
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Private Function} {@var{s} =} svmKernelParams (@var{fcn}, @var{scale}, @var{order})
##
## Build the @qcode{KernelParameters} structure of an SVM model.
##
## @var{fcn} is the kernel as the class recorded it, @var{scale} its scale and
## @var{order} the polynomial order.  The structure carries @qcode{Function}
## and @qcode{Scale}, and @qcode{Order} for a polynomial kernel only, which is
## the shape MATLAB reports.  MATLAB names the radial basis kernel
## @qcode{'gaussian'} whatever the caller spelled it, so @qcode{'rbf'} is
## renamed here; the kernel the fit was given is unchanged and remains in
## @qcode{ModelParameters}.
##
## @end deftypefn

function s = svmKernelParams (fcn, scale, order)

  if (strcmpi (fcn, 'rbf'))
    fcn = 'gaussian';
  endif
  s = struct ('Function', fcn, 'Scale', scale);
  if (strcmpi (fcn, 'polynomial'))
    s.Order = order;
  endif

endfunction
