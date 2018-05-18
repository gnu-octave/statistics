## Copyright (C) 2017 - Juan Pablo Carbajal
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

## Author: Juan Pablo Carbajal <ajuanpi+dev@gmail.com>

## -*- texinfo -*-
## @deftypefn {} {@var{pts} =} sigma_pts (@var{n})
## @deftypefnx {@var{pts} =} sigma_pts (@var{n}, @var{m})
## @deftypefnx {@var{pts} =} sigma_pts (@var{n}, @var{m}, @var{K})
## @deftypefnx {@var{pts} =} sigma_pts (@var{n}, @var{m}, @var{K}, @var{l})
## Calculates 2*@var{n}+1 sigma points in @var{n} dimensions.
##
## Sigma points are used in the unscented transfrom to estimate
## the result of applying a given nonlinear transformation to a probability
## distribution that is characterized only in terms of a finite set of statistics.
##
## If only the dimension @var{n} is given the resulting points have zero mean
## and identity covariance matrix.
## If the mean @var{m} or the covaraince matrix @var{K} are given, then the resulting points
## will have those statistics.
## The factor @var{l} scaled the points away from the mean. It is useful to tune
## the accuracy of the unscented transfrom.
##
## There is no unique way of computing sigma points, this function implements the
## algorithm described in section 2.6 "The New Filter" pages 40-41 of
##
## Uhlmann, Jeffrey (1995). "Dynamic Map Building and Localization: New Theoretical Foundations".
## Ph.D. thesis. University of Oxford.
##
## @end deftypefn

function pts = sigma_pts (n, m = [], K = [], l = 0)

  if isempty (K)
    K = eye (n);
  endif
  if isempty (m)
    m = zeros (1, n);
  endif

  if (n ~= length (m))
    error ("Dimension and size of mean vector don't match.")
  endif
  if any(n ~= size (K))
    error ("Dimension and size of covariance matrix don't match.")
  endif

  if isdefinite (K) <= 0
    error ("Covariance matrix should be positive definite.")
  endif

  pts      = zeros (2 * n + 1, n);
  pts(1,:) = m;

  K              = sqrtm ((n + l) * K);
  pts(2:n+1,:)   = bsxfun (@plus, m , K);
  pts(n+2:end,:) = bsxfun (@minus, m , K);

endfunction

%!demo
%! K      = [1 0.5; 0.5 1]; # covaraince matrix
%! # calculate and build associated ellipse
%! [R,S,~] = svd (K);
%! theta   = atan2 (R(2,1), R(1,1));
%! v       = sqrt (diag (S));
%! v       = v .* [cos(theta) sin(theta); -sin(theta) cos(theta)];
%! t       = linspace (0, 2*pi, 100).';
%! xe      = v(1,1) * cos (t) + v(2,1) * sin (t);
%! ye      = v(1,2) * cos (t) + v(2,2) * sin (t);
%!
%! figure(1); clf; hold on
%! # Plot ellipse and axes
%! line ([0 0; v(:,1).'],[0 0; v(:,2).'])
%! plot (xe,ye,'-r');
%!
%! col = 'rgb';
%! l     = [-1.8 -1 1.5];
%! for li = 1:3
%!  p     = sigma_pts (2, [], K, l(li));
%!  tmp   = plot (p(2:end,1), p(2:end,2), ['x' col(li)], ...
%!               p(1,1), p(1,2), ['o' col(li)]);
%!  h(li) = tmp(1);
%! endfor
%! hold off
%! axis image
%! legend (h, arrayfun (@(x) sprintf ("l:%.2g", x), l, "unif", 0));


%!test
%! p = sigma_pts (5);
%! assert (mean (p), zeros(1,5), sqrt(eps));
%! assert (cov (p), eye(5), sqrt(eps));

%!test
%! m = randn(1, 5);
%! p = sigma_pts (5, m);
%! assert (mean (p), m, sqrt(eps));
%! assert (cov (p), eye(5), sqrt(eps));

%!test
%! x = linspace (0,1,5);
%! K = exp (- (x.' - x).^2/ 0.5);
%! p = sigma_pts (5, [], K);
%! assert (mean (p), zeros(1,5), sqrt(eps));
%! assert (cov (p), K, sqrt(eps));

%!error sigma_pts(2,1);
%!error sigma_pts(2,[],1);
%!error sigma_pts(2,1,1);
%!error sigma_pts(2,[0.5 0.5],[-1 0; 0 0]);
