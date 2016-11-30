function [A,B,r,U,V] = canoncorr (X,Y)

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{A} @var{B} @var{r} @var{U} @var{V}] =} canoncorr (@var{X}, @var{Y})
## Canonical correlation analysis
##
## Given @var{X} (size @var{k}*@var{m}) and @var{Y} (@var{k}*@var{n}), returns projection matrices of canonical coefficients @var{A} (size @var{m}*@var{d}, where @var{d}=@code{min}(@var{m}, @var{n})) and @var{B} (size @var{m}*@var{d}); the canonical correlations @var{r} (1*@var{d}, arranged in decreasing order); and the canonical variables @var{U}, @var{V} (both @var{k}*@var{d}, with orthonormal columns).
##
## Reference: William H. Press (2011), Canonical Correlation Clarified by Singular Value Decomposition, http://numerical.recipes/whp/notes/CanonCorrBySVD.pdf
##
## @seealso{princomp}
## @end deftypefn

#	Copyright (C) 2016 by Nir Krakauer <mail@nirkrakauer.net>

#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; If not, see <http://www.gnu.org/licenses/>.


k = size (X, 1); #should also be size (Y, 1)
d = min (size(X, 2), size(Y, 2));

X = center (X);
Y = center (Y);

[Qx Rx] = qr (X, 0);
[Qy Ry] = qr (Y, 0);

[U S V] = svd (Qx' * Qy, 0);

A = Rx \ U(:, 1:d);
B = Ry \ V(:, 1:d);

#A, B are scaled to make the covariance matrices of the outputs U, V identity matrices
f = sqrt (k-1);
A .*= f;
B .*= f;

if isargout (3)
  r = max(0, min(diag(S), 1))';
endif
if isargout (4)
  U = X * A;
endif
if isargout (5)
  V = Y * B;
endif

%!shared X,Y,A,B,r,U,V
%! X = [1:10; sin(1:10); cos(1:10)]'; Y = [tan(1:10); tanh((1:10)/10)]';
%! [A,B,r,U,V] = canoncorr (X,Y);
%!assert (A, [-0.329229   0.072908; 0.074870   1.389318; -0.069302  -0.024109], 1E-6);
%!assert (B, [-0.017086  -0.398402; -4.475049  -0.824538], 1E-6);
%!assert (r, [0.99590   0.26754], 1E-5);
%!assert (U, center(X) * A, 10*eps);
%!assert (V, center(Y) * B, 10*eps);
%!assert (cov(U), eye(size(U, 2)), 10*eps);
%!assert (cov(V), eye(size(V, 2)), 10*eps);

