## Copyright (c) 2012 Juan Pablo Carbajal <carbajal@ifi.uzh.ch>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {[@var{m}, @var{K}] =} regress_gp (@var{x}, @var{y}, @var{Sp})
## @deftypefnx {statistics} {[@dots{} @var{yi} @var{dy}] =} regress_gp (@dots{}, @var{xi})
##
## Linear scalar regression using gaussian processes.
##
## It estimates the model @var{y} = @var{x}'*m for @var{x} R^D and @var{y} in R.
## The information about errors of the predictions (interpolation/extrapolation) is given
## by the covarianve matrix @var{K}. If D==1 the inputs must be column vectors,
## if D>1 then @var{x} is n-by-D, with n the number of data points. @var{Sp} defines
## the prior covariance of @var{m}, it should be a (D+1)-by-(D+1) positive definite matrix,
## if it is empty, the default is @code{Sp = 100*eye(size(x,2)+1)}.
##
## If @var{xi} inputs are provided, the model is evaluated and returned in @var{yi}.
## The estimation of the variation of @var{yi} are given in @var{dy}.
##
## Run @code{demo regress_gp} to see an examples.
##
## The function is a direc implementation of the formulae in pages 11-12 of
## Gaussian Processes for Machine Learning. Carl Edward Rasmussen and @
## Christopher K. I. Williams. The MIT Press, 2006. ISBN 0-262-18253-X.
## available online at @url{http://gaussianprocess.org/gpml/}.
##
## @seealso{regress}
## @end deftypefn

function [wm, K, yi, dy] = regress_gp (x, y, Sp=[], xi=[])

  if isempty(Sp)
    Sp = 100*eye(size(x,2)+1);
  end

  x  = [ones(1,size(x,1)); x'];

  ## Juan Pablo Carbajal <carbajal@ifi.uzh.ch>
  ## Note that in the book the equation (below 2.11) for the A reads
  ## A  = (1/sy^2)*x*x' + inv (Vp);
  ## where sy is the scalar variance of the of the residuals (i.e y = x' * w + epsilon)
  ## and epsilon is drawn from N(0,sy^2). Vp is the variance of the parameters w.
  ## Note that
  ## (sy^2 * A)^{-1} = (1/sy^2)*A^{-1} =  (x*x' + sy^2 * inv(Vp))^{-1};
  ## and that the formula for the w mean is
  ## (1/sy^2)*A^{-1}*x*y
  ## Then one obtains
  ## inv(x*x' + sy^2 * inv(Vp))*x*y
  ## Looking at the formula bloew we see that Sp = (1/sy^2)*Vp
  ## making the regression depend on only one parameter, Sp, and not two.
  A  = x*x' + inv (Sp);
  K  = inv (A);
  wm = K*x*y;

  yi =[];
  dy =[];
  if !isempty (xi);
    xi = [ones(size(xi,1),1) xi];
    yi = xi*wm;
    dy = diag (xi*K*xi');
  end

endfunction

%!demo
%! % 1D Data
%! x = 2*rand (5,1)-1;
%! y = 2*x -1 + 0.3*randn (5,1);
%!
%! % Points for interpolation/extrapolation
%! xi = linspace (-2,2,10)';
%!
%! [m K yi dy] = regress_gp (x,y,[],xi);
%!
%! plot (x,y,'xk',xi,yi,'r-',xi,bsxfun(@plus, yi, [-dy +dy]),'b-');

%!demo
%! % 2D Data
%! x = 2*rand (4,2)-1;
%! y = 2*x(:,1)-3*x(:,2) -1 + 1*randn (4,1);
%!
%! % Mesh for interpolation/extrapolation
%! [xi yi] = meshgrid (linspace (-1,1,10));
%!
%! [m K zi dz] = regress_gp (x,y,[],[xi(:) yi(:)]);
%! zi = reshape (zi, 10,10);
%! dz = reshape (dz,10,10);
%!
%! plot3 (x(:,1),x(:,2),y,'.g','markersize',8);
%! hold on;
%! h = mesh (xi,yi,zi,zeros(10,10));
%! set(h,'facecolor','none');
%! h = mesh (xi,yi,zi+dz,ones(10,10));
%! set(h,'facecolor','none');
%! h = mesh (xi,yi,zi-dz,ones(10,10));
%! set(h,'facecolor','none');
%! hold off
%! axis tight
%! view(80,25)

%!demo
%! % Projection over basis function
%! pp = [2 2 0.3 1];
%! n = 10;
%! x = 2*rand (n,1)-1;
%! y = polyval(pp,x) + 0.3*randn (n,1);
%!
%! % Powers
%! px = [sqrt(abs(x)) x x.^2 x.^3];
%!
%! % Points for interpolation/extrapolation
%! xi = linspace (-1,1,100)';
%! pxi = [sqrt(abs(xi)) xi xi.^2 xi.^3];
%!
%! Sp = 100*eye(size(px,2)+1);
%! Sp(2,2) = 1; # We don't believe the sqrt is present
%! [m K yi dy] = regress_gp (px,y,Sp,pxi);
%! disp(m)
%!
%! plot (x,y,'xk;Data;',xi,yi,'r-;Estimation;',xi,polyval(pp,xi),'g-;True;');
%! axis tight
%! axis manual
%! hold on
%! plot (xi,bsxfun(@plus, yi, [-dy +dy]),'b-');
%! hold off


# input validation
%!test
%! x = [1; 2; 3];
%! y = [4; 5; 6];
%! [m, K] = regress_gp(x, y);
%! assert(size(m), [2, 1]);
%! assert(size(K), [2, 2]);

%!test
%! x = [1; 2; 3];
%! y = [4; 5; 6];
%! Sp = [2, 0.5; 0.5, 1];
%! [m, K] = regress_gp(x, y, Sp);
%! assert(size(m), [2, 1]);
%! assert(size(K), [2, 2]);

%!test
%! x = [1; 2; 3];
%! y = [4; 5; 6];
%! xi = [4; 5; 6];
%! [m, K, yi, dy] = regress_gp(x, y, [], xi);
%! assert(size(m), [2, 1]);
%! assert(size(K), [2, 2]);
%! assert(size(yi), [3, 1]);
%! assert(size(dy), [3, 1]);

%!test
%! x = [1; 2; 3];
%! y = [4; 5; 6];
%! Sp = [2, 0.5; 0.5, 1];
%! xi = [4; 5; 6];
%! [m, K, yi, dy] = regress_gp(x, y, Sp, xi);
%! assert(size(m), [2, 1]);
%! assert(size(K), [2, 2]);
%! assert(size(yi), [3, 1]);
%! assert(size(dy), [3, 1]);

%!test
%! x = [1; 2; 3];
%! y = [4; 5; 6];
%! Sp = [2, 0.5; 0.5, 1];
%! xi = [4, 5, 6];
%!error id=Octave:undefined-function regress_gp(x, y, Sp, xi);

%!test
%! x = [1; 2; 3];
%! y = [4; 5; 6];
%! xi = [4; 5; 6];
%!error id=Octave:undefined-function regress_gp(x, y, xi);

%!test
%! x = [1, 2; 3, 4];
%! y = [5; 6];
%!error id=Octave:undefined-function regress_gp(x, y);

%!test
%! [wm, K, yi, dy] = regress_gp([], []);
%! assert(isempty(wm));

%!test
%! x = [1 2 3; 4 5 6];
%! y = [7; 8];
%! [wm, K, yi, dy] = regress_gp(x, y);
%! assert(size(wm), [4, 1])
%! assert(size(K), [4, 4])
%! assert(isempty(yi) && isempty(dy))
%! assert(length(K(:)), 16);

%!xtest
%! x = [1 2 3; 4 5 6];
%! y = [7; 8];
%! Sp = [1 0.5 0.3; 0.5 1 0.7; 0.3 0.7 2];
%! [wm, K, yi, dy] = regress_gp(x, y, Sp);
%! assert(size(wm), [3, 1])
%! assert(size(K), [3, 3])
%! assert(isempty(yi) && isempty(dy))
%! assert(det(K) > 0);
%! assert(length(K(:)), length(unique(K(:))))

%!xtest
%! x = [1 2 3; 4 5 6];
%! y = [7; 8];
%! Sp = [1 0.5; 0.5 1];
%! [wm, K, yi, dy] = regress_gp(x, y, Sp);
%! assert(size(wm), [3, 1])
%! assert(size(K), [3, 3])
%! assert(isempty(yi) && isempty(dy))
%! assert(det(K) > 0);
%! assert(length(K(:)), length(unique(K(:))))

# output validation
%!xtest
%! x = [1;2;3];
%! y = [2;3;4];
%! Sp = [100, 0, 0; 0, 100, 0; 0, 0, 100];
%! [m, K] = regress_gp(x, y, Sp);
%! error <operator +: nonconformant arguments (op1 is 2x2, op2 is 3x3)> size(m)
%! error <operator +: nonconformant arguments (op1 is 2x2, op2 is 3x3)> size(K)

%!xtest
%! x = [1;2;3];
%! y = [2;3;4];
%! Sp = [100, 0, 0; 0, 100, 0; 0, 0, 100];
%! [m, K, yi, dy] = regress_gp(x, y, Sp, x);
%! assert( size(yi), [3, 1]);
%! assert( size(dy), [3, 1]);

%!xtest
%! x = [1 2 3; 4 5 6];
%! y = [1; 2];
%! Sp = eye(3);
%! xi = [7 8 9; 10 11 12];
%! [wm, K, yi, dy] = regress_gp(x, y, Sp, xi);
%! assert(isempty(wm));

%!xtest
%! x = [1 2; 3 4];
%! y = [1; 2; 3];
%! Sp = eye(2);
%! xi = [];
%! [wm, K, yi, dy] = regress_gp(x, y, Sp, xi);
%! assert(isempty(wm));

%!xtest
%! x = [1 2 3; 4 5 6];
%! y = [1; 2];
%! Sp = zeros(3);
%! xi = [7 8 9; 10 11 12];
%! [wm, K, yi, dy] = regress_gp(x, y, Sp, xi);
%! assert(isempty(wm));

%!test
%! x = rand(10, 2);
%! y = rand(10, 1);
%! [wm, K, yi, dy] = regress_gp(x, y);
%! assert(size(wm), [3,1]);
%! assert(size(K), [3,3]);
%! assert(isvector(yi), false);
%! assert(numel(yi), 0);
%! assert(numel(dy), 0);

% Test with empty Sp
%!test
%! x = rand(10, 2);
%! y = rand(10, 1);
%! [wm, K, yi, dy] = regress_gp(x, y, []);
%! assert(size(wm), [3,1]);
%! assert(size(K), [3,3]);
%! assert(isvector(yi), false);
%! assert(numel(yi), 0);
%! assert(numel(dy), 0);

% Test with xi
%!test
%! x = rand(10, 2);
%! y = rand(10, 1);
%! xi = rand(5, 2);
%! [wm, K, yi, dy] = regress_gp(x, y, [], xi);
%! assert(size(wm), [3,1]);
%! assert(size(K), [3,3]);
%! assert(isvector(yi));
%! assert(size(yi), [5,1]);
%! assert(size(dy), [5,1]);

% Test valid input (D=1, n=1)
%!test
%! [m, K] = regress_gp(1, 1, [], []);
%! assert(m, [0.4975; 0.4975], 1e-4);
%! assert(K, [50.249, -49.751; -49.751, 50.249], 1e-3);

% Test valid input (D=2, n=2)
%!test
%! [m, K] = regress_gp([1,2; 2,3], [1,2]', [], []);
%! assert(m, [-0.3031; 0.6406; 0.3375], 1e-4);
%! assert(K, [36.877, 30.312, -32.811; 30.312, 35.940, -33.749; -32.811, -33.749, 33.440], 1e-3);

% Test valid input (D=2, n=2, with xi input)
%!test
%! [m, K, yi, dy] = regress_gp([1,2; 2,3], [1,2]', [], [1,2; 2,3]);
%! assert(yi, [1.0125; 1.9906], 1e-4);
%! assert(dy, [0.9562; 0.9812], 1e-3);















