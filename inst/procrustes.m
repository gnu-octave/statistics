## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/OR
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, OR (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{d} =} procrustes (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{d} =} procrustes (@var{X}, @var{Y}, @var{param1}, @var{value1}, @dots{})
## @deftypefnx {statistics} {[@var{d}, @var{Z}] =} procrustes (@dots{})
## @deftypefnx {statistics} {[@var{d}, @var{Z}, @var{transform}] =} procrustes (@dots{})
##
## Procrustes Analysis.
##
## @code{@var{d} = procrustes (@var{X}, @var{Y})} computes a linear
## transformation of the points in the matrix @var{Y} to best conform them to
## the points in the matrix @var{X} by minimizing the sum of squared errors, as
## the goodness of fit criterion, which is returned in @var{d} as a
## dissimilarity measure.  @var{d} is standardized by a measure of the scale of
## @var{X}, given by
## @itemize
## @item @qcode{sum (sum ((X - repmat (mean (X, 1), size (X, 1), 1)) .^ 2, 1))}
## @end itemize
## i.e., the sum of squared elements of a centered version of @var{X}.  However,
## if @var{X} comprises repetitions of the same point, the sum of squared errors
## is not standardized.
##
## @var{X} and @var{Y} must have the same number of points (rows) and
## @qcode{procrustes} matches the @math{i}-th point in @var{Y} to the
## @math{i}-th point in @var{X}.  Points in @var{Y} can have smaller dimensions
## (columns) than those in @var{X}, but not the opposite.  Missing dimensions in
## @var{Y} are added with padding columns of zeros as necessary to match the
## the dimensions in @var{X}.
##
## @code{[@var{d}, @var{Z}] = procrustes (@var{X}, @var{Y})} also returns the
## transformed values in @var{Y}.
##
## @code{[@var{d}, @var{Z}, @var{transform}] = procrustes (@var{X}, @var{Y})}
## also returns the transformation that maps @var{Y} to @var{Z}.
##
## @var{transform} is a structure with fields:
##
## @multitable @columnfractions 0.05 0.1 0.05 0.8
## @item @tab @qcode{c} @tab @tab the translation component
## @item @tab @qcode{T} @tab @tab the orthogonal rotation and reflection
## component
## @item @tab @qcode{b} @tab @tab the scale component
## @end multitable
##
## So that @code{@var{Z} = @var{transform}.@qcode{b} * @var{Y} *
## @var{transform}.@qcode{T} + @var{transform}.@qcode{c}}
##
## @qcode{procrustes} can take two optional parameters as Name-Value pairs.
##
## @code{[@dots{}] = procrustes (@dots{}, @qcode{"Scaling"}, @qcode{false})}
## computes a transformation that does not include scaling, that is
## @var{transform}.@qcode{b} = 1.  Setting @qcode{"Scaling"} to @qcode{true}
## includes a scaling component, which is the default.
##
## @code{[@dots{}] = procrustes (@dots{}, @qcode{"Reflection"}, @qcode{false})}
## computes a transformation that does not include a reflection component, that
## is @var{transform}.@qcode{T} = 1.  Setting @qcode{"Reflection"} to
## @qcode{true} forces the solution to include a reflection component in the
## computed transformation, that is @var{transform}.@qcode{T} = -1.
##
## @code{[@dots{}] = procrustes (@dots{}, @qcode{"Reflection"}, @qcode{"best"})}
## computes the best fit procrustes solution, which may or may not include a
## reflection component, which is the default.
##
## @seealso{cmdscale}
## @end deftypefn

function [d, Z, transform] = procrustes (X, Y, varargin)

  if (nargin < 1)
    error ("procrustes: contingency table is missing.");
  endif

  if (nargin > 5)
    error ("procrustes: too many input parameters.");
  endif

  ## Check X and Y for appropriate input
  if (isempty (X) || ! ismatrix (X) || ndims (X) != 2 || ...
      isempty (Y) || ! ismatrix (Y) || ndims (Y) != 2)
    error ("procrustes: X and Y must be 2-dimensional matrices.");
  endif
  if (any (isnan (X(:))) || any (isinf (X(:))) || iscomplex (X) || ...
      any (isnan (Y(:))) || any (isinf (Y(:))) || iscomplex (Y))
    error ("procrustes: values in X and Y must be real.");
  endif
  [Xp, Xd] = size (X);
  [Yp, Yd] = size (Y);
  if (Yp != Xp)
    error ("procrustes: X and Y must have equal number of rows.");
  elseif (Yd > Xd)
    error ("procrustes: X must have at least as many columns as Y.");
  endif

  ## Add defaults and parse optional arguments
  scaling = true;
  reflection = "best";
  if (nargin > 2)
    params = numel (varargin);
    if ((params / 2) != fix (params / 2))
      error ("procrustes: optional arguments must be in Name-Value pairs.")
    endif
    for idx = 1:2:params
      name = varargin{idx};
      value = varargin{idx+1};
      switch (lower (name))
        case "scaling"
          scaling = value;
          if (! (isscalar (scaling) && islogical (scaling)))
            error ("procrustes: invalid value for scaling.");
          endif
        case "reflection"
          reflection = value;
          if (! (strcmpi (reflection, "best") || islogical (reflection)))
            error ("procrustes: invalid value for reflection.");
          endif
        otherwise
          error ("procrustes: invalid name for optional arguments.");
      endswitch
    endfor
  endif

  ## Center at the origin.
  Xmu = mean (X, 1);
  Ymu = mean (Y, 1);
  X_0 = X - repmat (Xmu, Xp, 1);
  Y_0 = Y - repmat (Ymu, Xp, 1);

  ## Get centroid size and check for X or Y having identical points
  Xsumsq = sum (X_0 .^ 2, 1);
  Ysumsq = sum (Y_0 .^ 2, 1);
  constX = all (Xsumsq <= abs (eps (class (X)) * Xp * Xmu) .^ 2);
  constY = all (Ysumsq <= abs (eps (class (X)) * Xp * Ymu) .^ 2);
  Xsumsq = sum (Xsumsq);
  Ysumsq = sum (Ysumsq);

  if (! constX && ! constY)
    ## Scale to "centered" Frobenius norm.
    normX = sqrt (Xsumsq);
    normY = sqrt (Ysumsq);
    X_0 = X_0 / normX;
    Y_0 = Y_0 / normY;

    ## Fix dimension space (if necessary)
    if (Yd < Xd)
        Y_0 = [Y_0 zeros(Xp, Xd-Yd)];
    end

    ## Find optimal rotation matrix of Y
    A = X_0' * Y_0;
    [U, S, V] = svd (A);
    T = V * U';
    ## Handle reflection only if 'true' or 'false' was given
    if (! strcmpi (reflection, "best"))
      is_reflection = (det(T) < 0);
      ## Force a reflection if data and reflection option disagree
      if (reflection != is_reflection)
        V(:,end) = -V(:,end);
        S(end,end) = -S(end,end);
        T = V * U';
      endif
    endif

    ## Apply scaling (if requested)
    traceTA = sum (diag (S));
    if (scaling)
      b = traceTA * normX / normY;
      d = 1 - traceTA .^ 2;
      if (nargout > 1)
        Z = normX * traceTA * Y_0 * T + repmat (Xmu, Xp, 1);
      endif
    else
      b = 1;
      d = 1 + Ysumsq / Xsumsq - 2 * traceTA * normY / normX;
      if (nargout > 1)
        Z = normY * Y_0 * T + repmat (Xmu, Xp, 1);
      endif
    endif

    ## 3rd output argument
    if (nargout > 2)
      if (Yd < Xd)
        T = T(1:Yd,:);
      endif
      c = Xmu - b * Ymu * T;
      transform = struct ("T", T, "b", b, "c", repmat (c, Xp, 1));
    end

  ## Special cases
  elseif constX   # Identical points in X
    d = 0;
    Z = repmat (Xmu, Xp, 1);
    T = eye (Yd, Xd);
    transform = struct ("T", T, "b", 0, "c", Z);
  else            # Identical points in Y
    d = 1;
    Z = repmat (Xmu, Xp, 1);
    T = eye (Yd, Xd);
    transform = struct ("T", T, "b", 0, "c", Z);
  endif
endfunction

%!demo
%! ## Create some random points in two dimensions
%! n = 10;
%! randn ("seed", 1);
%! X = normrnd (0, 1, [n, 2]);
%!
%! ## Those same points, rotated, scaled, translated, plus some noise
%! S = [0.5, -sqrt(3)/2; sqrt(3)/2, 0.5]; # rotate 60 degrees
%! Y = normrnd (0.5*X*S + 2, 0.05, n, 2);
%!
%! ## Conform Y to X, plot original X and Y, and transformed Y
%! [d, Z] = procrustes (X, Y);
%! plot (X(:,1), X(:,2), "rx", Y(:,1), Y(:,2), "b.", Z(:,1), Z(:,2), "bx");

%!demo
%! ## Find Procrustes distance and plot superimposed shape
%!
%! X = [40 88; 51 88; 35 78; 36 75; 39 72; 44 71; 48 71; 52 74; 55 77];
%! Y = [36 43; 48 42; 31 26; 33 28; 37 30; 40 31; 45 30; 48 28; 51 24];
%! plot (X(:,1),X(:,2),"x");
%! hold on
%! plot (Y(:,1),Y(:,2),"o");
%! xlim ([0 100]);
%! ylim ([0 100]);
%! legend ("Target shape (X)", "Source shape (Y)");
%! [d, Z] = procrustes (X, Y)
%! plot (Z(:,1), Z(:,2), "s");
%! legend ("Target shape (X)", "Source shape (Y)", "Transformed shape (Z)");
%! hold off

%!demo
%! ## Apply Procrustes transformation to larger set of points
%!
%! ## Create matrices with landmark points for two triangles
%! X = [5, 0; 5, 5; 8, 5];   # target
%! Y = [0, 0; 1, 0; 1, 1];   # source
%!
%! ## Create a matrix with more points on the source triangle
%! Y_mp = [linspace(Y(1,1),Y(2,1),10)', linspace(Y(1,2),Y(2,2),10)'; ...
%!         linspace(Y(2,1),Y(3,1),10)', linspace(Y(2,2),Y(3,2),10)'; ...
%!         linspace(Y(3,1),Y(1,1),10)', linspace(Y(3,2),Y(1,2),10)'];
%!
%! ## Plot both shapes, including the larger set of points for the source shape
%! plot ([X(:,1); X(1,1)], [X(:,2); X(1,2)], "bx-");
%! hold on
%! plot ([Y(:,1); Y(1,1)], [Y(:,2); Y(1,2)], "ro-", "MarkerFaceColor", "r");
%! plot (Y_mp(:,1), Y_mp(:,2), "ro");
%! xlim ([-1 10]);
%! ylim ([-1 6]);
%! legend ("Target shape (X)", "Source shape (Y)", ...
%!         "More points on Y", "Location", "northwest");
%! hold off
%!
%! ## Obtain the Procrustes transformation
%! [d, Z, transform] = procrustes (X, Y)
%!
%! ## Use the Procrustes transformation to superimpose the more points (Y_mp)
%! ## on the source shape onto the target shape, and then visualize the results.
%! Z_mp = transform.b * Y_mp * transform.T + transform.c(1,:);
%! figure
%! plot ([X(:,1); X(1,1)], [X(:,2); X(1,2)], "bx-");
%! hold on
%! plot ([Y(:,1); Y(1,1)], [Y(:,2); Y(1,2)], "ro-", "MarkerFaceColor", "r");
%! plot (Y_mp(:,1), Y_mp(:,2), "ro");
%! xlim ([-1 10]);
%! ylim ([-1 6]);
%! plot ([Z(:,1); Z(1,1)],[Z(:,2); Z(1,2)],"ks-","MarkerFaceColor","k");
%! plot (Z_mp(:,1),Z_mp(:,2),"ks");
%! legend ("Target shape (X)", "Source shape (Y)", ...
%!         "More points on Y", "Transformed source shape (Z)", ...
%!         "Transformed additional points", "Location", "northwest");
%! hold off

%!demo
%! ## Compare shapes without reflection
%!
%! T = [33, 93; 33, 87; 33, 80; 31, 72; 32, 65; 32, 58; 30, 72; ...
%!      28, 72; 25, 69; 22, 64; 23, 59; 26, 57; 30, 57];
%! S = [48, 83; 48, 77; 48, 70; 48, 65; 49, 59; 49, 56; 50, 66; ...
%!      52, 66; 56, 65; 58, 61; 57, 57; 54, 56; 51, 55];
%! plot (T(:,1), T(:,2), "x-");
%! hold on
%! plot (S(:,1), S(:,2), "o-");
%! legend ("Target shape (d)", "Source shape (b)");
%! hold off
%! d_false = procrustes (T, S, "reflection", false);
%! printf ("Procrustes distance without reflection: %f\n", d_false);
%! d_true = procrustes (T, S, "reflection", true);
%! printf ("Procrustes distance with reflection: %f\n", d_true);
%! d_best = procrustes (T, S, "reflection", "best");
%! printf ("Procrustes distance with best fit: %f\n", d_true);

## Test input validation
%!error procrustes ();
%!error procrustes (1, 2, 3, 4, 5, 6);
%!error<procrustes: X and Y must be 2-dimensional matrices.> ...
%! procrustes (ones (2, 2, 2), ones (2, 2, 2));
%!error<procrustes: values in X and Y must be real.> ...
%! procrustes ([1, 2; -3, 4; 2, 3], [1, 2; -3, 4; 2, 3+i]);
%!error<procrustes: values in X and Y must be real.> ...
%! procrustes ([1, 2; -3, 4; 2, 3], [1, 2; -3, 4; 2, NaN]);
%!error<procrustes: values in X and Y must be real.> ...
%! procrustes ([1, 2; -3, 4; 2, 3], [1, 2; -3, 4; 2, Inf]);
%!error<procrustes: X and Y must have equal number of rows.> ...
%! procrustes (ones (10 ,3), ones (11, 3));
%!error<procrustes: X must have at least as many columns as Y.> ...
%! procrustes (ones (10 ,3), ones (10, 4));
%!error<procrustes: optional arguments must be in Name-Value pairs.> ...
%! procrustes (ones (10 ,3), ones (10, 3), "reflection");
%!error<procrustes: optional arguments must be in Name-Value pairs.> ...
%! procrustes (ones (10 ,3), ones (10, 3), true);
%!error<procrustes: invalid value for scaling.> ...
%! procrustes (ones (10 ,3), ones (10, 3), "scaling", 0);
%!error<procrustes: invalid value for scaling.> ...
%! procrustes (ones (10 ,3), ones (10, 3), "scaling", [true true]);
%!error<procrustes: invalid value for reflection.> ...
%! procrustes (ones (10 ,3), ones (10, 3), "reflection", 1);
%!error<procrustes: invalid value for reflection.> ...
%! procrustes (ones (10 ,3), ones (10, 3), "reflection", "some");
%!error<procrustes: invalid name for optional arguments.> ...
%! procrustes (ones (10 ,3), ones (10, 3), "param1", "some");
