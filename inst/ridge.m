## Copyright (C) 2023, Mohammed Azmat Khan <azmat.dev0@gmail.com>
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


%funciton to implement ridge regression



function b = ridge (y,X,k,flag)

  if  (nargin > 4 || nargin < 3)
    print_usage();
  endif

  if (! ismatrix (y))
    error ("ridge: y must be a numeric matrix");
  endif
  if (! ismatrix (X))
    error ("regress: X must be a numeric matrix");
  endif

  if (nargin < 4 || isempty(flag))
    unscale = false;
  elseif ( flag == 1)
    unscale = false;
  elseif (flag == 0)
    unscale = true;
  endif

  if (columns (y) != 1)
    error ("regress: y must be a column vector");
  endif

  if (rows (y) != rows (X))
    error ("regress: y and X must contain the same number of rows");
  endif

  % remove any missing values
  notnans = ! logical (sum (isnan ([y X]), 2));
  y = y(notnans);
  X = X(notnans,:);

  % normalising X to mean zero and Std deviation ones
  m = mean(X);
  stdx = std(X,0,1);

  z = (X - m) ./ stdx;

  % add pseudo observations
  Z_pseudo = [ X; sqrt(k(1)) * eye( rows(X))];
  Y_pseudo = [ y; zeros( rows(X), 1)];

  % coefficient
  b = Z_pseudo \ Y_pseudo;

  ks = numel(k);

  % compute the coefficient estimates for additional ridge parameters.
  if (ks >= 2)
    % adding a multiple of the identity matrix to the last p rows.

    % b is set to 0 for the current ridge parameter value
   b(ks) = 0;

   for i=2:ks
     Z_pseudo ( end - p + 1 : end) = sqrt( k(i) * eye( rows(X)));
     b(:,j) = Z_pseudo \ Y_pseudo;

   endfor
  endif

   % changing back the scale
   if( unscale)
   b = b ./ repmat( stdx', 1, ks);
   b = [mean(y) - m * b; b];
   endif

endfunction


% write tests and demos
