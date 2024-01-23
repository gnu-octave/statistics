## Copyright (C) 2023 Arun Giridhar
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {} einstein ()
## @deftypefnx {statistics} {@var{tiles} =} einstein (@var{a}, @var{b})
## @deftypefnx {statistics} {[@var{tiles}, @var{rhat}] =} einstein (@var{a}, @var{b})
## @deftypefnx {statistics} {[@var{tiles}, @var{rhat}, @var{that}] =} einstein (@var{a}, @var{b})
## @deftypefnx {statistics} {[@var{tiles}, @var{rhat}, @var{that}, @var{shat}] =} einstein (@var{a}, @var{b})
## @deftypefnx {statistics} {[@var{tiles}, @var{rhat}, @var{that}, @var{shat}, @var{phat}] =} einstein (@var{a}, @var{b})
## @deftypefnx {statistics} {[@var{tiles}, @var{rhat}, @var{that}, @var{shat}, @var{phat}, @var{fhat}] =} einstein (@var{a}, @var{b})
##
## Plots the tiling of the basic clusters of einstein tiles.
##
## Scalars @var{a} and @var{b} define the shape of the einstein tile.
## See Smith et al (2023) for details: @url{https://arxiv.org/abs/2303.10798}
##
## @itemize
## @item @var{tiles} is a structure containing the coordinates of the einstein
## tiles that are tiled on the plot.  Each field contains the tile coordinates
## of the corresponding clusters.
## @itemize
## @item @var{tiles}@qcode{.rhat} contains the reflected einstein tiles
## @item @var{tiles}@qcode{.that} contains the three-hat shells
## @item @var{tiles}@qcode{.shat} contains the single-hat clusters
## @item @var{tiles}@qcode{.phat} contains the paired-hat clusters
## @item @var{tiles}@qcode{.fhat} contains the fylfot clusters
## @end itemize
##
## @item @var{rhat} contains the coordinates of the first reflected tile
## @item @var{that} contains the coordinates of the first three-hat shell
## @item @var{shat} contains the coordinates of the first single-hat cluster
## @item @var{phat} contains the coordinates of the first paired-hat cluster
## @item @var{fhat} contains the coordinates of the first fylfot cluster
## @end itemize
##
## @end deftypefn

function [varargout] = einstein (a, b, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    print_usage;
  endif

  ## Check A and B for valid type and range
  if (! (isscalar (a) && isscalar (b) && isnumeric (a) && isnumeric (b)...
      && isreal (a) && isreal (b)))
    error ("einstein: A and B must real scalars.");
  end
  if (a <= 0 || a >= 1 || b <= 0 || b >= 1)
    error ("einstein: A and B must be within the open interval (0,1).");
  endif

  ## Get initial hat points
  single_hat = getpoly (a, b) * rotz (-150)([1,2],[1,2]);

  ## Make a reflected-hat
  reflecthat = [-1, 1] .* single_hat;

  ## Make a three-hat shell cluster
  three_hat1 = rotatehat (single_hat, -60);
  three_hat1 = three_hat1 + (reflecthat(1,:) - three_hat1(9,:));
  three_hat2 = three_hat1 + (reflecthat(6,:) - three_hat1(11,:));
  three_hat3 = rotatehat (three_hat1, 120);
  three_hat3 = three_hat3 + (reflecthat(5,:) - three_hat1(9,:));
  three_hat = [three_hat1, three_hat2, three_hat3];

  ## Translate another four-tile cluster
  translate = three_hat(5,[3,4]) - three_hat(12,[1,2]);
  all.rhat = translatehat (reflecthat, translate);
  all.that = translatehat (three_hat, translate);
  all.rhat = [reflecthat, all.rhat];
  all.that = [three_hat, all.that];

  ## Rotate and translate another four-tile cluster
  tmp_3hat = rotatehat (three_hat, 120);
  tmp_rhat = rotatehat (reflecthat, 120);
  translate = three_hat(12,[5,6]) - tmp_3hat(5,[3,4]);
  tmp_3hat = translatehat (tmp_3hat, translate);
  tmp_rhat = translatehat (tmp_rhat, translate);
  all.that = [all.that, tmp_3hat];
  all.rhat = [all.rhat, tmp_rhat];

  ## Plot four-tile clusters
  patch (all.that(:,[1:2:end]), all.that(:,[2:2:end]), "LineWidth", 2, ...
         "FaceColor", "c", "EdgeColor","k");
  patch (all.rhat(:,[1:2:end]), all.rhat(:,[2:2:end]), "LineWidth", 2, ...
         "FaceColor", "b", "EdgeColor","k");

  title (sprintf ("a = %4.2f  b = %4.2f", a, b), "FontSize",30)

  ## Make a single-hat cluster
  translate = three_hat(10,[3,4]) - single_hat(2,:);
  singlehat = translatehat (single_hat, translate);
  all.shat = singlehat;

  ## Plot single-tile cluster
  patch (all.shat(:,1), all.shat(:,2), "LineWidth", 2, ...
         "FaceColor", [0.95, 0.95, 0.95], "EdgeColor","k");
  axis ("equal")

  ## Make a paired-hat cluster
  paired_hat1 = all.shat + (all.rhat(9,[5,6]) - all.shat(5,:));
  paired_hat2 = three_hat(:,[3,4]) + (three_hat(1,[5,6]) - three_hat(3,[3,4]));
  paired_hat = [paired_hat1, paired_hat2];

  ## Rotate and translate another two paired-hat clusters
  tmp_phat = rotatehat (paired_hat, -120);
  translate = three_hat(4,[3,4]) - tmp_phat(4,[3,4]);
  tmp_phat = translatehat (tmp_phat, translate);
  all.phat = [paired_hat, tmp_phat];
  tmp_phat = rotatehat (paired_hat, -60);
  translate = all.that(10,[11,12]) - tmp_phat(11,[3,4]);
  tmp_phat = translatehat (tmp_phat, translate);
  all.phat = [all.phat, tmp_phat];

  ## Plot paired-tiles clusters
  patch (all.phat(:,[1:2:end]), all.phat(:,[2:2:end]), "LineWidth", 2, ...
         "FaceColor", [0.9, 0.9, 0.9], "EdgeColor","k");

  ## Make a fylfot cluster
  fylfot_hat = paired_hat;
  tmp_fhat = rotatehat (paired_hat, -120);
  translate = fylfot_hat(1,[3,4]) - tmp_fhat(3,[3,4]);
  tmp_fhat = translatehat (tmp_fhat, translate);
  fylfot_hat = [fylfot_hat, tmp_fhat];
  tmp_fhat = rotatehat (paired_hat, 120);
  translate = fylfot_hat(3,[3,4]) - tmp_fhat(1,[3,4]);
  tmp_fhat = translatehat (tmp_fhat, translate);
  fylfot_hat = [fylfot_hat, tmp_fhat];
  translate = all.that(13,[13,14]) - fylfot_hat(13,[7,8]);
  fylfot_hat = translatehat (fylfot_hat, translate);

  ## Translate another two fylfot clusters
  translate = all.that(4,[9,10]) - fylfot_hat(13,[3,4]);
  tmp_fhat = translatehat (fylfot_hat, translate);
  all.fhat = [fylfot_hat, tmp_fhat];
  translate = all.that(13,[1,2]) - fylfot_hat(4,[3,4]);
  tmp_fhat = translatehat (fylfot_hat, translate);
  all.fhat = [all.fhat, tmp_fhat];

  ## Plot fylfot clusters
  patch (all.fhat(:,[1:2:end]), all.fhat(:,[2:2:end]), "LineWidth", 2, ...
         "FaceColor", "r", "EdgeColor","k");

  if (nargout > 0)
    varargout{1} = all;
  endif
  if (nargout > 1)
    varargout{2} = reflecthat;
  endif
  if (nargout > 2)
    varargout{3} = three_hat;
  endif
  if (nargout > 3)
    varargout{4} = singlehat;
  endif
  if (nargout > 4)
    varargout{5} = paired_hat;
  endif
  if (nargout > 5)
    varargout{6} = fylfot_hat;
  endif

endfunction

## Rotates a cluster of hats.
function newhat = rotatehat (hat, degrees)
  rotM = rotz (degrees)([1,2],[1,2]);
  newhat = zeros (size (hat));
  for i=1:2:columns (hat)
    newhat(:,[i,i+1]) = hat(:,[i,i+1]) * rotM;
  endfor
endfunction

## Translates a cluster of hats.
function newhat = translatehat (hat, dist)
  nhats = columns (hat) / 2;
  newhat = hat + repmat (dist, 1, nhats);
endfunction

## Returns a unit vector given a direction angle in degrees
function ret = u (t)
  persistent angles = (0:30:330);
  persistent tbl = [cosd(angles); sind(angles)]';
  ret = tbl(angles == mod (t, 360), :);
endfunction

## Returns the hat polygon
function single_hat = getpoly (a, b)
  single_hat = zeros (14, 2);
  pos = 0;
  single_hat(++pos, :) = [0 0];
  t = 270; single_hat(pos+1, :) = single_hat(pos++, :) + a * u(t);
  t += 60; single_hat(pos+1, :) = single_hat(pos++, :) + a * u(t);
  t -= 90; single_hat(pos+1, :) = single_hat(pos++, :) + b * u(t);
  t += 60; single_hat(pos+1, :) = single_hat(pos++, :) + b * u(t);
  t += 90; single_hat(pos+1, :) = single_hat(pos++, :) + a * u(t);
  t -= 60; single_hat(pos+1, :) = single_hat(pos++, :) + a * u(t);
  t += 90; single_hat(pos+1, :) = single_hat(pos++, :) + b * u(t);
  t -= 60; single_hat(pos+1, :) = single_hat(pos++, :) + b * u(t);
  t += 90; single_hat(pos+1, :) = single_hat(pos++, :) + a * u(t);
  t += 60; single_hat(pos+1, :) = single_hat(pos++, :) + a * u(t) * 2;
  t += 60; single_hat(pos+1, :) = single_hat(pos++, :) + a * u(t);
  t -= 90; single_hat(pos+1, :) = single_hat(pos++, :) + b * u(t);
  t += 60; single_hat(pos+1, :) = single_hat(pos++, :) + b * u(t);
  ## t(14, :) == t(1, :) to within numerical roundoff
endfunction

%!demo
%! einstein (0.4, 0.6)

%!demo
%! einstein (0.2, 0.5)

%!demo
%! einstein (0.6, 0.1)

## Test plotting
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   tiles = einstein (0.4, 0.6);
%!   assert (isstruct (tiles), true);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Test input validation
%!error<Invalid call to einstein.  Correct usage is> einstein
%!error<Invalid call to einstein.  Correct usage is> einstein (0.5)
%!error<einstein: A and B must be within the open interval> einstein (0, 0.9)
%!error<einstein: A and B must be within the open interval> einstein (0.4, 1)
%!error<einstein: A and B must be within the open interval> einstein (-0.4, 1)
