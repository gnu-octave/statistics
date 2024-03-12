## Copyright (c) 2012 Juan Pablo Carbajal <carbajal@ifi.uzh.ch>
## Copyright (C) 2021 Stefano Guidoni <ilguido@users.sf.net>
## Copyright (C) 2022-2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {} dendrogram (@var{tree})
## @deftypefnx {statistics} {} dendrogram (@var{tree}, @var{p})
## @deftypefnx {statistics} {} dendrogram (@var{tree}, @var{prop}, @var{val})
## @deftypefnx {statistics} {} dendrogram (@var{tree}, @var{p}, @var{prop}, @var{val} )
## @deftypefnx {statistics} {@var{h} =} dendrogram (@dots{})
## @deftypefnx {statistics} {[@var{h}, @var{t}, @var{perm}] =} dendrogram (@dots{})
##
## Plot a dendrogram of a hierarchical binary cluster tree.
##
## Given @var{tree}, a hierarchical binary cluster tree as the output of
## @code{linkage}, plot a dendrogram of the tree. The number of leaves shown by
## the dendrogram plot is limited to @var{p}. The default value for @var{p} is
## 30. Set @var{p} to 0 to plot all leaves.
##
## The optional outputs are @var{h}, @var{t} and @var{perm}:
## @itemize @bullet
## @item @var{h} is a handle to the lines of the plot.
##
## @item @var{t} is the vector with the numbers assigned to each leaf.
## Each element of @var{t} is a leaf of @var{tree} and its value is the number
## shown in the plot.
## When the dendrogram plot is collapsed, that is when the number of shown
## leaves @var{p} is inferior to the total number of leaves, a single leaf of
## the plot can represent more than one leaf of @var{tree}: in that case
## multiple elements of @var{t} share the same value, that is the same leaf of
## the plot.
## When the dendrogram plot is not collapsed, each leaf of the plot is the leaf
## of @var{tree} with the same number.
##
## @item @var{perm} is the vector list of the leaves as ordered as in the plot.
## @end itemize
##
## Additional input properties can be specified by pairs of properties and
## values. Known properties are:
## @itemize @bullet
## @item @qcode{"Reorder"}
## Reorder the leaves of the dendrogram plot using a numerical vector of size n,
## the number of leaves. When @var{p} is smaller than @var{n}, the reordering
## cannot break the @var{p} groups of leaves.
##
## @item @qcode{"Orientation"}
## Change the orientation of the plot. Available values: @qcode{top} (default),
## @qcode{bottom}, @qcode{left}, @qcode{right}.
##
## @item @qcode{"CheckCrossing"}
## Check if the lines of a reordered dendrogram cross each other. Available
## values: @qcode{true} (default), @qcode{false}.
##
## @item @qcode{"ColorThreshold"}
## Not implemented.
##
## @item @qcode{"Labels"}
## Use a char, string or cellstr array of size @var{n} to set the label for each
## leaf; the label is dispayed only for nodes with just one leaf.
## @end itemize
##
## @seealso{cluster, clusterdata, cophenet, inconsistent, linkage, pdist}
## @end deftypefn

function [H, T, perm] = dendrogram (tree, varargin)

  [m, d] = size (tree);
  if ((d != 3) || (! isnumeric (tree)) || (! (max (tree(end, 1:2)) == m * 2)))
    error (strcat (["dendrogram: tree must be a matrix as generated"], ...
                   [" by the linkage function."]));
  endif

  pair_index = 1;

  ## Node count
  n = m + 1;

  ## Add default values
  P = 30;
  vReorder = [];
  csLabels = {};
  checkCrossing = 1;
  orientation = "top";
  if (nargin > 1)
    if (isnumeric (varargin{1}) && isscalar (varargin{1}))
      ## dendrogram (tree, P)
      P = varargin{1};
      pair_index++;
    endif

    ## dendrogram (..., Name, Value)
    while (pair_index < (nargin - 1))
      switch (lower (varargin{pair_index}))
        case "reorder"
          if (isvector (varargin{pair_index + 1}) &&
              isnumeric (varargin{pair_index + 1}) &&
              length (varargin{pair_index + 1}) == n )
            vReorder = varargin{pair_index + 1};
          else
            error (strcat (["dendrogram: 'reorder' must be a numeric"], ...
                           [" vector of size n, the number of leaves."]));
          endif
        case "checkcrossing"
          if (ischar (varargin{pair_index + 1}))
            switch (lower (varargin{pair_index + 1}))
              case "true"
                checkCrossing = 1;
              case "false"
                checkCrossing = 0;
              otherwise
                error ("dendrogram: unknown value '%s' for CheckCrossing.", ...
                       varargin{pair_index + 1});
            endswitch
          else
            error (strcat (["dendrogram: 'CheckCrossing' must be"], ...
                           [" either 'true' or 'false'."]));
          endif
        case "colorthreshold"
          warning ("dendrogram: property '%s' not implemented.",...
            varargin{pair_index});
        case "orientation"
          orientation = varargin{pair_index + 1}; # validity check below
        case "labels"
          if (ischar (varargin{pair_index + 1}) &&
              (isvector (varargin{pair_index + 1}) &&
               length (varargin{pair_index + 1}) == n) ||
              (ismatrix (varargin{pair_index + 1}) &&
               rows (varargin{pair_index + 1}) == n))
            csLabels = cellstr (varargin{pair_index + 1});
          elseif (iscellstr (varargin{pair_index + 1}) &&
                  length (varargin{pair_index + 1}) == n)
            csLabels = varargin{pair_index + 1};
          else
            error (strcat (["dendrogram: labels must be a char or"], ...
                           [" string or cellstr array of size n."]));
          endif
        otherwise
          error ("dendrogram: unknown property '%s'.", varargin{pair_index});
      endswitch

      pair_index += 2;
    endwhile
  endif

  ## MATLAB compatibility:
  ## P <= 0 to plot all leaves
  if (P < 1)
    P = n;
  endif

  if (n > P)
    level_0 = tree((n - P), 3);
  else
    P = n;
    level_0 = 0;
  endif

  vLeafPosition = zeros((n + m), 1);
  T = (1:n)';
  nodecnt = 1;

  ## main
  dendrogram_recursive (m, 0);

  ## T reordering
  ## MATLAB compatibility: when n > P, each node group is renamed with a number
  ## between 1 and P, according to the smallest node index of each group;
  ## the group with the node 1 is always group 1, while group 2 is the group
  ## with the smallest node index outside of group 1, and group 3 is the group
  ## with the smallest node index outside of groups 1 and 2...
  newT = 1 : (length (T));
  if (n > P)
    uniqueT = unique (T);
    minT = zeros (size (uniqueT));
    counter = 1;

    for i = 1:length (uniqueT)  # it should be exactly equal to P
      idcs = find (T == uniqueT(i));
      minT(i) = min (idcs);
    endfor

    minT = minT(find (minT > 0)); # to prevent a strange bug
    [minT, minTidcs] = sort (minT);
    uniqueT = uniqueT(minTidcs);
    for i = 1:length (uniqueT)
      idcs = find (T == uniqueT(i));
      newT(idcs) = counter++;
    endfor
  endif

  ## leaf reordering
  if (! isempty(vReorder))
    if (P < n)
      checkT = newT(vReorder(:));
      for i = 1 : P
        idcs = find (checkT == i);
        if (length (idcs) > 1)
          if (max (idcs) - min (idcs) >= length (idcs))
            error (strcat (["dendrogram: invalid reordering that"], ...
                           [" redefines the 'P' groups of leaves"]));
          endif
        endif
      endfor
      checkT = unique (checkT, "stable");
      vNewLeafPosition = zeros (n, 1);
      uT = unique (T, "stable");
      for i = 1:P
        vNewLeafPosition(uT(checkT(i))) = i;
      endfor
      vLeafPosition = vNewLeafPosition;
    else
      for i = 1:length (vReorder)
        vLeafPosition(vReorder(i)) = i;
      endfor
    endif
  endif

  ## figure
  x = [];

  ## ticks and tricks
  xticks = 1:P;
  perm = zeros (P, 1);
  for i = 1 : length (vLeafPosition)
    if (vLeafPosition(i) != 0)
      idcs = find (T == i);
      perm(vLeafPosition(i)) = newT(idcs(1));
    endif
  endfor
  T = newT; # this should be unnecessary for n <= P

  ## lines
  for i = (n - P + 1):m
    vLeafPosition(n + i) = mean (vLeafPosition(tree(i, 1:2), 1));
    x(end + 1,1:4) = [vLeafPosition(tree(i, 1:2))' tree(i, [3 3])];
    for j = 1 : 2
      x0 = 0;
      if (tree(i,j) > (2 * n - P))
        x0 = tree(tree(i, j) - n, 3);
      endif
      x(end + 1, 1:4) = [vLeafPosition(tree(i, [j j]))' x0 tree(i, 3)];
    endfor
  endfor

  ## plot stuff
  if (strcmp (orientation, "top"))
    H = line (x(:, 1:2)', x(:, 3:4)', "color", "blue");
    set (gca, "xticklabel", perm, "xtick", xticks);
  elseif (strcmp (orientation, "bottom"))
    H = line (x(:, 1:2)', x(:, 3:4)', "color", "blue");
    set (gca, "xticklabel", perm, "xtick", xticks, "xaxislocation", "top");
    axis ("ij");
  elseif (strcmp (orientation, "left"))
    H = line (x(:, 3:4)', x(:, 1:2)', "color", "blue");
    set (gca, "yticklabel", perm, "ytick", xticks, "xdir", "reverse",...
         "yaxislocation", "right");
  elseif (strcmp (orientation, "right"))
    H = line (x(:, 3:4)', x(:, 1:2)', "color", "blue");
    set (gca, "yticklabel", perm, "ytick", xticks);
  else
    close (H);
    error ("dendrogram: invalid orientation '%s'", orientation);
  endif

  ## labels
  if (! isempty (csLabels))
    csCurrent = cellstr (num2str (perm));

    for i = 1:n
      ## when there is just one leaf, use the named label for that leaf
      if (1 == length (find (T == i)))
        csCurrent(find (perm == i)) = csLabels(find (T == i));
      endif
    endfor

    switch (orientation)
      case {"top", "bottom"}
        xticklabels (csCurrent);
      case {"left", "right"}
        yticklabels (csCurrent);
    endswitch
  endif

  ## check crossings
  if (checkCrossing && ! isempty(vReorder))
    for j = 1:rows (x)
      if (x(j, 3) == x(j, 4)) # an horizontal line
        for i = 1:rows (x)
          if (x(i, 1) == x(i, 2) && ... # orthogonal lines
            (x(i, 1) > x(j, 1) && x(i, 1) < x(j, 2)) && ...
            (x(j, 3) > x(i, 3) && x(j, 3) < x(i, 4)))
            warning ("dendrogram: line intersection detected");
          endif
        endfor
      endif
    endfor
  endif

  ## dendrogram_recursive
  function dendrogram_recursive (k, cn)
    if (tree(k, 3) > level_0)
      for j = 1:2
        if (tree(k, j) > n)
          dendrogram_recursive (tree(k, j) - n, 0)
        else
          vLeafPosition(tree(k, j)) = nodecnt++;
          T(tree(k, j)) = tree(k, j);
        endif
      endfor
    else
      for j = 1:2
        if (cn == 0)
          cn = n + k;
          vLeafPosition(cn) = nodecnt++;
        endif

        if (tree(k, j) > n)
          dendrogram_recursive (tree(k, j) - n, cn)
        else
          T(tree(k, j)) = cn;
        endif
      endfor
    endif
  endfunction
endfunction


%!demo
%! ## simple dendrogram
%! y = [4, 5; 2, 6; 3, 7; 8, 9; 1, 10];
%! y(:,3) = 1:5;
%! dendrogram (y);
%! title ("simple dendrogram");

%!demo
%! ## another simple dendrogram
%! v = 2 * rand (30, 1) - 1;
%! d = abs (bsxfun (@minus, v(:, 1), v(:, 1)'));
%! y = linkage (squareform (d, "tovector"));
%! dendrogram (y);
%! title ("another simple dendrogram");

%!demo
%! ## collapsed tree, find all the leaves of node 5
%! X = randn (60, 2);
%! D = pdist (X);
%! y = linkage (D, "average");
%! subplot (2, 1, 1);
%! title ("original tree");
%! dendrogram (y, 0);
%! subplot (2, 1, 2);
%! title ("collapsed tree");
%! [~, t] = dendrogram (y, 20);
%! find(t == 5)

%!demo
%! ## optimal leaf order
%! X = randn (30, 2);
%! D = pdist (X);
%! y = linkage (D, "average");
%! order = optimalleaforder (y, D);
%! subplot (2, 1, 1);
%! title ("original leaf order");
%! dendrogram (y);
%! subplot (2, 1, 2);
%! title ("optimal leaf order");
%! dendrogram (y, "Reorder", order);

%!demo
%! ## horizontal orientation and labels
%! X = randn (8, 2);
%! D = pdist (X);
%! L = ["Snow White"; "Doc"; "Grumpy"; "Happy"; "Sleepy"; "Bashful"; ...
%!      "Sneezy"; "Dopey"];
%! y = linkage (D, "average");
%! dendrogram (y, "Orientation", "left", "Labels", L);
%! title ("horizontal orientation and labels");

## Test plotting
%!shared visibility_setting
%! visibility_setting = get (0, "DefaultFigureVisible");
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   y = [4, 5; 2, 6; 3, 7; 8, 9; 1, 10];
%!   y(:,3) = 1:5;
%!   dendrogram (y);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   y = [4, 5; 2, 6; 3, 7; 8, 9; 1, 10];
%!   y(:,3) = 1:5;
%!   dendrogram (y);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   v = 2 * rand (30, 1) - 1;
%!   d = abs (bsxfun (@minus, v(:, 1), v(:, 1)'));
%!   y = linkage (squareform (d, "tovector"));
%!   dendrogram (y);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   X = randn (30, 2);
%!   D = pdist (X);
%!   y = linkage (D, "average");
%!   order = optimalleaforder (y, D);
%!   subplot (2, 1, 1);
%!   title ("original leaf order");
%!   dendrogram (y);
%!   subplot (2, 1, 2);
%!   title ("optimal leaf order");
%!   dendrogram (y, "Reorder", order);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Test input validation
%!error dendrogram ();
%!error <tree must be .*> dendrogram (ones (2, 2), 1);
%!error <unknown property .*> dendrogram ([1 2 1], 1, "xxx", "xxx");
%!error <reorder.*> dendrogram ([1 2 1], "Reorder", "xxx");
%!error <reorder.*> dendrogram ([1 2 1], "Reorder", [1 2 3 4]);
%! fail ('dendrogram ([1 2 1], "Orientation", "north")', "invalid orientation .*")
