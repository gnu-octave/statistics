## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Function File} manovacluster (@var{stats})
## @deftypefnx {Function File} manovacluster (@var{stats}, @var{method})
## @deftypefnx {Function File} @var{h} = manovacluster (@var{stats})
## @deftypefnx {Function File} @var{h} = manovacluster (@var{stats}, @var{method})
##
## Cluster group means using manova1 output.
##
## @code{manovacluster (@var{stats})} draws a dendrogram showing the clustering
## of group means, calculated using the output STATS structure from
## @code{manova1} and applying the single linkage algorithm.  See the
## @code{dendrogram} function for more information about the figure.
##
## @code{manovacluster (@var{stats}, @var{method})} uses the @var{method}
## algorithm in place of single linkage.  The available methods are:
##
## @multitable @columnfractions 0.05 0.2 0.75
## @item @tab "single" @tab --- nearest distance
## @item @tab "complete" @tab --- furthest distance
## @item @tab "average" @tab --- average distance
## @item @tab "centroid" @tab --- center of mass distance
## @item @tab "ward" @tab --- inner squared distance
## @end multitable
##
## @code{@var{h} = manovacluster (@dots{})} returns a vector of line handles.
##
## @seealso{manova1}
## @end deftypefn

function h = manovacluster (stats, method)
  
  ## Check for valid input arguments
  narginchk (1, 2);
  if nargin > 1
    valid_methods = {"single", "complete", "average", "centroid", "ward"};
    if ! any (strcmpi (method, valid_methods))
      error ("manovacluster: invalid method.");
    endif
  else
    method = "single";
  end
  ## Get stats fields and create dendrogram
  dist = stats.gmdist;
  group_names = stats.gnames;
  [a, b] = meshgrid (1:length (dist));
  hh = dendrogram (linkage (dist(a < b)', method), 0);
  ## Fix tick labels on x-axis
  oldlab = get (gca, "XTickLabel");
  maxlen = max (cellfun ("length", group_names));
  newlab = repmat(" ", size (oldlab, 1), maxlen);
  ng = size (group_names, 1);
  for j = 1:size (oldlab, 1)
    k = str2num (oldlab(j,:));
    if (! isempty (k) & k > 0 & k <= ng)
      x = group_names{k,:};
      newlab(j,1:length(x)) = x;
    endif
  endfor
  set(gca, "XtickLabel", newlab);
  ## Return plot handles if requested
  if nargout > 0
    h = hh;
  endif
endfunction

%!demo
%! load carbig
%! X = [MPG Acceleration Weight Displacement];
%! [d, p, stats] = manova1 (X, Origin);
%! manovacluster (stats)