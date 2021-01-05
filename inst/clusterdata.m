## Copyright (C) 2021 Stefano Guidoni <ilguido@users.sf.net>
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
## @deftypefn {Function File} {@var{T} =} clusterdata (@var{X}, @var{cutoff})
## @deftypefnx {Function File} @
##   {@var{T} =} clusterdata (@var{X}, @var{Name}, @var{Value})
##
## Wrapper function for @code{linkage} and @code{cluster}.
##
## If @var{cutoff} is used, then @code{clusterdata} calls @code{linkage} and
## @code{cluster} with default value, using @var{cutoff} as a threshold value 
## for @code{cluster}. If @var{cutoff} is an integer and greater or equal to 2,
## then @var{cutoff} is interpreted as the maximum number of cluster desired
## and the 'MaxClust' option is used for @code{cluster}.
##
## If @var{cutoff} is not used, then @code{clusterdata} expects a list of pair
## arguments. Then you must specify either the 'Cutoff' or 'MaxClust' option
## for @code{cluster}. The method and metric used by @code{linkage}, are
## defined through the 'linkage' and 'distance' arguments.
##
## @end deftypefn
##
## @seealso{cluster,dendrogram,inconsistent,kmeans,linkage,pdist}

## Author: Stefano Guidoni <ilguido@users.sf.net>

function T = clusterdata (X, varargin)
  if (nargin < 2)
    print_usage ();
  else
    linkage_criterion = "single";
    distance_method = "euclidean";
    savememory = "off";
    clustering_method = [];
    criterion = "inconsistent";
    D = 2;
    
    if (isnumeric (varargin{1})) # clusterdata (X, cutoff)
      if (isinteger (varargin{1}) && (varargin{1} >= 2))
        clustering_method = "MaxClust";
      else
        clustering_method = "Cutoff";
      endif
      C = varargin{1};
    else # clusterdata (Name, Value)
      pair_index = 1;
      
      while (pair_index < (nargin - 1))
        switch (lower (varargin{pair_index}))
          case "criterion"
            criterion = varargin{pair_index + 1};
          case "cutoff"
            clustering_method = "Cutoff";
            C = varargin{pair_index + 1};
          case "depth"
            D = varargin{pair_index + 1};
          case "distance"
            distance_method = varargin{pair_index + 1};
          case "linkage"
            linkage_criterion = varargin{pair_index + 1};
          case "maxclust"
            clustering_method = "MaxClust";
            C = varargin{pair_index + 1};
          case "savememory"
            savememory = varargin{pair_index + 1};
          otherwise
            error ("clusterdata: unknown property %s", varargin{pair_index});
        endswitch
        
        pair_index += 2;
      endwhile  
    endif
  endif
  
  if (isempty (clustering_method))
    error ...
      (["clusterdata: you must specify either 'MaxClust' or 'Cutoff' when" ...
        "using name-value arguments"]);
  endif
  
  ## main body
  Z = linkage (X, linkage_criterion, distance_method, "savememory", savememory);
  if (strcmp (lower (clustering_method), "cutoff"))
    T = cluster (Z, clustering_method, C, "Criterion", criterion, "Depth", D);
  else
    T = cluster (Z, clustering_method, C);
  endif  
endfunction

## Test input validation
%!error clusterdata ()
%!error clusterdata (1)
%!error <unknown property .*> clusterdata ([1 1], "Bogus", 1)
%!error <specify .* 'MaxClust' or 'Cutoff' .*> clusterdata ([1 1], "Depth", 1)
## Demonstation
%!demo
%! X = [(randn (10, 2) * 0.25) + 1; (randn (20, 2) * 0.5) - 1];
%! T = clusterdata (X, "linkage", "ward", "MaxClust", 2);
%! scatter (X(:,1), X(:,2), 36, T, "filled");
