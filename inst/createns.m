## Copyright (C) 2025 Swayam Shah <swayamshah66@gmail.com>
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
## @deftypefn  {Function File} {@var{obj} =} createns (@var{X})
## @deftypefnx {Function File} {@var{obj} =} createns (@var{X}, @var{name}, @var{value}, @dots{})
##
## Create a nearest neighbor searcher object.
##
## @code{@var{obj} = createns (@var{X})} creates a nearest neighbor searcher
## object using the training data @var{X}. By default, it constructs an
## @code{ExhaustiveSearcher} object with the Euclidean distance metric.
##
## @code{@var{obj} = createns (@var{X}, @var{name}, @var{value}, @dots{})}
## allows customization of the searcher type and its properties through
## name-value pairs. The following name-value pair is supported to specify
## the searcher type:
##
## @multitable @columnfractions 0.18 0.02 0.8
## @headitem @var{Name} @tab @tab @var{Value}
##
## @item @qcode{"NSMethod"} @tab @tab Specifies the nearest neighbor search
## method. Possible values are:
##   @itemize @bullet
##   @item @qcode{"exhaustive"}: Creates an @code{ExhaustiveSearcher} object.
##   @item @qcode{"kdtree"}: Creates a @code{KDTreeSearcher} object.
##   @item @qcode{"hnsw"}: Creates an @code{hnswSearcher} object.
##   @end itemize
##   Default is @qcode{"exhaustive"}.
##
## @end multitable
##
## Additional name-value pairs depend on the selected @qcode{"NSMethod"} and
## are passed directly to the constructor of the corresponding class:
##
## @itemize @bullet
## @item For @qcode{"exhaustive"}, see @code{ExhaustiveSearcher} documentation
## for parameters like @qcode{"Distance"}, @qcode{"P"}, @qcode{"Scale"}, and
## @qcode{"Cov"}.
## @item For @qcode{"kdtree"}, see @code{KDTreeSearcher} documentation for
## parameters like @qcode{"Distance"}, @qcode{"P"}, and @qcode{"BucketSize"}.
## @item For @qcode{"hnsw"}, see @code{hnswSearcher} documentation for parameters
## like @qcode{"Distance"}, @qcode{"P"}, @qcode{"Scale"}, @qcode{"Cov"},
## @qcode{"M"}, @qcode{"efConstruction"}, and @qcode{"efSearch"}.
## @end itemize
##
## @strong{Input Arguments:}
## @itemize @bullet
## @item @var{X} - Training data, specified as an @math{NxP} numeric matrix
## where rows represent observations and columns represent features. Must be
## finite and numeric.
## @end itemize
##
## @strong{Output:}
## @itemize @bullet
## @item @var{obj} - A nearest neighbor searcher object of type
## @code{ExhaustiveSearcher}, @code{KDTreeSearcher}, or @code{hnswSearcher},
## depending on the specified @qcode{"NSMethod"}.
## @end itemize
##
## @strong{Examples:}
##
## @example
## ## Create an ExhaustiveSearcher with default parameters
## X = [1, 2; 3, 4; 5, 6];
## obj = createns (X);
##
## ## Create a KDTreeSearcher with Euclidean distance
## obj = createns (X, "NSMethod", "kdtree", "Distance", "euclidean");
##
## ## Create an hnswSearcher with Minkowski distance and custom parameters
## obj = createns (X, "NSMethod", "hnsw", "Distance", "minkowski", "P", 3, "M", 20);
## @end example
##
## @seealso{ExhaustiveSearcher, KDTreeSearcher, hnswSearcher, knnsearch, rangesearch}
## @end deftypefn

function obj = createns (X, varargin)
  ## Input validation
  if (nargin < 1)
    error ("createns: too few input arguments.");
  endif

  if (mod (numel (varargin), 2) != 0)
    error ("createns: Name-Value arguments must be in pairs.");
  endif

  if (! (isnumeric (X) && ismatrix (X) && all (isfinite (X)(:))))
    error ("createns: X must be a finite numeric matrix.");
  endif

  ## Set default NSMethod
  NSMethod = "exhaustive";

  ## Extract 'NSMethod' from varargin and remove it
  names = varargin(1:2:end);
  idx = find (strcmpi (names, "nsmethod"));
  if (! isempty (idx))
    if (length (idx) > 1)
      warning ("createns: multiple 'NSMethod' specified, using the last one.");
    endif
    NSMethod = varargin{2*idx(end)};
    if (! ischar (NSMethod))
      error ("createns: 'NSMethod' must be a string.");
    endif
    NSMethod = lower (NSMethod);
    remove_idx = [];
    for i = idx
      remove_idx = [remove_idx, 2*i-1, 2*i];
    endfor
    varargin(remove_idx) = [];
  endif

  ## Validate NSMethod value
  allowed_methods = {"exhaustive", "kdtree", "hnsw"};
  if (! any (strcmp (NSMethod, allowed_methods)))
    error ("createns: invalid 'NSMethod' value: '%s'.", NSMethod);
  endif

  ## Instantiate the appropriate searcher object
  switch (NSMethod)
    case "exhaustive"
      obj = ExhaustiveSearcher (X, varargin{:});
    case "kdtree"
      obj = KDTreeSearcher (X, varargin{:});
    case "hnsw"
      obj = hnswSearcher (X, varargin{:});
  endswitch

endfunction

## Test Cases

%!test
%! ## Default ExhaustiveSearcher
%! X = [1, 2; 3, 4; 5, 6];
%! obj = createns (X);
%! assert (isa (obj, "ExhaustiveSearcher"));
%! assert (obj.X, X);
%! assert (obj.Distance, "euclidean");

%!test
%! ## KDTreeSearcher with default parameters
%! X = [1, 2; 3, 4; 5, 6];
%! obj = createns (X, "NSMethod", "kdtree");
%! assert (isa (obj, "KDTreeSearcher"));
%! assert (obj.X, X);
%! assert (obj.Distance, "euclidean");

%!test
%! ## hnswSearcher with custom parameters
%! X = [1, 2; 3, 4; 5, 6];
%! obj = createns (X, "NSMethod", "hnsw", "M", 20, "efSearch", 30);
%! assert (isa (obj, "hnswSearcher"));
%! assert (obj.X, X);
%! assert (obj.M, 20);
%! assert (obj.efSearch, 30);

%!test
%! ## ExhaustiveSearcher with custom distance
%! X = [1, 2; 3, 4];
%! obj = createns (X, "NSMethod", "exhaustive", "Distance", "cityblock");
%! assert (isa (obj, "ExhaustiveSearcher"));
%! assert (obj.Distance, "cityblock");

%!error<createns: too few input arguments.>
%! createns ()
%!error<createns: Name-Value arguments must be in pairs.>
%! X = [1, 2; 3, 4]; createns (X, "NSMethod")
%!error<createns: X must be a finite numeric matrix.>
%! createns ([1; Inf; 3])
%!error<createns: 'NSMethod' must be a string.>
%! X = [1, 2; 3, 4]; createns (X, "NSMethod", 1)
%!error<createns: invalid 'NSMethod' value: 'invalid'.>
%! X = [1, 2; 3, 4]; createns (X, "NSMethod", "invalid")