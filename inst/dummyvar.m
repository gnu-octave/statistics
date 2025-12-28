## Copyright (C) 2025 Jayant Chauhan <0001jayant@gmail.com>
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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.
##
## -*- texinfo -*-
## @deftypefn  {statistics} {@var{D} =} dummyvar (@var{g})
##
## Create dummy variables (one-hot encoding) from a grouping variable.
##
## The input @var{g} must be a numeric vector of group indices or a
## @code{categorical} array.  The output @var{D} is a numeric matrix whose
## columns correspond to the distinct groups and whose rows correspond to the
## elements of @var{g}.
##
## For numeric inputs, the number of columns in @var{D} is equal to
## @code{max (@var{g})}, and column @math{k} corresponds to group @math{k}.
##
## For @code{categorical} inputs, the number and order of columns in @var{D}
## correspond to the categories returned by @code{categories (@var{g})}.
## Categories that are defined but not present in @var{g} produce columns of
## zeros.
##
## Elements of @var{g} that are @code{<undefined>} result in rows of
## @code{NaN} values in @var{D}.
##
## If @var{g} is a single-column table, the grouping variable is taken from that
## column.  For example:
##
## @example
## D = dummyvar (T.Group)
## @end example
##
## @seealso{tabulate, grpstats}
## @end deftypefn

function D = dummyvar (g)

  ## Table single-column extraction (non-fatal if 'table' not present)
  try
    if isa(g, "table")
      if (size(g,2) ~= 1)
        error ("dummyvar:InvalidInput", ...
               "dummyvar on a table expects a single-column input or call dummyvar(T.Var).");
      end
      g = g{:,1};
    end
  catch
    ## If table class not available, skip extraction and let later checks handle it.
  end

  ## --- CATEGORICAL branch ---
  if (exist ("categorical", "class") && isa (g, "categorical"))
    ## "Categorical grouping variable must have one column."
    s = size (g);
    if (numel (s) > 2 || s(2) ~= 1)
      error ("Categorical grouping variable must have one column.");
    end

    ## categories is a categorical class method; no standalone function check needed

    cats = categories (g);
    K = numel (cats);
    n = numel (g);

    if (n == 0)
      error ("Categorical grouping variable must have one column.");
      return;
    end

    ## Convert to numeric indices. 
    ## Unknown/undefined map to NaN or 0 depending on API.
    ## <undefined> results in NaN rows in output.

    try
      idx = double (g);    ## maps categories to 1..K, undefined -> NaN
    catch
      ## fallback: construct mapping manually (less efficient)
      ## Use grp2idx-like approach
      [~, ~, idx] = unique (g);
      idx = double (idx);
      ## Unique will not produce NaN for undefined; handle undefined explicitly
      undef_mask = isundefined (g);
      idx (undef_mask) = NaN;
    end

    ## Build matrix: rows with idx==NaN are rows of NaN(1,K)
    rows = (1:n)';
    D = zeros (n, K);

    ## Build using sparse (skip NaNs)
    valid = ~isnan (idx);
    if any (valid)
      S = sparse (rows(valid), idx(valid), 1, n, K);
      D (:) = full (S);
    end

    ## Replace rows where idx is NaN with NaN across all columns
    nan_rows = find (isnan (idx));
    if ~ isempty (nan_rows)
      D (nan_rows, :) = NaN;
    end

    D = double (D);
    return;
  end

  ## --- NUMERIC / LEGACY branch ---
  ## If the input is numeric (vector), create indicator columns for 1..max(g)
  if (isnumeric (g) && isvector (g))
     ## ensure column
    g = g(:);                
    if isempty (g)
      D = zeros (0, 0);
      return;
    end
    ## If g has non-integer values, we implicitly coerces to integer indices
    ## by using unique? Historically dummyvar expects group indices (positive ints).
    ## We'll follow this behavior: use max(g) as number of columns.
    ## Construct column count K = max(g)
    K = max (g);
    if ~ isreal (K) || K < 0
      error ("dummyvar:InvalidInput", "Numeric grouping must produce a positive integer number of groups.");
    end
    K = double (K);
    n = numel (g);
    ## Build sparse/dense
    rows = (1:n)';
    ## keep integer mapping
    idx = round (g);          
    valid = (idx >= 1) & (idx <= K) & ~isnan (idx);
    if any (valid)
      S = sparse (rows(valid), idx(valid), 1, n, K);
      D = full (S);
    else
      D = zeros (n, K);
    end
    D = double (D);
    return;
  end

  ## --- FALLBACK: unsupported input types ---
  error ("dummyvar:UnsupportedType", "dummyvar requires a numeric vector or a categorical array.");
end

## Test dummyvar behavior 

%!test
%! ## numeric grouping vector
%! g = [1;2;1;3;2];
%! D = dummyvar(g);
%! assert(isequal(D, [1 0 0; 0 1 0; 1 0 0; 0 0 1; 0 1 0]));

%!test
%! ## categorical with universe -> columns for each category in same order
%! g = categorical({'a';'b';'a'}, {'a','b','c'});
%! D = dummyvar(g);
%! assert(size(D,2) == numel(categories(g)));
%! assert(all(D(:,1) == [1;0;1]));
%! assert(all(D(:,2) == [0;1;0]));
%! assert(all(D(:,3) == [0;0;0]));

%!test
%! ## categorical with <undefined> -> row of NaNs
%! g = categorical({'a'; ''; 'b'}, {'a','b','c'});
%! D = dummyvar(g);
%! assert(all(isnan(D(2,:))));
%! assert(all(D(1,:) == [1 0 0]));
%! assert(all(D(3,:) == [0 1 0]));

%!test
%! ## empty categorical -> 
%! g = categorical({}, {'a','b'});
%! assert (throws (@() dummyvar(g)));

%!test
%! ## table column input
%! G = categorical({'a'; 'b'; 'a'}, {'a','b','c'});
%! T = table(G, [10;20;30], 'VariableNames', {'G','Val'});
%! D = dummyvar(T.G);
%! assert(size(D,2) == numel(categories(G)));

## Test input validation

%!error<Invalid call to dummyvar.  Correct usage is> dummyvar
%!error<Invalid call to dummyvar.  Correct usage is> dummyvar (1, 2)

%!error<Categorical grouping variable must have one column> ...
%! dummyvar (categorical ([], {'a','b'}))

%!error<Categorical grouping variable must have one column> ...
%! dummyvar (categorical ({'a','b'}, {'a','b'}))   ## row categorical

%!error<dummyvar on a table expects a single-column input> ...
%! dummyvar (table ([1;2], [3;4]))

%!error<dummyvar requires a numeric vector or a categorical array> ...
%! dummyvar (struct ("a", 1))
