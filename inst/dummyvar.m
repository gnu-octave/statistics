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

  if (nargin ~= 1)
    error ("Invalid call to dummyvar.  Correct usage is:\n\n  D = dummyvar (g)");
  end

  ## Table single-column extraction
  if (isa (g, "table"))
    if (size (g, 2) ~= 1)
      error ("dummyvar on a table expects a single-column input");
    end
    try
      g = g{:,1};
    catch
      ## table class exists but indexing failed — let later checks handle it
    end
  end


  ## --- CATEGORICAL branch ---
  if (exist ("categorical", "class") && isa (g, "categorical"))

    if (! isvector (g) || size (g,2) ~= 1)
      error ("Categorical grouping variable must have one column.");
    end

    cats = cellstr (categories (g));
    K = numel (cats);
    n = rows (g);

    if (n == 0)
      error ("Categorical grouping variable must have one column.");
    end

    g_str = cellstr (g(:));
    D = zeros (n, K);

    for i = 1:n
      if (isundefined (g(i)))
        D(i,:) = NaN;
      else
        for k = 1:K
          if (strcmp (g_str{i}, cats{k}))
            D(i,k) = 1;
            break;
          end
        endfor
      end
    endfor

    D = double (D);
    return;
  end

  ## --- NUMERIC branch ---
  if (isnumeric (g) && isvector (g))

    g = g(:);
    if (isempty (g))
      D = zeros (0, 0);
      return;
    end

    K = max (g);
    if (! isreal (K) || K < 0)
      error ("dummyvar:InvalidInput", ...
             "Numeric grouping must produce a positive integer number of groups.");
    end

    n = numel (g);
    D = zeros (n, K);
    idx = round (g);

    for i = 1:n
      if (! isnan (idx(i)) && idx(i) >= 1 && idx(i) <= K)
        D(i, idx(i)) = 1;
      end
    endfor

    D = double (D);
    return;
  end

  error ("dummyvar:UnsupportedType", ...
         "dummyvar requires a numeric vector or a categorical array.");
end

## Test dummyvar behavior 

%!test
%! ## numeric grouping vector
%! g = [1;2;1;3;2];
%! D = dummyvar(g);
%! assert(isequal(D, [1 0 0; 0 1 0; 1 0 0; 0 0 1; 0 1 0]));

%!test
%! g = categorical({'a';'b';'a'}, {'a','b','c'});
%! D = dummyvar(g);
%! cats = categories(g);
%! g_str = cellstr(g);
%!
%! for k = 1:numel(cats)
%!   mask = strcmp(g_str, cats{k});
%!   assert(all(D(mask, k) == 1));
%!   assert(all(D(!mask, k) == 0));
%! endfor


%!test
%! g = categorical({'a'; ''; 'b'}, {'a','b','c'});
%! D = dummyvar(g);
%! assert(all(isnan(D(2,:))));
%! assert(sum(D(1,:) == 1) == 1);
%! assert(sum(D(3,:) == 1) == 1);

%!test
%! G = categorical({'a'; 'b'; 'a'}, {'a','b','c'});
%! T = table(G, [10;20;30], 'VariableNames', {'G','Val'});
%! D = dummyvar(T.G);
%! assert(size(D,2) == numel(categories(G)));

%!test
%! ## table column input
%! G = categorical({'a'; 'b'; 'a'}, {'a','b','c'});
%! T = table(G, [10;20;30], 'VariableNames', {'G','Val'});
%! D = dummyvar(T.G);
%! assert(size(D,2) == numel(categories(G)));

## Test input validation

%!error<Invalid call to dummyvar.  Correct usage is> dummyvar
%!error<function called with too many inputs> dummyvar (1, 2)

%!error<Categorical grouping variable must have one column> ...
%! dummyvar (categorical ([], {'a','b'}))

%!error<Categorical grouping variable must have one column> ...
%! dummyvar (categorical ({'a','b'}, {'a','b'}))   ## row categorical

%!error<dummyvar on a table expects a single-column input> ...
%! dummyvar (table ([1;2], [3;4]))

%!error<dummyvar requires a numeric vector or a categorical array> ...
%! dummyvar (struct ("a", 1))
