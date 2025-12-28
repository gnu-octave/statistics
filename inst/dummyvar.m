function D = dummyvar (g)
  %#ok<*STRETCH>  % for compatibility with different Octave linters

  % Table single-column extraction (non-fatal if 'table' not present)
  try
    if isa(g, "table")
      if (size(g,2) ~= 1)
        error ("dummyvar:InvalidInput", ...
               "dummyvar on a table expects a single-column input or call dummyvar(T.Var).");
      end
      g = g{:,1};
    end
  catch
    % If table class not available, skip extraction and let later checks handle it.
  end

  % --- CATEGORICAL branch ---
  if (exist ("categorical", "class") && isa (g, "categorical"))
    % Match MATLAB error for empty categorical shape: MATLAB complains
    % "Categorical grouping variable must have one column."
    % In MATLAB that arises when categorical is 0x0 or not a single column.
    s = size (g);
    if (numel (s) > 2 || s(2) ~= 1)
      error ("Categorical grouping variable must have one column.");
    end

    % categories and double mapping must exist
    if ~ (exist ("categories", "file") || exist ("categories", "builtin"))
      error ("datatypes:Missing", "datatypes: 'categories' not found on path; load datatypes.");
    end

    cats = categories (g);
    K = numel (cats);
    n = numel (g);

    % For empty (0x1) categorical, MATLAB errors out above; if not thrown,
    % keep consistent handling. (We already enforced one column.)
    if (n == 0)
      % If we got here, treat as MATLAB does (above check should have errored),
      % but to be safe produce an empty 0xK double.
      D = zeros (0, K);
      return;
    end

    % Convert to numeric indices. Unknown/undefined map to NaN or 0 depending on API.
    % In MATLAB, <undefined> results in NaN rows in output.
    try
      idx = double (g);    % maps categories to 1..K, undefined -> NaN
    catch
      % fallback: construct mapping manually (less efficient)
      % Use grp2idx-like approach
      [~, ~, idx] = unique (g);
      idx = double (idx);
      % Unique will not produce NaN for undefined; handle undefined explicitly
      undef_mask = isundefined (g);
      idx (undef_mask) = NaN;
    end

    % Build matrix: rows with idx==NaN are rows of NaN(1,K)
    rows = (1:n)';
    D = zeros (n, K);

    % Build using sparse (skip NaNs)
    valid = ~isnan (idx);
    if any (valid)
      S = sparse (rows(valid), idx(valid), 1, n, K);
      D (:) = full (S);
    end

    % Replace rows where idx is NaN with NaN across all columns (MATLAB semantics)
    nan_rows = find (isnan (idx));
    if ~ isempty (nan_rows)
      D (nan_rows, :) = NaN;
    end

    D = double (D);
    return;
  end

  % --- NUMERIC / LEGACY branch ---
  % If the input is numeric (vector), create indicator columns for 1..max(g)
  if (isnumeric (g) && isvector (g))
    g = g(:);                 % ensure column
    if isempty (g)
      D = zeros (0, 0);
      return;
    end
    % If g has non-integer values, MATLAB implicitly coerces to integer indices
    % by using unique? Historically dummyvar expects group indices (positive ints).
    % We'll follow MATLAB's numeric behavior: use max(g) as number of columns.
    % Construct column count K = max(g)
    K = max (g);
    if ~ isreal (K) || K < 0
      error ("dummyvar:InvalidInput", "Numeric grouping must produce a positive integer number of groups.");
    end
    K = double (K);
    n = numel (g);
    % Build sparse/dense
    rows = (1:n)';
    idx = round (g);          % keep integer mapping
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

  % --- FALLBACK: unsupported input types ---
  error ("dummyvar:UnsupportedType", "dummyvar requires a numeric vector or a categorical array.");
end

## Test dummyvar behavior (MATLAB-compatible)

%!test
%! % numeric grouping vector
%! g = [1;2;1;3;2];
%! D = dummyvar(g);
%! assert(isequal(D, [1 0 0; 0 1 0; 1 0 0; 0 0 1; 0 1 0]));

%!test
%! % categorical with universe -> columns for each category in same order
%! g = categorical({'a';'b';'a'}, {'a','b','c'});
%! D = dummyvar(g);
%! assert(size(D,2) == numel(categories(g)));
%! assert(all(D(:,1) == [1;0;1]));
%! assert(all(D(:,2) == [0;1;0]));
%! assert(all(D(:,3) == [0;0;0]));

%!test
%! % categorical with <undefined> -> row of NaNs
%! g = categorical({'a'; ''; 'b'}, {'a','b','c'});
%! D = dummyvar(g);
%! assert(all(isnan(D(2,:))));
%! assert(all(D(1,:) == [1 0 0]));
%! assert(all(D(3,:) == [0 1 0]));

%!test
%! % empty categorical -> MATLAB-style error
%! g = categorical({}, {'a','b'});
%! assert (throws (@() dummyvar(g)));

%!test
%! % table column input
%! G = categorical({'a'; 'b'; 'a'}, {'a','b','c'});
%! T = table(G, [10;20;30], 'VariableNames', {'G','Val'});
%! D = dummyvar(T.G);
%! assert(size(D,2) == numel(categories(G)));
