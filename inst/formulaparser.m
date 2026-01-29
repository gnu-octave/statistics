function varargout = formulaparser (varargin)

  if (nargin < 1)
    error ("formulaparser: Input formula string is required.");
  elseif (nargin > 3)
    error ("formulaparser: Too many input arguments.");
  endif

  formula_str = varargin{1};
  if (nargin > 1)
    mode = varargin{2};
  else
    mode = "tokenize";
  endif

  switch (lower (mode))
    case "tokenize"
      varargout{1} = run_lexer (formula_str);

    case "parse"
      tokens = run_lexer (formula_str);
      varargout{1} = run_parser (tokens);

    case "expand"
      tokens = run_lexer (formula_str);
      tree = run_parser (tokens);
      varargout{1} = run_expander (tree);

    case "matrix"
      tokens = run_lexer (formula_str);
      tree = run_parser (tokens);
      expanded = run_expander (tree);
      varargout{1} = run_schema_builder (expanded);

    case "model_matrix"
      ## Process Formula
      tokens = run_lexer (formula_str);
      tree = run_parser (tokens);

      ## Auto-wrap in Formula (~) if strictly a term list
      if (! strcmp (tree.type, "OPERATOR") || ! strcmp (tree.value, "~"))
        wrapper.type = "OPERATOR";
        wrapper.value = "~";
        wrapper.left = [];
        wrapper.right = tree;
        tree = wrapper;
      endif

      expanded = run_expander (tree);
      schema = run_schema_builder (expanded);

      ## Compile with Data
      if (nargin < 3)
        error ("formulaparser: 'model_matrix' mode requires a Data Table.");
      endif

      data_table = varargin{3};
      [X, y, names] = run_model_matrix_builder (schema, data_table);

      varargout{1} = X;
      if (nargout > 1), varargout{2} = y; endif
      if (nargout > 2), varargout{3} = names; endif

    otherwise
      error ("formulaparser: Unknown mode: %s", mode);
  endswitch

endfunction

## lexer
function tokens = run_lexer (formula_str)

  if (isempty (formula_str))
    tokens = struct ("type", {}, "value", {}, "pos", {});
    return;
  endif

  str = char (formula_str);
  n = length (str);
  tokens(n) = struct ("type", "", "value", "", "pos", 0);
  tok_idx = 0;
  i = 1;

  while (i <= n)
    c = str(i);
    start_pos = i;

    if (isspace (c))
      i = i + 1;
      continue;
    endif

    ## Identifiers (Factors)
    if (isletter (c))
      val = c;
      i = i + 1;
      while (i <= n)
        next_c = str(i);
        if (isletter (next_c) || (next_c >= '0' && next_c <= '9') ...
            || next_c == '_')
          val = [val, next_c];
          i = i + 1;
        else
          break;
        endif
      endwhile
      tok_idx = tok_idx + 1;
      tokens(tok_idx) = create_token ("IDENTIFIER", val, start_pos);
      continue;
    endif

    ## Numbers (Implicit Intercept '1' or others)
    if (c >= '0' && c <= '9')
      val = c;
      i = i + 1;
      while (i <= n)
        next_c = str(i);
        if (next_c >= '0' && next_c <= '9')
          val = [val, next_c];
          i = i + 1;
        elseif (next_c == '.')
          ## Simple float support
          if (i < n && str(i+1) >= '0' && str(i+1) <= '9')
            val = [val, next_c];
            i = i + 1;
          else
            break;
          endif
        else
          break;
        endif
      endwhile
      tok_idx = tok_idx + 1;
      tokens(tok_idx) = create_token ("NUMBER", val, start_pos);
      continue;
    endif

    ## Operators
    type = "";
    val = c;
    skip = 0;

    switch (c)
      case '-'
        if (i < n)
          if (i < n && str(i+1) == '/')
            type = "OP_MINUS_MARGIN"; val = "-/"; skip = 1;
          elseif (i < n && str(i+1) == '*')
            type = "OP_MINUS_CLEAN"; val = "-*"; skip = 1;
          else
            type = "OP_MINUS";
          endif
        else
          type = "OP_MINUS";
        endif
      case '*'
        type = "OP_CROSS";
      case '/'
        type = "OP_NEST";
      case '+'
        type = "OP_PLUS";
      case {'.', ':'}
        type = "OP_DOT";
      case '~'
        type = "SEPARATOR";
      case '('
        type = "LPAREN";
      case ')'
        type = "RPAREN";
      otherwise
        error ("formulaparser: Unexpected character '%s' at position %d", c, i);
    endswitch

    tok_idx = tok_idx + 1;
    tokens(tok_idx) = create_token (type, val, start_pos);
    i = i + 1 + skip;
  endwhile

  tokens = tokens(1:tok_idx);
  tokens(end+1) = create_token ("EOF", "EOF", i);

endfunction

function t = create_token (type, val, pos)
  t.type = type;
  t.value = val;
  t.pos = pos;
endfunction

## parser
function [tree, curr] = run_parser (tokens, curr, prec_limit)

  if (nargin < 2), curr = 1; endif
  if (nargin < 3), prec_limit = 0; endif

  n = length (tokens);
  if (curr > n)
    error ("formulaparser: Unexpected End Of Formula.");
  endif

  t = tokens(curr);
  curr = curr + 1;

  ## Basic Units
  if (strcmp (t.type, "IDENTIFIER") || strcmp (t.type, "NUMBER"))
    tree.type = t.type;
    tree.value = t.value;
    tree.left = [];
    tree.right = [];

  elseif (strcmp (t.type, "LPAREN"))
    [tree, curr] = run_parser (tokens, curr, 0);
    if (curr <= n && strcmp (tokens(curr).type, "RPAREN"))
      curr = curr + 1;
    else
      error ("formulaparser: Mismatched Parentheses. Missing ')'.");
    endif

  elseif (strcmp (t.type, "EOF"))
    error ("formulaparser: Unexpected End Of Formula.");
  else
    error ("formulaparser: Syntax Error. Unexpected token: '%s'", t.value);
  endif

  ## Operator Handling with Correct Wilkinson Precedence
  while (curr <= n)
    op_type = tokens(curr).type;
    op_prec = 0;

    ## Precedence Levels (Higher # = Tighter Binding)
    if (strcmp (op_type, "OP_DOT"))
      op_prec = 50;
    elseif (strcmp (op_type, "OP_NEST"))
      op_prec = 40;
    elseif (strcmp (op_type, "OP_CROSS"))
      op_prec = 30;
    elseif (strcmp (op_type, "OP_PLUS"))
      op_prec = 20;
    elseif (strncmp (op_type, "OP_MINUS", 8))
      op_prec = 10;
    elseif (strcmp (op_type, "SEPARATOR"))
      op_prec = 5;
    else
      break;
    endif

    if (op_prec <= prec_limit)
      break;
    endif

    op_val = tokens(curr).value;
    curr = curr + 1;

    ## Parse RHS with current precedence
    [right, curr] = run_parser (tokens, curr, op_prec);

    new_node.type = "OPERATOR";
    new_node.value = op_val;
    new_node.left = tree;
    new_node.right = right;
    tree = new_node;
  endwhile

endfunction

## expander
function result = run_expander (node)

  if (isempty (node))
    result = {};
    return;
  endif

  ## Terminals
  if (strcmp (node.type, "IDENTIFIER"))
    result = {{node.value}};
    return;
  elseif (strcmp (node.type, "NUMBER"))
    ## "1" implies intercept (empty set of factors)
    if (strcmp (node.value, "1"))
      result = {{}};
    else
      result = {{node.value}};
    endif
    return;
  endif

  if (strcmp (node.type, "OPERATOR"))
    ## Formula Separator
    if (strcmp (node.value, "~"))
      result.response = run_expander (node.left);

      ## Implicit Intercept Logic
      add_intercept = true;
      if (! isempty (node.right) && strcmp (node.right.type, "OPERATOR") ...
          && (strcmp (node.right.value, "-") ...
          || strcmp (node.right.value, "-/") ...
          || strcmp (node.right.value, "-*")))

        if (! isempty (node.right.right) ...
            && strcmp (node.right.right.type, "NUMBER") ...
            && strcmp (node.right.right.value, "1"))
          add_intercept = false;
        endif
      endif

      model_raw = run_expander (node.right);

      if (add_intercept)
        result.model = list_union ({{}}, model_raw);
      else
        result.model = model_raw;
      endif
      return;
    endif

    lhs = run_expander (node.left);
    rhs = run_expander (node.right);

    switch (node.value)
      case "+"
        result = list_union (lhs, rhs);

      case {".", ":"}
        result = list_product (lhs, rhs);

      case "*"  ## Crossing: A + B + A.B
        interaction = list_product (lhs, rhs);
        step1 = list_union (lhs, rhs);
        result = list_union (step1, interaction);

      case "/"  ## Nesting: L + FAC(L).M
        fac_L = get_fac (lhs);
        interaction = list_product (fac_L, rhs);
        result = list_union (lhs, interaction);

      case "-"  ## Simple Deletion
        result = list_difference (lhs, rhs, "exact");

      case "-*" ## Delete term + Higher Order Interactions
        result = list_difference (lhs, rhs, "clean");

      case "-/" ## Delete terms where T is marginal
        result = list_difference (lhs, rhs, "margin");

      otherwise
        error ("formulaparser: Unknown operator '%s'", node.value);
    endswitch
    return;
  endif

  error ("formulaparser: Corrupt Tree.");

endfunction

## set operations.
function C = list_union (A, B)
  raw_list = [A, B];
  C = simplify_term_list (raw_list);
endfunction

function C = list_product (A, B)
  C = {};
  idx = 1;
  for i = 1:length (A)
    for j = 1:length (B)
      ## Dot product merges factor sets: {A}.{B} -> {A,B}
      new_term = union (A{i}, B{j});
      C{idx} = new_term;
      idx = idx + 1;
    endfor
  endfor
  C = simplify_term_list (C);
endfunction

function C = list_difference (S, T, mode)
  if (isempty (S)), C = {}; return; endif
  if (isempty (T)), C = S; return; endif

  C = {};
  strS = terms_to_strings (S);
  strT = terms_to_strings (T);

  keep_mask = true (size (S));

  for i = 1:length (S)
    term_s = S{i};
    s_str = strS{i};

    for j = 1:length (T)
      term_t = T{j};
      t_str = strT{j};

      match = false;

      switch (mode)
        case "exact"
          if (strcmp (s_str, t_str))
            match = true;
          endif

        case "clean" ## Delete T and any S where T is a subset of S
          if (strcmp (s_str, t_str) || is_subset (term_t, term_s))
            match = true;
          endif

        case "margin" ## Delete S where T is subset of S (but not T itself)
          if (! strcmp (s_str, t_str) && is_subset (term_t, term_s))
            match = true;
          endif
      endswitch

      if (match)
        keep_mask(i) = false;
        break;
      endif
    endfor
  endfor

  C = S(keep_mask);
endfunction

function fac = get_fac (term_list)
  all_factors = {};
  for i = 1:length (term_list)
    all_factors = union (all_factors, term_list{i});
  endfor
  fac = {all_factors};
endfunction

function is_sub = is_subset (small_set, large_set)
  is_sub = all (ismember (small_set, large_set));
endfunction

function clean_list = simplify_term_list (raw_list)
  if (isempty (raw_list))
    clean_list = {}; return;
  endif
  str_sigs = terms_to_strings (raw_list);
  [~, unique_idx] = unique (str_sigs);
  clean_list = raw_list(sort (unique_idx));
endfunction

function strs = terms_to_strings (term_list)
  strs = cell (size (term_list));
  for i = 1:length (term_list)
    if (isempty (term_list{i}))
      strs{i} = "1";
    else
      sorted_factors = sort (term_list{i});
      strs{i} = strjoin (sorted_factors, ":");
    endif
  endfor
endfunction

## schema builder
function schema = run_schema_builder (expanded)

  ## Handle Struct vs Cell
  if (isstruct (expanded))
    if (isfield (expanded, "model"))
      rhs_terms = expanded.model;
    else
      rhs_terms = expanded.rhs;
    endif
    if (isfield (expanded, "response"))
      lhs_term = expanded.response;
    else
      lhs_term = expanded.lhs;
    endif
  else
    rhs_terms = expanded;
    lhs_term = {};
  endif

  function out = flatten_recursive (in_val)
    out = {};
    if (ischar (in_val) || isstring (in_val))
      out = {char(in_val)};
    elseif (iscell (in_val))
      for k = 1:numel (in_val)
        out = [out, flatten_recursive(in_val{k})];
      endfor
    endif
  endfunction

  ## Extract Variables
  all_vars = {};
  if (! isempty (lhs_term))
    all_vars = [all_vars, flatten_recursive(lhs_term)];
  endif

  cleaned_rhs = cell (length (rhs_terms), 1);
  for i = 1:length (rhs_terms)
    term_vars = flatten_recursive (rhs_terms{i});
    final_term_vars = {};
    for j = 1:length (term_vars)
      parts = strsplit (term_vars{j}, ":");
      final_term_vars = [final_term_vars, parts];
    endfor
    cleaned_rhs{i} = final_term_vars;
    all_vars = [all_vars, final_term_vars];
  endfor

  all_vars = unique (all_vars);
  ## Remove intercept marker from var list
  all_vars(strcmp (all_vars, "1")) = [];

  schema.VariableNames = all_vars;

  ## Identify Response
  schema.ResponseIdx = [];
  if (! isempty (lhs_term))
    flat_lhs = flatten_recursive (lhs_term);
    if (! isempty (flat_lhs))
      [found, idx] = ismember (flat_lhs{1}, all_vars);
      if (found), schema.ResponseIdx = idx; endif
    endif
  endif

  ## Build Terms Matrix (Vars x Terms)
  n_vars = length (all_vars);
  n_terms = length (cleaned_rhs);
  terms_mat = zeros (n_terms, n_vars);

  for i = 1:n_terms
    vars_in_this_term = cleaned_rhs{i};

    ## Check for Intercept (empty or "1")
    if (isempty (vars_in_this_term) ...
        || (length (vars_in_this_term) == 1 ...
        && strcmp (vars_in_this_term{1}, "1")))
      continue;
    endif

    [found, idx] = ismember (vars_in_this_term, all_vars);
    if (any (! found))
      error ("formulaparser: Unknown variable in term definition.");
    endif
    terms_mat(i, idx) = 1;
  endfor

  ## Wilkinson Sorting: Order by Order (Main effects, then 2-way, etc.)
  term_orders = sum (terms_mat, 2);
  M = [term_orders, terms_mat];

  [~, unique_idx] = unique (M, "rows");
  terms_mat = terms_mat(unique_idx, :);

  [~, sort_idx] = sortrows ([sum(terms_mat, 2), terms_mat]);
  schema.Terms = terms_mat(sort_idx, :);

endfunction

## model matrix builder.
function [X, y, col_names] = run_model_matrix_builder (schema, data)

  req_vars = schema.VariableNames;

  ## Data Validation & Masking
  if (isempty (req_vars))
    fnames = fieldnames (data);
    n_total = length (data.(fnames{1}));
    valid_mask = true (n_total, 1);
  else
    n_total = length (data.(req_vars{1}));
    valid_mask = true (n_total, 1);
    for i = 1:length (req_vars)
      col = data.(req_vars{i});
      if (isnumeric (col))
        valid_mask = valid_mask & ! isnan (col);
      endif
    endfor
  endif

  if (! isempty (schema.ResponseIdx))
    y_name = req_vars{schema.ResponseIdx};
    y_col = data.(y_name);
    if (isnumeric (y_col))
      valid_mask = valid_mask & ! isnan (y_col);
    endif
  endif

  n_rows = sum (valid_mask);

  ## Process Predictors (Numeric vs Categorical)
  var_info = struct ();
  for i = 1:length (req_vars)
    vname = req_vars{i};
    raw = data.(vname);

    if (iscell (raw)), raw = raw(valid_mask);
    else, raw = raw(valid_mask, :); endif

    if (isnumeric (raw))
      var_info.(vname).type = "numeric";
      var_info.(vname).data = raw;
    else
      if (! iscellstr (raw) && ! isstring (raw)), raw = cellstr (raw); endif
      [u, ~, idx] = unique (raw);
      var_info.(vname).type = "categorical";
      var_info.(vname).levels = u;
      var_info.(vname).indices = idx;
    endif
  endfor

  ## Build Design Matrix X
  X = [];
  col_names = {};

  ## Check for Intercept (row of zeros in Terms matrix)
  intercept_row_idx = find (sum (schema.Terms, 2) == 0);
  has_intercept = ! isempty (intercept_row_idx);

  n_terms = size (schema.Terms, 1);

  for i = 1:n_terms
    term_row = schema.Terms(i, :);
    vars_idx = find (term_row);

    ## Intercept Term
    if (isempty (vars_idx))
      X = [X, ones(n_rows, 1)];
      col_names = [col_names; "(Intercept)"];
      continue;
    endif

    current_block = ones (n_rows, 1);
    current_names = {""};

    for v = vars_idx
      vname = req_vars{v};
      info = var_info.(vname);

      if (strcmp (info.type, "numeric"))
        current_block = current_block .* info.data;
        for k = 1:length (current_names)
          if (isempty (current_names{k}))
            current_names{k} = vname;
          else
            current_names{k} = [current_names{k}, ":", vname];
          endif
        endfor
      else
        ## Categorical: Reference Coding (Corner Point)
        n_lev = length (info.levels);

        ## Drop first level if Intercept exists to avoid rank deficiency
        if (has_intercept)
          start_lev = 2;
          n_cols = n_lev - 1;
        else
          start_lev = 1;
          n_cols = n_lev;
        endif

        dummies = zeros (n_rows, n_cols);
        dum_names = {};

        for L = start_lev:n_lev
          col_idx = L - start_lev + 1;
          dummies(:, col_idx) = (info.indices == L);
          dum_names = [dum_names; sprintf("%s_%s", vname, ...
                       char (info.levels{L}))];
        endfor

        ## Cartesian Product of current block and new dummies
        next_block = [];
        next_names = {};

        for c1 = 1:size (current_block, 2)
          for c2 = 1:size (dummies, 2)
            next_block = [next_block, current_block(:, c1) .* dummies(:, c2)];

            n1 = current_names{c1};
            n2 = dum_names{c2};

            if (isempty (n1))
              next_names = [next_names; n2];
            else
              next_names = [next_names; [n1, ":", n2]];
            endif
          endfor
        endfor
        current_block = next_block;
        current_names = next_names;
      endif
    endfor

    X = [X, current_block];
    col_names = [col_names; current_names];
  endfor

  ## Extract Response
  y = [];
  if (! isempty (schema.ResponseIdx))
    y_name = req_vars{schema.ResponseIdx};
    raw_y = data.(y_name);
    if (iscell (raw_y)), y = raw_y(valid_mask);
    else, y = raw_y(valid_mask, :); endif
  endif

endfunction