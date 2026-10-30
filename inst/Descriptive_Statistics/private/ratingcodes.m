## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Private Function} {[@var{codes}, @var{cats}, @var{errmsg}] =} ratingcodes (@var{R}, @var{cats})
##
## Code the ratings of several raters against one set of categories.
##
## @var{R} is a cell array holding the ratings of each rater, one vector
## each, all of the same kind and length: numeric or logical, categorical, or
## text (a character matrix with one rating per row, a cell array of
## character vectors, or a string array).  @var{codes} has one column per
## rater, holding the index of each rating in @var{cats}, and @qcode{NaN}
## where the rating is missing: @qcode{NaN}, an undefined category, an empty
## text or a missing string.
##
## Given @var{cats} empty, the categories are the sorted values of numeric
## ratings, the categories of categorical ratings in their own order, unused
## ones included, and the sorted texts of text ratings.  Otherwise @var{cats}
## lists them, in the order to be used; numeric for numeric ratings, text
## or categorical for the others.  @var{cats} is returned as a numeric
## column or a cell array of character vectors.  @var{errmsg} is empty, or
## says what is wrong with the input.
## @end deftypefn

function [codes, cats, errmsg] = ratingcodes (R, cats)

  codes = [];
  errmsg = "";

  ## The kind of each rater's ratings, and their number
  k = numel (R);
  kind = cell (1, k);
  n = zeros (1, k);
  for j = 1:k
    r = R{j};
    if ((isnumeric (r) || islogical (r)) && isreal (r) && isvector (r))
      kind{j} = 'numeric';
      n(j) = numel (r);
    elseif (iscategorical (r) && isvector (r))
      kind{j} = 'categorical';
      n(j) = numel (r);
    elseif (ischar (r) && ndims (r) == 2)
      kind{j} = 'text';
      n(j) = rows (r);
    elseif ((iscellstr (r) || isstring (r)) && isvector (r))
      kind{j} = 'text';
      n(j) = numel (r);
    else
      errmsg = "ratings must be numeric, logical, categorical or text vectors.";
      return;
    endif
  endfor
  if (! all (strcmp (kind, kind{1})))
    errmsg = "ratings of every rater must be of the same kind.";
    return;
  endif
  if (any (n != n(1)))
    errmsg = "every rater must rate the same number of subjects.";
    return;
  endif

  ## All ratings in one column, rater after rater, with what is missing
  switch (kind{1})
    case 'numeric'
      all_r = zeros (0, 1);
      for j = 1:k
        all_r = [all_r; double(R{j}(:))];
      endfor
      miss = isnan (all_r);
    case 'categorical'
      all_r = R{1}(:);
      for j = 2:k
        all_r = [all_r; R{j}(:)];
      endfor
      if (isempty (cats))
        cats = categories (all_r);
      endif
      miss = isundefined (all_r);
      all_r = cellstr (all_r);
    case 'text'
      all_r = cell (0, 1);
      for j = 1:k
        r = R{j};
        if (isstring (r))
          m = ismissing (r(:));
          r = cellstr (r(:));
          r(m) = {''};
        else
          r = cellstr (r);
        endif
        all_r = [all_r; r(:)];
      endfor
      miss = cellfun (@isempty, all_r);
  endswitch

  ## The categories, given or found
  if (strcmp (kind{1}, 'numeric'))
    if (isempty (cats))
      cats = unique (all_r(! miss));
    elseif (! ((isnumeric (cats) || islogical (cats)) && isreal (cats)
               && isvector (cats) && ! any (isnan (cats))))
      errmsg = "'Categories' must be a numeric vector for numeric ratings.";
      return;
    endif
    cats = double (cats(:));
  else
    if (isempty (cats))
      cats = unique (all_r(! miss));
    elseif (iscategorical (cats) || isstring (cats) || ischar (cats))
      cats = cellstr (cats);
    elseif (! iscellstr (cats))
      errmsg = "'Categories' must be text for text or categorical ratings.";
      return;
    endif
    cats = cats(:);
  endif
  if (numel (unique (cats)) != numel (cats))
    errmsg = "'Categories' must not repeat a category.";
    return;
  endif

  [tf, idx] = ismember (all_r, cats);
  if (any (! tf & ! miss))
    errmsg = "a rating is not one of 'Categories'.";
    return;
  endif
  codes = idx;
  codes(miss) = NaN;
  codes = reshape (codes, n(1), k);

endfunction
