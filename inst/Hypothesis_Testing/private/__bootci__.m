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
## @deftypefn  {Private Function} {@var{ci} =} __bootci__ (@var{stat}, @var{S}, @var{alpha}, @var{B}, @var{scheme})
##
## Bias-corrected and accelerated (BCa) bootstrap confidence interval of a
## statistic of one or more samples.
##
## @var{S} is a cell array of column vectors, the samples, and @var{stat} a
## function handle taking them as separate arguments, @code{@var{stat}
## (@var{S}@{:@})}.  Each argument it is given may be a matrix whose columns are
## resampled copies of that sample, and it returns a row vector holding the
## statistic of each column.  @var{alpha} sets the level of the two-sided
## interval, @math{100 (1 - @var{alpha})} percent, and @var{B} is the number of
## replicates.  @var{scheme} says how the samples are resampled:
##
## @table @asis
## @item @qcode{'rows'}
## The samples have the same length and are resampled by the same indices, as
## the rows of a matrix: one sample, paired samples, or the blocks of a
## design.
## @item @qcode{'strata'}
## Each sample is resampled on its own, keeping its size.
## @item @qcode{'pooled'}
## The observations of all samples are resampled together, each keeping its
## sample, so the sizes vary between replicates; @var{stat} is then given
## column vectors, one replicate at a time.  A replicate in which a sample
## comes out empty is left out.
## @end table
##
## The acceleration is estimated by the jackknife: over the rows for
## @qcode{'rows'}, over the pooled observations for @qcode{'pooled'}, and over
## each sample for @qcode{'strata'}, the influence values of each weighted by
## the size of its sample.  @var{ci} is @code{[NaN, NaN]} where the statistic
## of the samples or every replicate is @qcode{NaN}.
##
## @end deftypefn

function ci = __bootci__ (stat, S, alpha, B, scheme)

  ns = cellfun (@numel, S(:))';
  t0 = stat (S{:});
  bs = NaN (B, 1);
  ## Replicates in blocks of about a million values
  step = max (1, floor (1e6 / sum (ns)));
  switch (scheme)
    case 'rows'
      n = ns(1);
      for b1 = 1:step:B
        b = b1:min (b1 + step - 1, B);
        I = randi (n, n, numel (b));
        R = cellfun (@(v) pick (v, I), S, 'UniformOutput', false);
        bs(b) = stat (R{:});
      endfor
      u = {jackInfluence(@(i) jackRows (stat, S, i), n)};
    case 'strata'
      for b1 = 1:step:B
        b = b1:min (b1 + step - 1, B);
        R = cell (size (S));
        for g = 1:numel (S)
          R{g} = pick (S{g}, randi (ns(g), ns(g), numel (b)));
        endfor
        bs(b) = stat (R{:});
      endfor
      ## Jackknife influence values of each sample, weighted by its size
      u = cell (1, numel (S));
      for g = 1:numel (S)
        u{g} = (ns(g) - 1) * jackInfluence (@(i) jackStratum (stat, S, g, i), ...
                                            ns(g)) / ns(g);
      endfor
    case 'pooled'
      v = vertcat (S{:});
      g = repelem ((1:numel (S))', ns(:));
      n = numel (v);
      for b1 = 1:step:B
        b = b1:min (b1 + step - 1, B);
        I = randi (n, n, numel (b));
        for k = 1:numel (b)
          gb = g(I(:,k));
          if (all (accumarray (gb, 1, [numel(S), 1]) > 0))
            bs(b(k)) = stat (splitPooled (v(I(:,k)), gb, numel (S)){:});
          endif
        endfor
      endfor
      jk = NaN (n, 1);
      for i = 1:n
        k = [1:i-1, i+1:n]';
        if (all (accumarray (g(k), 1, [numel(S), 1]) > 0))
          jk(i) = stat (splitPooled (v(k), g(k), numel (S)){:});
        endif
      endfor
      jk = jk(! isnan (jk));
      u = {mean(jk) - jk};
  endswitch
  bs = bs(! isnan (bs));
  if (isempty (bs) || isnan (t0))
    ci = [NaN, NaN];
    return;
  endif

  ## Bias correction and acceleration
  z0 = norminv (mean (bs < t0) + mean (bs == t0) / 2);
  u = vertcat (u{:});
  u = u(isfinite (u));
  den = sum (u .^ 2);
  if (den > 0)
    a = sum (u .^ 3) / (6 * den ^ 1.5);
  else
    a = 0;
  endif
  zq = norminv ([alpha / 2, 1 - alpha / 2]);
  p = normcdf (z0 + (z0 + zq) ./ (1 - a * (z0 + zq)));
  ci = quantile (bs, p(:), 1, 5)';

endfunction

## The elements of the vector V at the indices I, in the shape of I
function out = pick (v, I)
  out = reshape (v(I), size (I));
endfunction

## The pooled values V split back into the NS samples their labels G name
function out = splitPooled (v, g, ns)
  out = cell (1, ns);
  for k = 1:ns
    out{k} = v(g == k);
  endfor
endfunction

## The statistic with every sample's rows left out as the columns of I say
function out = jackRows (stat, S, I)
  R = cellfun (@(v) pick (v, I), S, 'UniformOutput', false);
  out = stat (R{:});
endfunction

## The statistic with sample G left out as the columns of I say, and every
## other sample whole
function out = jackStratum (stat, S, g, I)
  R = cellfun (@(v) repmat (v, 1, columns (I)), S, 'UniformOutput', false);
  R{g} = pick (S{g}, I);
  out = stat (R{:});
endfunction

## The jackknife deviations mean (jk) - jk of a statistic over N
## observations, FCN taking a matrix of indices whose columns each leave one
## observation out, in blocks
function u = jackInfluence (fcn, n)
  jk = zeros (n, 1);
  step = max (1, floor (1e6 / n));
  for i1 = 1:step:n
    i = i1:min (i1 + step - 1, n);
    ## Column k holds 1:n without i(k)
    K = repmat ((1:n)', 1, numel (i));
    K(sub2ind (size (K), i, 1:numel (i))) = 0;
    K = reshape (K(K > 0), n - 1, numel (i));
    jk(i) = fcn (K);
  endfor
  u = mean (jk) - jk;
endfunction
