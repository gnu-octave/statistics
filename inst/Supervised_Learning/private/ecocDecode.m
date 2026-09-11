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
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Private Function} {[@var{NegLoss}, @var{errmsg}] =} ecocDecode (@var{S}, @var{M}, @var{lossname}, @var{decoding}, @var{range})
##
## Decode the scores of the binary learners of an error correcting code.
##
## @var{S} is the @math{NxL} matrix of the scores each binary learner gave
## the class it was trained to call +1, @var{M} the @math{KxL} coding matrix,
## and @var{NegLoss} the @math{NxK} negated average loss, the largest entry of
## a row naming the class predicted.
##
## @var{range} is the interval the learners score on, @code{[-Inf, Inf]} for a
## margin and @code{[0, 1]} for a posterior.  A score is read as
## @math{u = M(k,j) s} on the first and @math{u = M(k,j) (2s-1)} on the
## second, which is the whole of why a loss suits one range and not the other.
##
## Every loss is scaled so that it is exactly @math{0.5} at @math{u = 0}, so a
## class a column leaves out costs the same wherever it appears.  Under
## @qcode{'lossbased'} the average runs over every column, those left out
## included, and under @qcode{'lossweighted'} it is weighted by
## @math{|M(k,j)|}, which drops them.
##
## @var{errmsg} is empty when the scores were decoded, and otherwise the body
## of a message the caller raises under its own name.
##
## @seealso{CompactClassificationECOC, fitcecoc}
## @end deftypefn

function [NegLoss, errmsg] = ecocDecode (S, M, lossname, decoding, range)

  NegLoss = [];
  errmsg = '';

  posterior = isequal (range, [0, 1]);
  if (posterior)
    allowed = {'hamming', 'quadratic'};
    spread = " [0,1]";
  else
    allowed = {'binodeviance', 'exponential', 'hamming', 'hinge', ...
               'linear', 'logit'};
    spread = " (-Inf,+Inf)";
  endif
  if (! any (strcmp (lossname, allowed)))
    ## strcat strips trailing whitespace from every argument, so the space
    ## before the range rides on SPREAD rather than standing on its own.
    errmsg = strcat ("you cannot use '", lossname, "' loss for binary", ...
                     " learners with response in the range", spread, ".");
    return;
  endif

  ## A posterior is read across the same interval a margin is, which is what
  ## the documented form of the quadratic loss already carries.
  if (posterior)
    F = 2 * S - 1;
  else
    F = S;
  endif

  K = rows (M);
  L = columns (M);
  N = rows (F);
  NegLoss = zeros (N, K);
  W = abs (M);

  for k = 1:K
    U = F .* M(k,:);
    G = binaryLoss (U, lossname);
    if (strcmp (decoding, 'lossbased'))
      NegLoss(:,k) = - sum (G, 2) / L;
    else
      NegLoss(:,k) = - (G * W(k,:)') / sum (W(k,:));
    endif
  endfor

endfunction

## The binary losses, each scaled so that it is 0.5 at u = 0.  That scaling is
## what the halves and the 2*log (2) are for, and it is the only thing that
## makes a column a class sits out cost the same as any other.  Note that the
## logit takes -u where the binomial deviance takes -2u; the two are otherwise
## the same expression.
function G = binaryLoss (U, lossname)

  switch (lossname)
    case 'binodeviance'
      G = log (1 + exp (-2 * U)) / (2 * log (2));
    case 'exponential'
      G = exp (-U) / 2;
    case 'hamming'
      G = (1 - sign (U)) / 2;
    case 'hinge'
      G = max (0, 1 - U) / 2;
    case 'linear'
      G = (1 - U) / 2;
    case 'logit'
      G = log (1 + exp (-U)) / (2 * log (2));
    case 'quadratic'
      G = (1 - U) .^ 2 / 2;
  endswitch

endfunction
