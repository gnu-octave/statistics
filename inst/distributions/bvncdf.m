## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} @var{p} = bvncdf (@var{x}, @var{mu}, @var{sigma})
## @deftypefnx {statistics} @var{p} = bvncdf (@var{x}, [], @var{sigma})
##
## Bivariate normal cumulative distribution function (CDF).
##
## @code{@var{p} = bvncdf (@var{x}, @var{mu}, @var{sigma})} will compute the
## bivariate normal cumulative distribution function of @var{x} given a mean
## @var{mu}, which must be a scalar, and a 2x2 @var{sigma} covariance matrix,
## which must be positive definite.
##
## @seealso{mvncdf}
## @end deftypefn

## Code adapted from Thomas H. Jørgensen's work in BVNcdf.m function retrieved
## from https://www.tjeconomics.com/code/

function p = bvncdf (x, mu, sigma)
  ## Check input arguments and add defaults
  if (size (x, 2) != 2)
    error (strcat (["bvncdf: X must be an Nx2 matrix with each variable"], ...
                   [" a column vector."]));
  endif
  if (isempty (mu))
    mu = [0,0];
  endif
  if (numel (sigma) == 1)
    sigma = sigma * ones (2,2);
    sigma(1,1) = 1;
    sigma(2,2) = 1;
  elseif (numel (sigma) != 4)
    error (strcat (["bvncdf: the covariance matrix must be either a"], ...
                   [" scalar or a 2x2 matrix."]));
  endif
  ## Test for symmetric positive definite covariance matrix
  [~, err] = chol (sigma);
  if (err != 0)
    error (strcat (["bvncdf: the covariance matrix is not positive"], ...
                   [" definite and/or symmetric."]));
  endif
  dh = (x(:,1) - mu(:,1)) / sqrt (sigma(1,1));
  dk = (x(:,2) - mu(:,2)) / sqrt (sigma(2,2));
  r = sigma(1,2) / sqrt (sigma(1,1) * sigma(2,2));
  p = NaN (size (dh));
  p(dh == Inf & dk == Inf)   = 1;
  p(dk == Inf) = 0.5 * erfc (- dh(dk == Inf) / sqrt (2));
  p(dh == Inf) = 0.5 * erfc (- dh(dk == Inf) / sqrt (2));
  p(dh == -Inf | dk == -Inf) = 0;
  ind = (dh > -Inf & dh < Inf & dk > -Inf & dk < Inf);
  ## For p(x1 < dh, x2 < dk, r)
  if (sum (ind) > 0)
    p(ind) = calculate_bvncdf (-dh(ind), -dk(ind), r);
  endif
endfunction

function p = calculate_bvncdf (dh,dk,r)
  if (abs (r) < 0.3)
    lg = 3;
    ## Gauss Legendre points and weights, n =  6
    w = [0.1713244923791705, 0.3607615730481384, 0.4679139345726904];
    x = [0.9324695142031522, 0.6612093864662647, 0.2386191860831970];
  elseif (abs (r) < 0.75)
    lg = 6;
    ## Gauss Legendre points and weights, n = 12
    w = [.04717533638651177, 0.1069393259953183, 0.1600783285433464, ...
         0.2031674267230659, 0.2334925365383547, 0.2491470458134029];
    x = [0.9815606342467191, 0.9041172563704750, 0.7699026741943050, ...
         0.5873179542866171, 0.3678314989981802, 0.1252334085114692];
  else
    lg = 10;
    ## Gauss Legendre points and weights, n = 20
    w = [.01761400713915212, .04060142980038694, .06267204833410906, ...
         .08327674157670475, 0.1019301198172404, 0.1181945319615184, ...
         0.1316886384491766, 0.1420961093183821, 0.1491729864726037, ...
         0.1527533871307259];
    x = [0.9931285991850949, 0.9639719272779138, 0.9122344282513259, ...
         0.8391169718222188, 0.7463319064601508, 0.6360536807265150, ...
         0.5108670019508271, 0.3737060887154196, 0.2277858511416451, ...
         0.07652652113349733];
  endif
  dim1 = ones (size (dh, 1), 1);
  dim2 = ones (1, lg);
  hk = dh .* dk;
  bvn = dim1 * 0;
  phi_dh =  0.5 * erfc (dh / sqrt (2));
  phi_dk =  0.5 * erfc (dk / sqrt (2));
  if (abs(r) < 0.925)
    hs = (dh .* dh + dk .* dk) / 2;
    asr = asin (r);
    sn1 = sin (asr * (1 - x) / 2);
    sn2 = sin (asr * (1 + x) / 2);
    bvn = sum ((dim1 * w) .* exp (((dim1 * sn1) .* (hk * dim2) - ...
              hs * dim2) ./ (1 - dim1 * (sn1 .^ 2))) + ...
              (dim1 * w) .* exp (((dim1 * sn2) .* (hk * dim2) - ...
              hs * dim2) ./ (1 - dim1 * (sn2 .^ 2))), 2) * ...
              asr / (4 * pi) + phi_dh .* phi_dk;
  else
    twopi = 2 * pi;
    if r < 0
      dk = -dk;
      hk = -hk;
    endif
    if abs(r) < 1
      as = (1 - r) * (1 + r);
      a = sqrt (as);
      bs = (dh - dk) .^ 2;
      c = (4 - hk) / 8;
      d = (12 - hk) / 16;
      asr = - (bs ./ as + hk) / 2;
      ind = asr > -100;
      bvn(ind) = a * exp (asr(ind)) .* (1 - (c(ind) .* (bs(ind) - as)) ...
                  .* (1 - d(ind) .* bs(ind) / 5) /3 ...
                  + (c(ind) .* d(ind)) .* as .^ 2 / 5 );
      ind = hk > -100;
      b = sqrt (bs);
      phi_ba = 0.5 * erfc ((b/a) / sqrt (2));
      sp = sqrt (twopi) * phi_ba;
      bvn(ind) = bvn(ind) - (exp (-hk(ind) / 2) .* sp(ind)) ...
                 .* b(ind) .* (1 - c(ind) .* bs(ind) ...
                 .* (1 - d(ind) .* bs(ind) / 5) /3);
      a = a/2;
      for is = -1:2:1
        xs = (a + a * is * x) .^ 2;
        rs = sqrt (1 - xs);
        asr1 = - ((bs * dim2) ./ (dim1 * xs) + hk * dim2) / 2;
        ind1 = (asr1 > -100);
        sp1 = (1 + (c * dim2) .* (dim1 * xs) .* ...
              (1 + (d * dim2) .* (dim1 * xs)));
        ep1 = exp (- (hk * dim2) .* (1 - dim1 * rs) ./ ...
                  (2 * (1 + dim1 * rs))) ./ (dim1 * rs);
        bvn = bvn + sum (a .* (dim1 * w) .* exp (asr1 .* ind1) ...
                           .* (ep1 .* ind1 - sp1 .* ind1), 2);
      endfor
      bvn = -bvn/twopi;
    endif
    if (r > 0)
      tmp = max (dh, dk);
      bvn =  bvn + 0.5 * erfc (tmp / sqrt(2));
    elseif (r < 0)
      phi_dh =  0.5 * erfc (dh / sqrt (2));
      phi_dk =  0.5 * erfc (dk / sqrt (2));
      bvn = - bvn + max (0, phi_dh - phi_dk);
    endif
  endif
  p = max (0, min (1, bvn));
endfunction

%!demo
%! mu = [1, -1];
%! sigma = [0.9, 0.4; 0.4, 0.3];
%! [X1, X2] = meshgrid (linspace (-1, 3, 25)', linspace (-3, 1, 25)');
%! x = [X1(:), X2(:)];
%! p = bvncdf (x, mu, sigma);
%! Z = reshape (p, 25, 25);
%! surf (X1, X2, Z);
%! title ("Bivariate Normal Distribution");
%! ylabel "X1"
%! xlabel "X2"

%!test
%! mu = [1, -1];
%! sigma = [0.9, 0.4; 0.4, 0.3];
%! [X1,X2] = meshgrid (linspace (-1, 3, 25)', linspace (-3, 1, 25)');
%! x = [X1(:), X2(:)];
%! p = bvncdf (x, mu, sigma);
%! p_out = [0.00011878988774500, 0.00034404112322371, ...
%!          0.00087682502191813, 0.00195221905058185, ...
%!          0.00378235566873474, 0.00638175749734415, ...
%!          0.00943764224329656, 0.01239164888125426, ...
%!          0.01472750274376648, 0.01623228313374828]';
%! assert (p([1:10]), p_out, 1e-16);
%!test
%! mu = [1, -1];
%! sigma = [0.9, 0.4; 0.4, 0.3];
%! [X1,X2] = meshgrid (linspace (-1, 3, 25)', linspace (-3, 1, 25)');
%! x = [X1(:), X2(:)];
%! p = bvncdf (x, mu, sigma);
%! p_out = [0.8180695783608276, 0.8854485749482751, ...
%!          0.9308108777385832, 0.9579855743025508, ...
%!          0.9722897881414742, 0.9788150170059926, ...
%!          0.9813597788804785, 0.9821977956568989, ...
%!          0.9824283794464095, 0.9824809345614861]';
%! assert (p([616:625]), p_out, 2e-16);
%!error bvncdf (randn (25,3), [], [1, 1; 1, 1]);
%!error bvncdf (randn (25,2), [], [1, 1; 1, 1]);
%!error bvncdf (randn (25,2), [], ones (3, 2));
