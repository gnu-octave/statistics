## Copyright (C) 2013 Nir Krakauer <nkrakauer@ccny.cuny.edu>
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
## @deftypefn {Function File} {@var{h}, @var{p}, @var{stats} =} runstest (@var{x}, @var{v})
## Runs test for detecting serial correlation in the vector @var{x}.
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{x} is the vector of given values.
## @item
## @var{v} is the value to subtract from @var{x} to get runs (defaults to @code{median(x)})
## @end itemize
##
## @subheading Return values
##
## @itemize @bullet
## @item
## @var{h} is true if serial correlation is detected at the 95% confidence level (two-tailed), false otherwise.
## @item
## @var{p} is the probablity of obtaining a test statistic of the magnitude found under the null hypothesis of no serial correlation.
## @item
## @var{stats} is the structure containing as fields the number of runs @var{nruns}; the numbers of positive and negative values of @code{x - v}, @var{n1} and @var{n0}; and the test statistic @var{z}.
## 
## @end itemize
##
## Note: the large-sample normal approximation is used to find @var{h} and @var{p}. This is accurate if @var{n1}, @var{n0} are both greater than 10.
##
## Reference: 
## NIST Engineering Statistics Handbook, 1.3.5.13. Runs Test for Detecting Non-randomness, http://www.itl.nist.gov/div898/handbook/eda/section3/eda35d.htm
##
## @seealso{}
## @end deftypefn

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>
## Description: Runs test for detecting serial correlation

function [h, p, stats] = runstest (x, x2)

  # Check arguments
  if (nargin < 1)
    print_usage;
  endif
  
  if nargin > 1 && isnumeric(x2)
    v = x2;
  else
    v = median(x);
  endif
  
  x = x(~isnan(x)); #delete missing values
  x = sign(x - v);
  x = x(x ~= 0); #delete any zeros
  
  R = sum((x(1:(end-1)) .* x(2:end)) < 0) + 1; #number of runs
  
  #expected number of runs for an iid sequence
  n1 = sum(x > 0);
  n2 = sum(x < 0);
  R_bar = 1 + 2*n1*n2/(n1 + n2);
  
  #standard deviation of number of runs for an iid sequence
  s_R = sqrt(2*n1*n2*(2*n1*n2 - n1 - n2)/((n1 + n2)^2 * (n1 + n2 - 1)));
 
  #desired significance level
  alpha = 0.05;
   
  Z = (R - R_bar) / s_R; #test statistic
 
  p = 2 * normcdf(-abs(Z));

  h = p < alpha;

  if nargout > 2
    stats.nruns = R;
    stats.n1 = n1;
    stats.n0 = n2;
    stats.z = Z;
  endif
  
endfunction



%!test
%! data = [-213       -564       -35       -15       141       115       -420       -360       203       -338       -431       194       -220       -513       154       -125       -559       92       -21       -579       -52       99       -543       -175       162       -457       -346       204       -300       -474       164       -107       -572       -8       83       -541       -224       180       -420       -374       201       -236       -531       83       27       -564       -112       131       -507       -254       199       -311       -495       143       -46       -579       -90       136       -472       -338       202       -287       -477       169       -124       -568       17       48       -568       -135       162       -430       -422       172       -74       -577       -13       92       -534       -243       194       -355       -465       156       -81       -578       -64       139       -449       -384       193       -198       -538       110       -44       -577       -6       66       -552       -164       161       -460       -344       205       -281       -504       134       -28       -576       -118       156       -437       -381       200       -220       -540       83       11       -568       -160       172       -414       -408       188       -125       -572       -32       139       -492       -321       205       -262       -504       142       -83       -574       0       48       -571       -106       137       -501       -266       190       -391       -406       194       -186       -553       83       -13       -577       -49       103       -515       -280       201       300       -506       131       -45       -578       -80       138       -462       -361       201       -211       -554       32       74       -533       -235       187       -372       -442       182       -147       -566       25       68       -535       -244       194       -351       -463       174       -125       -570       15       72       -550       -190       172       -424       -385       198       -218       -536       96]; #NIST beam deflection data, http://www.itl.nist.gov/div898/handbook/eda/section4/eda425.htm
%! [h, p, stats] = runstest (data);
%! expected_h = true;
%! expected_p = 0.0070646;
%! expected_z = 2.6938;
%! assert (h, expected_h);
%! assert (p, expected_p, 1E-6);
%! assert (stats.z, expected_z, 1E-4);
