## Copyright (C) 2003 Andy Adler
##
## This file is part of Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, write to the Free
## Software Foundation, 59 Temple Place - Suite 330, Boston, MA
## 02111-1307, USA.

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{pval}, @var{f}, @var{df_b}, @var{df_w}] =} anovan (@var{data}, @var{grps})
## Perform a multi-way analysis of variance (ANOVA).  The goal is to test
## whether the population means of data taken from @var{k} different
## groups are all equal.
##
## Data is a single vector @var{data} with groups specified by
## a corresponding matrix of group labels @var{grps}, where @var{grps}
## has the same number of rows as @var{data}. For example, if
## @var{data} = [1.1;1.2]; @var{grps}= [1,2,1; 1,5,2];
## then data point 1.1 was measured under conditions 1,2,1 and
## data point 1.2 was measured under conditions 1,5,2.
## Note that groups do not need to be sequentially numbered.
##
## Under the null of constant means, the statistic @var{f} follows an F
## distribution with @var{df_b} and @var{df_w} degrees of freedom.
##
## The p-value (1 minus the CDF of this distribution at @var{f}) is
## returned in @var{pval}.
##
## If no output argument is given, the standard one-way ANOVA table is
## printed.
##
## NOTE: this function has not yet been tested with > 2 ways
## @end deftypefn

## Author: Andy Adler <adler@site.uottawa.ca>
## Based on code by: KH <Kurt.Hornik@ci.tuwien.ac.at>

function [pval, f, df_b, df_w] = anovan (data, grps)

    if nargin <= 1
        usage ("anovan (data, grps)");
    end

    if ~isvector (data)
          error ("anova: for `anova (data, grps)', data must be a vector");
    endif

    nd = size (grps,1); # number of data points
    nw = size (grps,2); # number of anova "ways"
    if (~ isvector (data) || (length(data) ~= nd))
      error ("anova: grps must be a matrix of the same number of rows as data");
    endif

    [g,grp_map]   = relabel_groups (grps);
    max_interact  = length(grp_map);
    ng = length(grp_map);
    int_tbl       = interact_tbl (nw, ng, max_interact );
    [gn, gs, gss] = raw_sums(data, g, ng, int_tbl);

    # The Mean squared error is the data - avg for each possible measurement
    sel = select_pat( ones(1,nw), ng, nw);
    SSE= sum( gss(sel) ) - sum( gs(sel).^2 ./ gn(sel) );  
    DFE= sum( gn(sel) -1 );
    MSE= SSE/DFE;
    [SS, DF, MS, F]= factor_sums( gn, gs, gss, [0 0], ng, nw);
    [SS, DF, MS, F]= factor_sums( gn, gs, gss, [0 1], ng, nw); SS,DF
    [SS, DF, MS, F]= factor_sums( gn, gs, gss, [1 0], ng, nw); SS,DF
    [SS, DF, MS, F]= factor_sums( gn, gs, gss, [1 1], ng, nw); SS,DF


  total_mean = mean (group_mean);
  SSB = sum (group_count .* (group_mean - total_mean) .^ 2);
  SST = sumsq (reshape (data, nd, 1) - total_mean);
  SSW = SST - SSB;
  df_b = k - 1;
  df_w = nd - k;
  v_b = SSB / df_b;
  v_w = SSW / df_w;
  f = v_b / v_w;
  pval = 1 - f_cdf (f, df_b, df_w);

  if (nargout == 0)
    printf ('%d-way ANOVA Table:\n\n', nw);
    printf ("Source of Variation   Sum of Squares    df  Empirical Var\n");
    printf ("*********************************************************\n");
    printf ("Between Groups       %15.4f  %4d  %13.4f\n", SSB, df_b, v_b);
    printf ("Within Groups        %15.4f  %4d  %13.4f\n", SSW, df_w, v_w);
    printf ("---------------------------------------------------------\n");
    printf ("Total                %15.4f  %4d\n", SST, nd - 1);
    printf ("\n");
    printf ("Test Statistic f     %15.4f\n", f);
    printf ("p-value              %15.4f\n", pval);
    printf ("\n");
  endif

endfunction

# Create interaction vectors
#
# ng is the number of ANOVA groups
function str=interaction_vectors(varname, ng, max_interact)

# use assoc array to hold sumsqr of each group. Separate groups with _.
# Thus [2,3,4], is added to d.2_3_4_, d.0_3_4_, d.2_0_4_, d.0_0_4_,
#                           d.2_0_3_, d.0_3_0_, d.2_0_0_, d.0_0_0_

    combin= 2^ng;
    static interact_tbl=[];
    if isempty( interact_tbl )
        interact_tbl= zeros( combin, ng);
        idx= (0:combin-1)';
        for i=1:ng;
           interact_tbl(:,i) = ( rem(idx,2^i) >= 2^(i-1) ); 
        end

        # find elements with more than max_interact 1's
        idx = ( sum(interact_tbl') > max_interact );
        interact_tbl(idx,:) =[];
    end

    str= cell{ 2^ng };
    for i= 1:ng

    end

endfunction

function [g,grp_map] = relabel_groups(grps)
    grp_vec= vec(grps);
    s= sort (grp_vec);
    uniq = 1+[0;find(diff(s))];
    # mapping from new grps to old groups
    grp_map = s(uniq);
    # create new group g
    ngroups= length(uniq);
    g= zeros(size(grp_vec));
    for i = 1:ngroups
        g( find( grp_vec== grp_map(i) ) ) = i;
    end
    g= reshape(g, size(grps));
endfunction

# Create interaction table
#
# Input: 
#    nw            number of "ways"
#    ng            number of ANOVA groups
#    max_interact  maximum number of interactions to consider
#                  default is nw
function int_tbl =interact_tbl(nw, ng, max_interact)
    if nargin<3
        max_interact= nw;
    end
    combin= 2^nw;
    interact_tbl= zeros( combin, nw);
    idx= (0:combin-1)';
    for i=1:nw;
       interact_tbl(:,i) = ( rem(idx,2^i) >= 2^(i-1) ); 
    end

    # find elements with more than max_interact 1's
    idx = ( sum(interact_tbl',1) > max_interact );
    int_tbl(idx,:) =[];
    combin= size(interact_tbl,1); # update value

    #scale interact_tbl 
    # use ng+1 to map combinations of groups to integers
    # this would be lots easier with a hash data structure
    int_tbl = interact_tbl .* (ones(combin,1) * (ng+1).^(0:nw-1) );
endfunction 

# Calculate sums for each combination
#
# Input: 
#    g             relabelled grouping matrix
#    ng            number of ANOVA groups
#    max_interact
#
# Output (virtual (ng+1)x(nw) matrices):
#    gn            number of data sums in each group
#    gs            sum of data in each group
#    gss           sumsqr of data in each group
function    [gn, gs, gss] = raw_sums(data, g, ng, int_tbl);
    nw=    size(g,2);
    ndata= size(g,1);
    gn= gs= gss=  zeros((ng+1)^nw, 1);
    for i=1:ndata
        # need offset by one for indexing
        datapt= data(i);
        idx = 1+ int_tbl*g(i,:)';
        gn(idx)  +=1;
        gs(idx)  +=datapt;
        gss(idx) +=datapt^2;
    end
endfunction

# Calcualte the various factor sums
# Input:  
#    gn            number of data sums in each group
#    gs            sum of data in each group
#    gss           sumsqr of data in each group
#    select        binary vector of factor for this "way"?
#    ng            number of ANOVA groups
#    nw            number of ways

function [SS, DF, MS, F]= factor_sums( gn, gs, gss, select, ng, nw);
   if all(select == 0)
       SS= gs(1)^2/gn(1);
       return;
   end

   sub_factor_sum=0;
   ff = find(select);
   lff= length(ff);
   for i= 1:2^lff-1
       remove= find( rem( floor( i * 2.^(-lff+1:0) ), 2) );
       sel1= select;
       sel1( ff( remove ) )=0;
       sub_factor_sum+= factor_sums(gn,gs,gss,sel1,ng,nw);
   end

   sel= select_pat( select, ng, nw);
   SS=  sum( gs(sel).^2 ./ gn(sel) ) - sub_factor_sum;
   DF=  sum( gn(sel) -1 );
   MS=  SS/DF;
endfunction

# Calcualte the various factor sums
# Input:  
#    select        binary vector of factor for this "way"?
#    ng            number of ANOVA groups
#    nw            number of ways
function sel= select_pat( select, ng, nw);
   # if select(i) is zero, remove nonzeros
   # if select(i) is zero, remove zero terms for i
   static field=[];

   if length(select) ~= nw;
       error("length of select must be = nw");
   end
   ng1= ng+1;

   if isempty(field)
       # expand 0:(ng+1)^nw in base ng+1
       field= (0:(ng1)^nw-1)'* ng1.^(-nw+1:0);
       field= rem( floor( field), ng1);
       # select zero or non-zero elements
       field= field>0;
   end
   sel= find( all( field == ones(ng1^nw,1)*select(:)', 2) );
endfunction

data=[7  9  9  8 12 10 ...
      9  8 10 11 13 13 ...
      9 10 10 12 10 12]';
grp = [1,1; 1,1; 1,2; 1,2; 1,3; 1,3;
       2,1; 2,1; 2,2; 2,2; 2,3; 2,3;
       3,1; 3,1; 3,2; 3,2; 3,3; 3,3];

