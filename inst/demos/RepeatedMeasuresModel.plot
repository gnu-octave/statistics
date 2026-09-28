%!demo
%! ## 1. Every subject's measurements
%!
%! load fisheriris
%! t = table (species, meas(:,1), meas(:,2), meas(:,3), meas(:,4), ...
%!            'VariableNames', {'species', 'meas1', 'meas2', 'meas3', 'meas4'});
%! rm = fitrm (t, 'meas1-meas4 ~ species');
%!
%! ## One line per flower through its four measurements
%! plot (rm);

%!demo
%! ## 2. Lines coloured by group
%!
%! load fisheriris
%! t = table (species, meas(:,1), meas(:,2), meas(:,3), meas(:,4), ...
%!            'VariableNames', {'species', 'meas1', 'meas2', 'meas3', 'meas4'});
%! rm = fitrm (t, 'meas1-meas4 ~ species');
%!
%! ## `Group` colours the lines by species, so a profile that one species
%! ## follows and the others do not stands out
%! plot (rm, 'Group', 'species', 'Marker', 'none');
