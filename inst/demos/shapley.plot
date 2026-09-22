%!demo
%! ## 1. Explaining one prediction
%!
%! load fisheriris
%! mdl = fitrtree (meas(:,2:4), meas(:,1));
%! s = shapley (mdl, 'QueryPoints', meas(1,2:4), ...
%!              'NumObservationsToSample', 'all');
%!
%! ## One value per predictor, summing to the deviation of the prediction
%! ## from the average prediction
%! s.Shapley
%!
%! ## A bar per predictor, the most important at the top
%! plot (s);

%!demo
%! ## 2. Which predictors matter over many points
%!
%! load fisheriris
%! mdl = fitrtree (meas(:,2:4), meas(:,1));
%! s = shapley (mdl, 'QueryPoints', meas(1:5:150,2:4), ...
%!              'NumObservationsToSample', 'all');
%!
%! ## Over several query points the bars hold the mean of the absolute
%! ## values, so the chart reads as an importance plot rather than as an
%! ## explanation of any one point
%! plot (s);

%!demo
%! ## 3. The lesser predictors gathered into one bar
%!
%! load fisheriris
%! mdl = fitrtree (meas(:,2:4), meas(:,1));
%! s = shapley (mdl, 'QueryPoints', meas(1:5:150,2:4), ...
%!              'NumObservationsToSample', 'all');
%!
%! ## Whatever `NumImportantPredictors` leaves out is summed into one
%! ## further bar, so nothing goes missing from the picture
%! plot (s, 'NumImportantPredictors', 1);

%!demo
%! ## 4. A classifier is explained one class at a time
%!
%! load fisheriris
%! mdl = fitctree (meas, species);
%! s = shapley (mdl, 'QueryPoints', meas([1, 60, 120],:), ...
%!              'NumObservationsToSample', 'all');
%!
%! ## Over several query points every class is drawn, in the model's own
%! ## order, and a legend names them
%! plot (s);
