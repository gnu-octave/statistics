%!demo
%! ## 1. What the simple model says about the query point
%!
%! load fisheriris
%! mdl = fitrtree (meas(:,2:4), meas(:,1));
%! ex = fit (lime (mdl, 'NumSyntheticData', 2000), meas(1,2:4), 3);
%!
%! ## One bar per predictor the simple model was fitted on, holding its
%! ## weight, so a long bar is a predictor the prediction turns on here
%! plot (ex);

%!demo
%! ## 2. A narrower explanation
%!
%! load fisheriris
%! mdl = fitrtree (meas(:,2:4), meas(:,1));
%!
%! ## Asking for fewer predictors gives a simple model that is easier to
%! ## read and further from the one it explains
%! ex = fit (lime (mdl, 'NumSyntheticData', 2000), meas(1,2:4), 1);
%! plot (ex);

%!demo
%! ## 3. A tree is drawn by what it splits on
%!
%! load fisheriris
%! mdl = fitrtree (meas(:,2:4), meas(:,1));
%! ex = fit (lime (mdl, 'NumSyntheticData', 2000), meas(1,2:4), 2, ...
%!           'SimpleModelType', 'tree');
%!
%! ## A tree has no weights to show, so the bars hold predictor importance
%! ## and the chart is titled after what the simple model is
%! plot (ex);
