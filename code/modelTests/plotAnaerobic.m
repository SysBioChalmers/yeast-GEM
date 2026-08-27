function [metrics, results] = plotAnaerobic(modelAn, doPlot)
% plotAnaerobic
%   Compares predicted exchange rates of the main fermentation products
%   against the measurements of Sjoberg et al. (2024), and plots each
%   prediction relative to its measured value.
%
%   The measurements are read from
%   data/physiology/exchange_data_anaerobic.tsv, together with the glucose
%   uptake rate they were measured at, which is fixed before solving.
%
%   modelAn     an anaerobic model, as returned by anaerobicModel
%   doPlot      whether to draw the bar plot into the current axes
%               (optional, default true)
%
%   metrics     struct summarising the comparison:
%               meanRelativeError   mean of |predicted-measured|/measured
%               maxRelativeError    largest such deviation
%               fractionWithinError fraction of predictions falling within
%                                   measured +/- stdev
%               nMeasurements       number of measurements compared
%               ammoniumExchange    flux through r_1115
%               ATPase              flux through r_0227
%               ammoniumPerATPase   r_1115 / r_0227, measured to be ~1
%   results     table with one row per measurement: rxnID, rxnName,
%               measured, stdev, unit, predicted, relativeError and
%               withinError
%
%   No single R2 is reported. The measurements mix units, the product
%   rates being in mmol/gDW/h and growth in 1/h, so pooling them into one
%   coefficient of determination would not mean anything. The relative
%   deviation per measurement is comparable across units.
%
%   Sjoberg et al. (2024) Evaluation of enzyme-constrained genome-scale
%   model through metabolic engineering of anaerobic co-production of
%   2,3-butanediol and glycerol by Saccharomyces cerevisiae. Metabolic
%   Engineering 82. doi:10.1016/j.ymben.2024.01.007
%
% Usage: [metrics, results] = plotAnaerobic(modelAn, doPlot)

if nargin < 2
    doPlot = true;
end

funcDir  = fileparts(which(mfilename));
dataFile = fullfile(funcDir,'..','..','data','physiology','exchange_data_anaerobic.tsv');
results  = readtable(dataFile,'FileType','text','Delimiter','\t');

% All rows share one glucose uptake rate; check rather than silently using
% the first of several different values.
qGlc = results.glucoseUptake(1);
if any(results.glucoseUptake ~= qGlc)
    error('The measurements give more than one glucose uptake rate')
end

modelAn = setParam(modelAn,'eq','r_1714',-qGlc);
solution = solveLP(modelAn,1);
if isempty(solution.x)
    error('Model cannot be solved at a glucose uptake rate of %g', qGlc)
end

predicted = zeros(height(results),1);
for i = 1:height(results)
    idx = find(strcmp(modelAn.rxns, results.rxnID{i}));
    if isempty(idx)
        error('Reaction ''%s'' is not in the model', results.rxnID{i})
    end
    predicted(i) = abs(solution.x(idx));
end
results.predicted     = predicted;
results.relativeError = (predicted - results.measured) ./ results.measured;
results.withinError   = abs(predicted - results.measured) <= results.stdev;

metrics.meanRelativeError   = mean(abs(results.relativeError));
metrics.maxRelativeError    = max(abs(results.relativeError));
metrics.fractionWithinError = mean(results.withinError);
metrics.nMeasurements       = height(results);

% Ammonium uptake per unit of plasma membrane ATPase flux. Nitrogen has to
% be taken up against the proton gradient that the ATPase maintains, so the
% two are coupled; the measured ratio is close to 1.
metrics.ammoniumExchange  = solution.x(strcmp(modelAn.rxns,'r_1115'));
metrics.ATPase            = solution.x(strcmp(modelAn.rxns,'r_0227'));
metrics.ammoniumPerATPase = metrics.ammoniumExchange / metrics.ATPase;

if doPlot
    measured = results.measured;
    bar(measured./measured,'FaceAlpha',0.5)
    hold on
    bar(results.predicted./measured,'FaceAlpha',0.5)
    er = errorbar(1:height(results), measured./measured, ...
                  results.stdev./measured, results.stdev./measured);
    er.Color     = [0 0 0];
    er.LineStyle = 'none';
    legend({'data','simulation'})
    ylabel('Relative value')
    xticks(1:height(results))
    xticklabels(strrep(results.rxnName,' exchange',''))
    hold off
end
end
