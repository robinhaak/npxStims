function results_gratingpatches(record)
%RESULTS_GRATINGPATCHES shows results of moving dots with different speeds
%
%  RESULTS_GRATINGPATCHES(RECORD)
%     it loads the measures field of the record into a global measures
%
% 2023, Alexander Heimel

global measures globalrecord %#ok<GVMIS>
evalin('base','global measures globalrecord');

globalrecord = record;

set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesFontSize', 10);

%interpolate
vecX_pix_interp = linspace(record.sStimuli.vecX_pix(1),record.sStimuli.vecX_pix(end),16);
vecY_pix_interp = linspace(record.sStimuli.vecY_pix(1),record.sStimuli.vecY_pix(end),9);


selected_measures = select_measures_by_channel( record.measures, record, 'intIndex');

for m = 1:length(selected_measures)
    measures = selected_measures(m);
    
    figure('Name',['Grating ' num2str(measures.intIndex)],'NumberTitle','off');
    hold on;
    title(['Channel: ' num2str(measures.intIndex)]);
    imagesc(vecX_pix_interp,vecY_pix_interp,measures.matAvgResp);
    set(gca, 'YDir','reverse');
    cb = colorbar;
    cb.Label.String = 'spks/s';
    axis image off
end % m


measures = record.measures;
logmsg('''measures'' and ''globalrecord'' available in workspace');