function results_gratingpatches(record)
%RESULTS_GRATINGPATCHES shows results of moving dots with different speeds
%
%  RESULTS_GRATINGPATCHES(RECORD)
%     it loads the measures field of the record into a global measures
%
% 2023, Alexander Heimel

global measures globalrecord 
evalin('base','global measures globalrecord');

globalrecord = record;

sParams = RH_defaultParameters( );

set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesFontSize', 10);

%interpolate
vecX_pix_interp = linspace(record.sStimuli.vecX_pix(1),record.sStimuli.vecX_pix(end),16);
vecY_pix_interp = linspace(record.sStimuli.vecY_pix(1),record.sStimuli.vecY_pix(end),9);

% ADHOC SHOULD BE READ IN FROM STIMULUS
vecX_pix_interp = vecX_pix_interp - 1920/2;


selected_measures = select_measures_by_channel( record.measures, record, 'intIndex');

for m = 1:length(selected_measures)
    measures = selected_measures(m);
    
    if ~measures.boolResponsive
        continue
    end
    
    
    figure('Name',['Grating ' num2str(measures.intIndex)],'NumberTitle','off');
    
    % Figure receptive field 
    subplot(2,2,1);
    hold on;
    title(['Channel: ' num2str(measures.intIndex)]);
    imagesc(vecX_pix_interp,vecY_pix_interp,measures.matAvgResp);
    colormap gray
    %set(gca, 'YDir','reverse');
    cb = colorbar('location','southoutside');
    cb.Label.String = 'spks/s';
    axis image off
    if isfield(measures,'dblXRFLeftFromOnsetFromMovingDots_pix')
       plot(measures.dblXRFLeftFromOnsetFromMovingDots_pix*[1 1],ylim,'--','Color',sParams.clrLeft);
       plot(measures.dblXRFLeftFromMovingDots_pix*[1 1],ylim,'-','Color',sParams.clrLeft);
       plot(measures.dblXRFRightFromOnsetFromMovingDots_pix*[1 1],ylim,'--','Color',sParams.clrRight);
       plot(measures.dblXRFRightFromMovingDots_pix*[1 1],ylim,'-','Color',sParams.clrRight);
        
    end
    contour(vecX_pix_interp,vecY_pix_interp,measures.matSignificant, [0.5 0.5], 'Color',sParams.clrPatches, 'LineWidth', 2);
    plot(measures.dblXRFLeft_pix*[1 1],ylim,'-','Color',sParams.clrPatches);
    plot(measures.dblXRFRight_pix*[1 1],ylim,'-','Color',sParams.clrPatches);
    
    
    subplot(2,2,3);
    histogram(measures.vecPeakLocationSpikeT,'BinWidth',0.010,'FaceColor',0.7*[1 1 1],'EdgeColor',0.7*[1 1 1]);
    hold on;
    plot(measures.dblPeakTime*[1 1],ylim,'-','color',sParams.clrPatches);
    plot(measures.dblOnsetTime*[1 1],ylim,'--','color',sParams.clrPatches);

    plot(measures.dblDeltaTLeftFromOnsetFromMovingDots*[1 1],ylim,'--','color',sParams.clrLeft);
    plot(measures.dblDeltaTRightFromOnsetFromMovingDots*[1 1],ylim,'-.','color',sParams.clrRight);
    plot(measures.dblDeltaTLeftFromMovingDots*[1 1],ylim,'-','color',sParams.clrLeft);
    plot(measures.dblDeltaTRightFromMovingDots*[1 1],ylim,'-','color',sParams.clrRight);

    axis square 
    box off
    xlabel('Time (s)');
    ylabel('Spikes per bin (all reps)');
    
end % m


measures = record.measures;
logmsg('''measures'' and ''globalrecord'' available in workspace');