function results_dotspeeds(record)
%RESULTS_DOTSPEEDS show results of moving dots with different speeds
%
%  RESULTS_DOTSPEEDS(RECORD)
%     it loads the measures field of the record into a global measures
%
% 2022-2023, Alexander Heimel

global measures globalrecord %#ok<GVMIS>
evalin('base','global measures globalrecord');

globalrecord = record;

set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesFontSize', 10);

vecSpeed_pix = record.sStimuli.vecSpeed_pix;

sParams = RH_defaultParameters(  );

indLeft = find([record.sStimuli.vecDirection]==0);
indRight = find([record.sStimuli.vecDirection]==180);

selected_measures = select_measures_by_channel( record.measures, record, 'intIndex');
for m = 1:length(selected_measures)
    measures = selected_measures(m);
    if ~measures.boolResponsive
        % silently skipping, because there are too few responses
        logmsg(['Skipping ' num2str(measures.intIndex) ' because of no response'])
        continue
    end
%     if isempty(measures.lmLeft) && isempty(measures.lmRight)
%         logmsg(['Skipping ' num2str(measures.intIndex) ' because of no fit'])
%         continue
%     end

    figure('Name',['Dots ' num2str(measures.intIndex)],'NumberTitle','off');
    subplot(3,2,1)
    plot(record.sStimuli.vecSpeed_deg(indLeft),measures.vecPeakRate(indLeft),'.-','Color',sParams.clrLeft);
    hold on;
    plot(record.sStimuli.vecSpeed_deg(indRight),measures.vecPeakRate(indRight),'.-','Color',sParams.clrRight);
    set(gca,'xscale','log')
    xlabel('Speed (dps)');
    ylabel('Peak rate (sp/s)');
    legend('Left','Right','Location','Best');
    
    
    subplot(3,2,2) % Peak time vs inverse speed
    hold on;
    vecInvSpeed_pix = 1./record.sStimuli.vecSpeed_pix;
    plot(vecInvSpeed_pix(indLeft),measures.vecPeakTime(indLeft),'o','Color',sParams.clrLeft);
    plot(vecInvSpeed_pix(indLeft),measures.vecOnsetTime(indLeft),'x','Color',sParams.clrLeft);
    plot(vecInvSpeed_pix(indRight),measures.vecPeakTime(indRight),'o','Color',sParams.clrRight);
    plot(vecInvSpeed_pix(indRight),measures.vecOnsetTime(indRight),'x','Color',sParams.clrRight);
    
    if ~isempty(measures.lmLeft) %&& ~isnan(measures.lmLeft.Coefficients.pValue(1))
        % plot full regression line
        b = measures.lmLeft.Coefficients.Estimate(1);
        a = measures.lmLeft.Coefficients.Estimate(2);
        x = vecInvSpeed_pix(indLeft);
        plot(x,a*x + b,'-','Color',sParams.clrLeft)
    end
    if ~isempty(measures.lmLeftFromOnset) %&& ~isnan(measures.lmLeft.Coefficients.pValue(1))
        % plot full regression line
        b = measures.lmLeftFromOnset.Coefficients.Estimate(1);
        a = measures.lmLeftFromOnset.Coefficients.Estimate(2);
        x = vecInvSpeed_pix(indLeft);
        plot(x,a*x + b,'--','Color',sParams.clrLeft)
    end
    if ~isempty(measures.lmRight) %&& ~isnan(measures.lmRight.Coefficients.pValue(1))
        % plot full regression line
        b = measures.lmRight.Coefficients.Estimate(1);
        a = measures.lmRight.Coefficients.Estimate(2);
        x = vecInvSpeed_pix(indRight);
        plot(x,a*x + b,'-','Color',sParams.clrRight)
    end
    if ~isempty(measures.lmRightFromOnset) %&& ~isnan(measures.lmLeft.Coefficients.pValue(1))
        % plot full regression line
        b = measures.lmRightFromOnset.Coefficients.Estimate(1);
        a = measures.lmRightFromOnset.Coefficients.Estimate(2);
        x = vecInvSpeed_pix(indLeft);
        plot(x,a*x + b,'--','Color',sParams.clrRight)
    end
    legend off
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    
    xlabel('1/Speed (spp)');
    ylabel('Peak time (s)');
    
    subplot(3,2,3) % Rate vs Stimulus position
    hold on
    indLeft = 2;
    indRight = indLeft + 6;


   vecStimPos_pix = record.sStimuli.vecStimStartX_pix(indLeft) + measures.cellTime{indLeft}*vecSpeed_pix(indLeft);
   plot(vecStimPos_pix,measures.cellRate{indLeft},'-','Color',sParams.clrLeft);
   vecStimPos_pix = record.sStimuli.vecStimStartX_pix(indRight) - measures.cellTime{indRight}*vecSpeed_pix(indRight);
   plot(vecStimPos_pix,measures.cellRate{indRight},'-','Color',sParams.clrRight);

    
    
    dblStimPosAtResponsePeak_pix =  record.sStimuli.vecStimStartX_pix(indLeft) +  measures.vecPeakTime(indLeft)*vecSpeed_pix(indLeft);
    plot(dblStimPosAtResponsePeak_pix*[1 1],ylim,'-','Color',sParams.clrLeft);
    dblStimPosAtResponsePeak_pix =  record.sStimuli.vecStimStartX_pix(indRight) -  measures.vecPeakTime(indRight)*vecSpeed_pix(indRight);
    plot(dblStimPosAtResponsePeak_pix*[1 1],ylim,'-','Color',sParams.clrRight);
    dblStimPosAtResponseOnset_pix =  record.sStimuli.vecStimStartX_pix(indLeft) +  measures.vecOnsetTime(indLeft)*vecSpeed_pix(indLeft);
    plot(dblStimPosAtResponseOnset_pix*[1 1],ylim,'--','Color',sParams.clrLeft);
    dblStimPosAtResponseOnset_pix =  record.sStimuli.vecStimStartX_pix(indRight) -  measures.vecOnsetTime(indRight)*vecSpeed_pix(indRight);
    plot(dblStimPosAtResponseOnset_pix*[1 1],ylim,'--','Color',sParams.clrRight);

    xlabel('Stim position (pix)');
    ylabel(['Rate stim ' num2str(indLeft) ' (sp/s)']);
    %legend('Left','Right');
    xlim([-1200 1200])
    
    subplot(3,2,4) % Rate vs Stimulus position
    hold on
    indLeft = 5;
    indRight = indLeft + 6;
    vecStimPos_pix = record.sStimuli.vecStimStartX_pix(indLeft) + measures.cellTime{indLeft}*vecSpeed_pix(indLeft);
    plot(vecStimPos_pix,measures.cellRate{indLeft},'-','Color',sParams.clrLeft)
    
    vecStimPos_pix = record.sStimuli.vecStimStartX_pix(indRight) - measures.cellTime{indRight}*vecSpeed_pix(indRight);
    plot(vecStimPos_pix,measures.cellRate{indRight},'-','Color',sParams.clrRight)
    
    dblStimPosAtResponsePeak_pix =  record.sStimuli.vecStimStartX_pix(indLeft) +  measures.vecPeakTime(indLeft)*vecSpeed_pix(indLeft);
    plot(dblStimPosAtResponsePeak_pix*[1 1],ylim,'-','Color',sParams.clrLeft);
    dblStimPosAtResponsePeak_pix =  record.sStimuli.vecStimStartX_pix(indRight) -  measures.vecPeakTime(indRight)*vecSpeed_pix(indRight);
    plot(dblStimPosAtResponsePeak_pix*[1 1],ylim,'-','Color',sParams.clrRight);
    dblStimPosAtResponseOnset_pix =  record.sStimuli.vecStimStartX_pix(indLeft) +  measures.vecOnsetTime(indLeft)*vecSpeed_pix(indLeft);
    plot(dblStimPosAtResponseOnset_pix*[1 1],ylim,'--','Color',sParams.clrLeft);
    dblStimPosAtResponseOnset_pix =  record.sStimuli.vecStimStartX_pix(indRight) -  measures.vecOnsetTime(indRight)*vecSpeed_pix(indRight);
    plot(dblStimPosAtResponseOnset_pix*[1 1],ylim,'--','Color',sParams.clrRight);
    
    xlabel('Stim position (pix)');
    ylabel(['Rate stim ' num2str(indLeft) ' (sp/s)']);
    %legend('Left','Right');
    xlim([-1200 1200])

    subplot(3,2,3) % Corrected cumulative spikes
    hold on
    indLeft = 2;
    indRight = indLeft + 6;

    % Left
    vecStimPos_pix = record.sStimuli.vecStimStartX_pix(indLeft) + measures.cellTime{indLeft}*vecSpeed_pix(indLeft);
    vecCorCurSpike = (1:length(vecStimPos_pix))' - (vecStimPos_pix - vecStimPos_pix(1)) * length(vecStimPos_pix) / (vecStimPos_pix(end)-vecStimPos_pix(1));
    plot(vecStimPos_pix,vecCorCurSpike,'-','Color',sParams.clrLeft);
    indPeakTime = find(measures.cellTime{indLeft}>=measures.vecPeakTime(indLeft),1);
    plot(vecStimPos_pix(indPeakTime),vecCorCurSpike(indPeakTime),'o','Color',sParams.clrLeft,'MarkerFaceColor',sParams.clrLeft);
    indOnsetTime = find(measures.cellTime{indLeft}>=measures.vecOnsetTime(indLeft),1);
    plot(vecStimPos_pix(indOnsetTime),vecCorCurSpike(indOnsetTime),'o','Color',sParams.clrLeft);

    % Right
    vecStimPos_pix = record.sStimuli.vecStimStartX_pix(indRight) - measures.cellTime{indRight}*vecSpeed_pix(indRight);
    vecCorCurSpike = (1:length(vecStimPos_pix))' - (vecStimPos_pix - vecStimPos_pix(1)) * length(vecStimPos_pix) / (vecStimPos_pix(end)-vecStimPos_pix(1));
    plot(vecStimPos_pix,(1:length(vecStimPos_pix))' - (vecStimPos_pix - vecStimPos_pix(1)) * length(vecStimPos_pix) / (vecStimPos_pix(end)-vecStimPos_pix(1)) ,'-','Color',sParams.clrRight);
    indPeakTime = find(measures.cellTime{indRight}>=measures.vecPeakTime(indRight),1);
    plot(vecStimPos_pix(indPeakTime),vecCorCurSpike(indPeakTime),'o','Color',sParams.clrRight,'MarkerFaceColor',sParams.clrRight);
    indOnsetTime = find(measures.cellTime{indRight}>=measures.vecOnsetTime(indRight),1);
    plot(vecStimPos_pix(indOnsetTime),vecCorCurSpike(indOnsetTime),'o','Color',sParams.clrRight);

    plot(measures.dblXRFLeft_pix*[1 1],ylim,'-','Color',sParams.clrLeft);
    plot(measures.dblXRFLeftFromOnset_pix*[1 1],ylim,'--','Color',sParams.clrLeft);
    plot(measures.dblXRFRight_pix*[1 1],ylim,'-','Color',sParams.clrRight);
    plot(measures.dblXRFRightFromOnset_pix*[1 1],ylim,'--','Color',sParams.clrRight);

    xlabel('Stim position (pix)');
    ylabel(['CorCurSpikes stim ' num2str(indLeft) ' (sp/s)']);
    xlim([-1200 1200])


    
    subplot(3,2,4) % Corrected cumulative spikes
    hold on
    indLeft = 3;
    indRight = indLeft + 6;

    % Left
    vecStimPos_pix = record.sStimuli.vecStimStartX_pix(indLeft) + measures.cellTime{indLeft}*vecSpeed_pix(indLeft);
    vecCorCurSpike = (1:length(vecStimPos_pix))' - (vecStimPos_pix - vecStimPos_pix(1)) * length(vecStimPos_pix) / (vecStimPos_pix(end)-vecStimPos_pix(1));
    plot(vecStimPos_pix,vecCorCurSpike,'-','Color',sParams.clrLeft);
    indPeakTime = find(measures.cellTime{indLeft}>=measures.vecPeakTime(indLeft),1);
    plot(vecStimPos_pix(indPeakTime),vecCorCurSpike(indPeakTime),'o','Color',sParams.clrLeft,'MarkerFaceColor',sParams.clrLeft);
    indOnsetTime = find(measures.cellTime{indLeft}>=measures.vecOnsetTime(indLeft),1);
    plot(vecStimPos_pix(indOnsetTime),vecCorCurSpike(indOnsetTime),'o','Color',sParams.clrLeft);

    % Right
    vecStimPos_pix = record.sStimuli.vecStimStartX_pix(indRight) - measures.cellTime{indRight}*vecSpeed_pix(indRight);
    vecCorCurSpike = (1:length(vecStimPos_pix))' - (vecStimPos_pix - vecStimPos_pix(1)) * length(vecStimPos_pix) / (vecStimPos_pix(end)-vecStimPos_pix(1));
    plot(vecStimPos_pix,(1:length(vecStimPos_pix))' - (vecStimPos_pix - vecStimPos_pix(1)) * length(vecStimPos_pix) / (vecStimPos_pix(end)-vecStimPos_pix(1)) ,'-','Color',sParams.clrRight);
    indPeakTime = find(measures.cellTime{indRight}>=measures.vecPeakTime(indRight),1);
    plot(vecStimPos_pix(indPeakTime),vecCorCurSpike(indPeakTime),'o','Color',sParams.clrRight,'MarkerFaceColor',sParams.clrRight);
    indOnsetTime = find(measures.cellTime{indRight}>=measures.vecOnsetTime(indRight),1);
    plot(vecStimPos_pix(indOnsetTime),vecCorCurSpike(indOnsetTime),'o','Color',sParams.clrRight);

    plot(measures.dblXRFLeft_pix*[1 1],ylim,'-','Color',sParams.clrLeft);
    plot(measures.dblXRFLeftFromOnset_pix*[1 1],ylim,'--','Color',sParams.clrLeft);
    plot(measures.dblXRFRight_pix*[1 1],ylim,'-','Color',sParams.clrRight);
    plot(measures.dblXRFRightFromOnset_pix*[1 1],ylim,'--','Color',sParams.clrRight);

    xlabel('Stim position (pix)');
    ylabel(['CorCurSpikes stim ' num2str(indLeft) ' (sp/s)']);
    xlim([-1200 1200])


    subplot(3,2,5) % Corrected cumulative spikes
    hold on
    indLeft = 4;
    indRight = indLeft + 6;

    % Left
    vecStimPos_pix = record.sStimuli.vecStimStartX_pix(indLeft) + measures.cellTime{indLeft}*vecSpeed_pix(indLeft);
    vecCorCurSpike = (1:length(vecStimPos_pix))' - (vecStimPos_pix - vecStimPos_pix(1)) * length(vecStimPos_pix) / (vecStimPos_pix(end)-vecStimPos_pix(1));
    plot(vecStimPos_pix,vecCorCurSpike,'-','Color',sParams.clrLeft);
    indPeakTime = find(measures.cellTime{indLeft}>=measures.vecPeakTime(indLeft),1);
    plot(vecStimPos_pix(indPeakTime),vecCorCurSpike(indPeakTime),'o','Color',sParams.clrLeft,'MarkerFaceColor',sParams.clrLeft);
    indOnsetTime = find(measures.cellTime{indLeft}>=measures.vecOnsetTime(indLeft),1);
    plot(vecStimPos_pix(indOnsetTime),vecCorCurSpike(indOnsetTime),'o','Color',sParams.clrLeft);

    % Right
    vecStimPos_pix = record.sStimuli.vecStimStartX_pix(indRight) - measures.cellTime{indRight}*vecSpeed_pix(indRight);
    vecCorCurSpike = (1:length(vecStimPos_pix))' - (vecStimPos_pix - vecStimPos_pix(1)) * length(vecStimPos_pix) / (vecStimPos_pix(end)-vecStimPos_pix(1));
    plot(vecStimPos_pix,(1:length(vecStimPos_pix))' - (vecStimPos_pix - vecStimPos_pix(1)) * length(vecStimPos_pix) / (vecStimPos_pix(end)-vecStimPos_pix(1)) ,'-','Color',sParams.clrRight);
    indPeakTime = find(measures.cellTime{indRight}>=measures.vecPeakTime(indRight),1);
    plot(vecStimPos_pix(indPeakTime),vecCorCurSpike(indPeakTime),'o','Color',sParams.clrRight,'MarkerFaceColor',sParams.clrRight);
    indOnsetTime = find(measures.cellTime{indRight}>=measures.vecOnsetTime(indRight),1);
    plot(vecStimPos_pix(indOnsetTime),vecCorCurSpike(indOnsetTime),'o','Color',sParams.clrRight);

    plot(measures.dblXRFLeft_pix*[1 1],ylim,'-','Color',sParams.clrLeft);
    plot(measures.dblXRFLeftFromOnset_pix*[1 1],ylim,'--','Color',sParams.clrLeft);
    plot(measures.dblXRFRight_pix*[1 1],ylim,'-','Color',sParams.clrRight);
    plot(measures.dblXRFRightFromOnset_pix*[1 1],ylim,'--','Color',sParams.clrRight);

    xlabel('Stim position (pix)');
    ylabel(['CorCurSpikes stim ' num2str(indLeft) ' (sp/s)']);
    xlim([-1200 1200])

    
%     % Panel XRF
%     subplot(3,2,5);
%     hold on
%     plot(vecInvSpeed_pix(1:6),measures.vecPeakXRF_pix,'.-k');
%     plot(vecInvSpeed_pix(1:6),measures.vecMeanXRF_pix,'.--k');
%     plot(xlim,measures.dblXRFLeft_pix*[1 1],'-','Color',sParams.clrLeft);
%     plot(xlim,measures.dblXRFRight_pix*[1 1],'-','Color',sParams.clrRight);
%     plot(xlim,measures.dblXRFLeftFromOnset_pix*[1 1],'--','Color',sParams.clrLeft);
%     plot(xlim,measures.dblXRFRightFromOnset_pix*[1 1],'--','Color',sParams.clrRight);
%     ylim([-0.5 0.5]*1920)
%     xlabel('1/Speed (spp)');
%     ylabel('x (pix)');
%     legend('Peak','Mean','Location','Best');
    
    % Panel Delta t
    subplot(3,2,6); 
    hold on
    plot(vecInvSpeed_pix(1:6),measures.vecPeakDeltaT,'.-k');
    plot(vecInvSpeed_pix(1:6),measures.vecMeanDeltaT,'.--k');
    plot(xlim,measures.dblDeltaTLeft*[1 1],'-','Color',sParams.clrLeft);
    plot(xlim,measures.dblDeltaTRight*[1 1],'-','Color',sParams.clrRight);
    plot(xlim,measures.dblDeltaTLeftFromOnset*[1 1],'--','Color',sParams.clrLeft);
    plot(xlim,measures.dblDeltaTRightFromOnset*[1 1],'-.','Color',sParams.clrRight);
    xlabel('1/Speed (spp)');
    ylabel('\Deltat (s)');
    ylim([-2 2]);
    drawnow
end % m

measures = record.measures;
logmsg('''measures'' and ''globalrecord'' available in workspace');
