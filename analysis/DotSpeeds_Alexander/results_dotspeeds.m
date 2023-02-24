function results_dotspeeds(record)
%RESULTS_DOTSPEEDS show results of moving dots with different speeds
%
%  RESULTS_DOTSPEEDS(RECORD)
%     it loads the measures field of the record into a global measures
%
% 2022-2023, Alexander Heimel

<<<<<<< HEAD
global measures %#ok<GVMIS>
evalin('base','global measures');
=======
global measures globalrecord %#ok<GVMIS>
evalin('base','global measures globalrecord');

globalrecord = record;
>>>>>>> 4f50c0e1182137ef01f3280bd5dbd8de9686c074

set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesFontSize', 10);

vecSpeed_pix = record.sStimuli.vecSpeed_pix;

<<<<<<< HEAD
for m = 1:length(record.measures)
    measures = record.measures(m);


    if measures.dblZetaP(1)>0.1 % no response to slowest stimulus
        continue
    end

    figure('Name',['Ch. ' num2str(measures.intCh)],'NumberTitle','off');
=======

indLeft = find([record.sStimuli.vecDirection]==0);
indRight = find([record.sStimuli.vecDirection]==180);

clrLeft = [1 0 0];
clrRight = [0 0 1];

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
>>>>>>> 4f50c0e1182137ef01f3280bd5dbd8de9686c074
    subplot(3,2,1)
    plot(record.sStimuli.vecSpeed_deg(indLeft),measures.vecPeakRate(indLeft),'.-','Color',clrLeft);
    hold on;
    plot(record.sStimuli.vecSpeed_deg(indRight),measures.vecPeakRate(indRight),'.-','Color',clrRight);
    set(gca,'xscale','log')
    xlabel('Speed (dps)');
    ylabel('Peak rate (sp/s)');
    legend('Left','Right','Location','Best');
<<<<<<< HEAD


=======
    
    
>>>>>>> 4f50c0e1182137ef01f3280bd5dbd8de9686c074
    subplot(3,2,2) % Peak time vs inverse speed
    vecInvSpeed_pix = 1./record.sStimuli.vecSpeed_pix;
    plot(vecInvSpeed_pix(indLeft),measures.vecPeakTime(indLeft),'o','Color',clrLeft);
    hold on;
<<<<<<< HEAD
    plot(vecInvSpeed_pix(7:12),measures.vecPeakTime(7:12),'.-b');
    set(gca,'xscale','log')
    xlabel('1/Speed (spp)');
    ylabel('Peak time (s)');

    subplot(3,2,3) % Rate vs Stimulus position
=======
    plot(vecInvSpeed_pix(indRight),measures.vecPeakTime(indRight),'o','Color',clrRight);
    
    if ~isempty(measures.lmLeft) %&& ~isnan(measures.lmLeft.Coefficients.pValue(1))
        %h = plot(measures.lmLeft);
        %set(h(1),'color',clrLeft);
        %set(h(2),'color',clrLeft);
        %set(h(3),'color',clrLeft);

        % plot full regression line
        b = measures.lmLeft.Coefficients.Estimate(1);
        a = measures.lmLeft.Coefficients.Estimate(2);
        x = vecInvSpeed_pix(indLeft);
        plot(x,a*x + b,'-','Color',clrLeft)
    end
    if ~isempty(measures.lmRight) %&& ~isnan(measures.lmRight.Coefficients.pValue(1))
        %h = plot(measures.lmRight);
        %set(h(1),'color',clrRight)
        %set(h(2),'color',clrRight)
        %set(h(3),'color',clrRight)

                % plot full regression line
        b = measures.lmRight.Coefficients.Estimate(1);
        a = measures.lmRight.Coefficients.Estimate(2);
        x = vecInvSpeed_pix(indRight);
        plot(x,a*x + b,'-','Color',clrRight)

    end
    legend off
    xlabel('1/Speed (spp)');
    ylabel('Peak time (s)');
    
    subplot(3,2,3) % Rate vs Stimulus position
    hold on
>>>>>>> 4f50c0e1182137ef01f3280bd5dbd8de9686c074
    ind1 = 2;
    ind2 = ind1 + 6;

    vecStimPos_pix = record.sStimuli.vecStimStartX_pix(ind1) + measures.cellTime{ind1}*vecSpeed_pix(ind1);
    plot(vecStimPos_pix,measures.cellRate{ind1},'-','Color',clrLeft);
    vecStimPos_pix = record.sStimuli.vecStimStartX_pix(ind2) - measures.cellTime{ind2}*vecSpeed_pix(ind2);
    plot(vecStimPos_pix,measures.cellRate{ind2},'-','Color',clrRight);

    
    dblStimPosAtResponsePeak_pix =  record.sStimuli.vecStimStartX_pix(ind1) +  measures.vecPeakTime(ind1)*vecSpeed_pix(ind1);
    plot(dblStimPosAtResponsePeak_pix*[1 1],ylim,'-','Color',clrLeft);
    dblStimPosAtResponsePeak_pix =  record.sStimuli.vecStimStartX_pix(ind2) -  measures.vecPeakTime(ind2)*vecSpeed_pix(ind2);
    plot(dblStimPosAtResponsePeak_pix*[1 1],ylim,'-','Color',clrRight);
    dblStimPosAtResponseOnset_pix =  record.sStimuli.vecStimStartX_pix(ind1) +  measures.vecOnsetTime(ind1)*vecSpeed_pix(ind1);
    plot(dblStimPosAtResponseOnset_pix*[1 1],ylim,'--','Color',clrLeft);
    dblStimPosAtResponseOnset_pix =  record.sStimuli.vecStimStartX_pix(ind2) -  measures.vecOnsetTime(ind2)*vecSpeed_pix(ind2);
    plot(dblStimPosAtResponseOnset_pix*[1 1],ylim,'--','Color',clrRight);

    xlabel('Stim position (pix)');
    ylabel(['Rate stim ' num2str(ind1) ' (sp/s)']);
    %legend('Left','Right');
<<<<<<< HEAD

    subplot(3,2,4) % Rate vs Stimulus position
=======
    
    subplot(3,2,4) % Rate vs Stimulus position
    hold on
>>>>>>> 4f50c0e1182137ef01f3280bd5dbd8de9686c074
    ind1 = 5;
    ind2 = ind1 + 6;
    vecStimPos_pix = record.sStimuli.vecStimStartX_pix(ind1) + measures.cellTime{ind1}*vecSpeed_pix(ind1);
    plot(vecStimPos_pix,measures.cellRate{ind1},'-','Color',clrLeft)
    
    vecStimPos_pix = record.sStimuli.vecStimStartX_pix(ind2) - measures.cellTime{ind2}*vecSpeed_pix(ind2);
    plot(vecStimPos_pix,measures.cellRate{ind2},'-','Color',clrRight)
    
    dblStimPosAtResponsePeak_pix =  record.sStimuli.vecStimStartX_pix(ind1) +  measures.vecPeakTime(ind1)*vecSpeed_pix(ind1);
    plot(dblStimPosAtResponsePeak_pix*[1 1],ylim,'-','Color',clrLeft);
    dblStimPosAtResponsePeak_pix =  record.sStimuli.vecStimStartX_pix(ind2) -  measures.vecPeakTime(ind2)*vecSpeed_pix(ind2);
    plot(dblStimPosAtResponsePeak_pix*[1 1],ylim,'-','Color',clrRight);
    dblStimPosAtResponseOnset_pix =  record.sStimuli.vecStimStartX_pix(ind1) +  measures.vecOnsetTime(ind1)*vecSpeed_pix(ind1);
    plot(dblStimPosAtResponseOnset_pix*[1 1],ylim,'--','Color',clrLeft);
    dblStimPosAtResponseOnset_pix =  record.sStimuli.vecStimStartX_pix(ind2) -  measures.vecOnsetTime(ind2)*vecSpeed_pix(ind2);
    plot(dblStimPosAtResponseOnset_pix*[1 1],ylim,'--','Color',clrRight);
    
    xlabel('Stim position (pix)');
    ylabel(['Rate stim ' num2str(ind1) ' (sp/s)']);
    xlim([-record.intScreenWidth_pix/2 record.intScreenWidth_pix/2]);
    %legend('Left','Right');
<<<<<<< HEAD

=======
    
    % Panel XRF
>>>>>>> 4f50c0e1182137ef01f3280bd5dbd8de9686c074
    subplot(3,2,5);
    hold on
    plot(vecInvSpeed_pix(1:6),measures.vecPeakXRF_pix,'.-k');
    plot(vecInvSpeed_pix(1:6),measures.vecMeanXRF_pix,'.--k');
    plot(xlim,measures.dblXRFLeft_pix*[1 1],'-','Color',clrLeft);
    plot(xlim,measures.dblXRFRight_pix*[1 1],'-','Color',clrRight);
    plot(xlim,measures.dblXRFLeftFromOnset_pix*[1 1],'--','Color',clrLeft);
    plot(xlim,measures.dblXRFRightFromOnset_pix*[1 1],'--','Color',clrRight);
    ylim([-0.5 0.5]*1920)
    xlabel('1/Speed (spp)');
    ylabel('x (pix)');
    set(gca,'xscale','log')

    legend('Peak','Mean','Location','Best');
<<<<<<< HEAD

    subplot(3,2,6);
=======
    
    % Panel Delta t
    subplot(3,2,6); 
>>>>>>> 4f50c0e1182137ef01f3280bd5dbd8de9686c074
    hold on
    plot(vecInvSpeed_pix(1:6),measures.vecPeakDeltaT,'.-k');
    plot(vecInvSpeed_pix(1:6),measures.vecMeanDeltaT,'.--k');
    plot(xlim,measures.dblDeltaTLeft*[1 1],'-','Color',clrLeft);
    plot(xlim,measures.dblDeltaTRight*[1 1],'-','Color',clrRight);
    plot(xlim,measures.dblDeltaTLeftFromOnset*[1 1],'--','Color',clrLeft);
    plot(xlim,measures.dblDeltaTRightFromOnset*[1 1],'--','Color',clrRight);
    xlabel('1/Speed (spp)');
    ylabel('\Deltat (s)');
<<<<<<< HEAD
    set(gca,'xscale','log')

    %legend('Peak','All spikes','Location','Best');

=======
    ylim([-2 2]);
    
>>>>>>> 4f50c0e1182137ef01f3280bd5dbd8de9686c074
end % m

measures = record.measures;
logmsg('''measures'' and ''globalrecord'' available in workspace');
