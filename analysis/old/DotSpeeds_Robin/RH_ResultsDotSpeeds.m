function RH_ResultsDotSpeeds(record)
%RH_RESULTSDOTSPEEDS show results of moving dots with different speeds
%
%  RH_RESULTSDOTSPEEDS(RECORD)
%     it loads the measures field of the record into a global measures
%
% 2022, Alexander Heimel, Robin Haak

global measures %#ok<GVMIS> 
evalin('base','global measures');

set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesFontSize', 10);

vecSpeed_pix = record.sStimuli.vecSpeed_pix;
%0= rightwards (starting LEFT), 180= leftwards (starting RIGHT)
indLeft = record.sStimuli.stimID(record.sStimuli.vecDirection==0);
indRight = record.sStimuli.stimID(record.sStimuli.vecDirection==180);
if sum(diff(abs(record.sStimuli.vecSpeed_pix(indLeft)))<0)>0 || sum(diff(abs(record.sStimuli.vecSpeed_pix(indRight)))<0)>0
    error('Order of vecSpeed_pix in record.sStimuli is weird! Please check');
end
for m = 1:length(record.measures)
    measures = record.measures(m);

%     if measures.dblZetaP(1)>0.1 % no response to slowest stimulus % also add other direction?
%         continue
%     end

    if ~any(measures.dblZetaP < 0.1)
        continue
    end
        
    figure('Name',['Cluster. ' num2str(measures.intClu)],'NumberTitle','off');
    subplot(3,2,1)
    plot(record.sStimuli.vecSpeed_deg(indRight),measures.vecPeakRate(indRight),'.-r');
    hold on; %peak rate vs speed
    plot(record.sStimuli.vecSpeed_deg(indLeft),measures.vecPeakRate(indLeft),'.-b');
    set(gca,'xscale','log')
    xlabel('Speed (dps)');
    ylabel('Peak rate (sp/s)');
    legend('Right','Left','Location','Best');
    %fixfig;

    subplot(3,2,2) %peak time vs inverse speed
    vecInvSpeed_pix = 1./record.sStimuli.vecSpeed_pix;
    plot(vecInvSpeed_pix(indRight),measures.vecPeakTime(indRight),'.-r');
    hold on;
    plot(vecInvSpeed_pix(indLeft),measures.vecPeakTime(indLeft),'.-b');
    xlabel('1/Speed (spp)');
    ylabel('Peak time (s)');
    %fixfig;
   
    subplot(3,2,3) %rate vs Stimulus position 
    ind1 = 2;
    ind2 = 8;
    vecStimPos_pix = -record.intScreenWidth_pix/2 + measures.cellT{ind1}*vecSpeed_pix(ind1);
    plot( vecStimPos_pix,measures.cellRate{ind1},'-r')
    hold on
    vecStimPos_pix = record.intScreenWidth_pix/2 - measures.cellT{ind2}*vecSpeed_pix(ind2);
    plot(vecStimPos_pix,measures.cellRate{ind2},'-b')
    xlabel('Stim position (pix)');
    ylabel(['Rate stim ' num2str(ind1) ' (sp/s)']);
    %legend('Right','Left');
    %fixfig;
    
    subplot(3,2,4) %rate vs Stimulus position 
    ind1 = 5;
    ind2 = 11;
    vecStimPos_pix = -record.intScreenWidth_pix/2 + measures.cellT{ind1}*vecSpeed_pix(ind1);
    plot(vecStimPos_pix,measures.cellRate{ind1},'-r')
    hold on
    vecStimPos_pix = record.intScreenWidth_pix/2 - measures.cellT{ind2}*vecSpeed_pix(ind2);
    plot(vecStimPos_pix,measures.cellRate{ind2},'-b')
    xlabel('Stim position (pix)');
    ylabel(['Rate stim ' num2str(ind1) ' (sp/s)']);
    %legend('Left','Right');
    %fixfig;
    
    subplot(3,2,5); %RF location vs inverse speed
    hold on
    plot(vecInvSpeed_pix(indRight),measures.vecPeakXRF_pix,'.-k');
    plot(vecInvSpeed_pix(indRight),measures.vecMeanXRF_pix,'.--k');
    ylim([-0.5 0.5]*1920)
    xlabel('1/Speed (spp)');
    ylabel('x (pix)');
    legend('Peak','Mean','Location','Best');
    %fixfig;

    subplot(3,2,6); %DeltaT vs inverse speed
    hold on
    plot(vecInvSpeed_pix(indRight),measures.vecPeakDeltaT,'.-k');
    plot(vecInvSpeed_pix(indRight),measures.vecMeanDeltaT,'.--k');
    %plot(vecInvSpeed_pix(indRight),measures.vecDeltaTAllSpikes,'.--k');
    yline(0,'r');
    xlabel('1/Speed (spp)');
    ylabel('\Deltat (s)');
    ylim([-10 10]);
    %legend('Peak','All spikes','Location','Best');
    %fixfig;

end % m

measures = record.measures;