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

vecSpeed_pix = record.sStimuli.sAllDots.vecSpeed_pix;

for m = 1:length(record.measures)
    measures = record.measures(m);
   
    if measures.dblZetaP(1)>0.1 % no response to slowest stimulus
        continue
    end
        
    figure('Name',['Cluster. ' num2str(measures.intClu)],'NumberTitle','off');
    subplot(3,2,1)
    plot(record.sStimuli.sAllDots.vecSpeed_deg(1:6),measures.vecPeakRate(1:6),'.-r');
    hold on;
    plot(record.sStimuli.sAllDots.vecSpeed_deg(7:12),measures.vecPeakRate(7:12),'.-b');
    set(gca,'xscale','log')
    xlabel('Speed (dps)');
    ylabel('Peak rate (sp/s)');
    legend('Left','Right','Location','Best');

    
    subplot(3,2,2) % Peak time vs inverse speed
    vecInvSpeed_pix = 1./record.sStimuli.sAllDots.vecSpeed_pix;
    plot(vecInvSpeed_pix(1:6),measures.vecPeakTime(1:6),'.-r');
    hold on;
    plot(vecInvSpeed_pix(7:12),measures.vecPeakTime(7:12),'.-b');
    xlabel('1/Speed (spp)');
    ylabel('Peak time (s)');
   
    subplot(3,2,3) % Rate vs Stimulus position 
    ind1 = 2;
    ind2 = 8;
    vecStimPos_pix = -record.intScreenWidth_pix/2 + measures.cellT{ind1}*vecSpeed_pix(ind1);
    plot( vecStimPos_pix,measures.cellRate{ind1},'-r')
    hold on
    vecStimPos_pix = record.intScreenWidth_pix/2 - measures.cellT{ind2}*vecSpeed_pix(ind2);
    plot(vecStimPos_pix,measures.cellRate{ind2},'-b')
    xlabel('Stim position (pix)');
    ylabel(['Rate stim ' num2str(ind1) ' (sp/s)']);
    %legend('Left','Right');
    
    subplot(3,2,4) % Rate vs Stimulus position 
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
    
    subplot(3,2,5);
    hold on
    plot(vecInvSpeed_pix(1:6),measures.vecPeakXRF_pix,'.-k');
    plot(vecInvSpeed_pix(1:6),measures.vecMeanXRF_pix,'.--k');
    ylim([-0.5 0.5]*1920)
    xlabel('1/Speed (spp)');
    ylabel('x (pix)');
    legend('Peak','Mean','Location','Best');
    
    subplot(3,2,6);
    hold on
    plot(vecInvSpeed_pix(1:6),measures.vecPeakDeltaT,'.-k');
    plot(vecInvSpeed_pix(1:6),measures.vecMeanDeltaT,'.--k');
    %plot(vecInvSpeed_pix(1:6),measures.vecDeltaTAllSpikes,'.--k');
    xlabel('1/Speed (spp)');
    ylabel('\Deltat (s)');
    %legend('Peak','All spikes','Location','Best');

end % m

measures = record.measures;