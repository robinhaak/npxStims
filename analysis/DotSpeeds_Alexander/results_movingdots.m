function results_movingdots(record)
%RESULTS_MOVINGDOTS show results of moving dots stimuli
%
%  RESULTS_MOVINGDOTS(RECORD)
%     it loads the measures field of the record into a global measures
%
% 2023, Alexander Heimel

%% Preamble

global measures globalrecord %#ok<GVMIS>
evalin('base','global measures globalrecord');

globalrecord = record;

set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesFontSize', 10);

%%
sParams = RH_defaultParameters(  );

selected_units = select_measures_by_channel( record.measures, record, 'intIndex');
for m = 1:length(selected_units)
    measures = selected_units(m);
    if ~measures.boolResponsive
        % silently skipping, because there are too few responses
        logmsg(['Skipping ' num2str(measures.intIndex) ' because of no response'])
        continue
    end

    switch(record.sStimuli.strStimSet)
        case 'flashing_dots'
            subplots_flashing_dots(measures, record, sParams);
        case 'dot_speeds'
            subplots_dot_speeds( measures, record, sParams);
        case 'dot_diffhist'
            subplots_dot_diffhist( measures, record, sParams);
        otherwise
            logmsg(['Results not implemented for ' record.sStimuli.strStimSet] );
    end

    drawnow
end % m


%% Postamble
measures = record.measures;
logmsg('''measures'' and ''globalrecord'' available in workspace');


end

function subplots_dot_speeds( measures, record, sParams)
indLeft = find([record.sStimuli.vecDirection]==0);
indRight = find([record.sStimuli.vecDirection]==180);

figure('Name',['Dots ' num2str(measures.intIndex)],'NumberTitle','off');

subplot(3,2,1)
hold on;
plot(record.sStimuli.vecSpeed_deg(indLeft),measures.vecPeakRate(indLeft),'o-','Color',sParams.clrLeft,'MarkerFaceColor',sParams.clrLeft);
plot(record.sStimuli.vecSpeed_deg(indRight),measures.vecPeakRate(indRight),'o-','Color',sParams.clrRight,'MarkerFaceColor',sParams.clrRight);
set(gca,'xscale','log')
xlabel('Speed (dps)');
ylabel('Peak rate (sp/s)');
legend('Left','Right','Location','Best');


subplot(3,2,2) % Peak time vs inverse speed
hold on;
vecInvSpeed_pix = 1./record.sStimuli.vecSpeed_pix;
plot(vecInvSpeed_pix(indLeft),measures.vecPeakTime(indLeft),'o','Color',sParams.clrLeft,'MarkerFaceColor',sParams.clrLeft);
plot(vecInvSpeed_pix(indLeft),measures.vecOnsetTime(indLeft),'o','Color',sParams.clrLeft);
plot(vecInvSpeed_pix(indRight),measures.vecPeakTime(indRight),'o','Color',sParams.clrRight,'MarkerFaceColor',sParams.clrRight);
plot(vecInvSpeed_pix(indRight),measures.vecOnsetTime(indRight),'x','Color',sParams.clrRight);

x = linspace(min(vecInvSpeed_pix(indLeft)),max(vecInvSpeed_pix(indLeft)),1000);
if ~isempty(measures.lmLeft) %&& ~isnan(measures.lmLeft.Coefficients.pValue(1))
    % plot full regression line
    b = measures.lmLeft.Coefficients.Estimate(1);
    a = measures.lmLeft.Coefficients.Estimate(2);
    plot(x,a*x + b,'-','Color',sParams.clrLeft)
end
if ~isempty(measures.lmLeftFromOnset) %&& ~isnan(measures.lmLeft.Coefficients.pValue(1))
    % plot full regression line
    b = measures.lmLeftFromOnset.Coefficients.Estimate(1);
    a = measures.lmLeftFromOnset.Coefficients.Estimate(2);
    plot(x,a*x + b,'--','Color',sParams.clrLeft)
end

x = linspace(min(vecInvSpeed_pix(indRight)),max(vecInvSpeed_pix(indRight)),1000);
if ~isempty(measures.lmRight) %&& ~isnan(measures.lmRight.Coefficients.pValue(1))
    % plot full regression line
    b = measures.lmRight.Coefficients.Estimate(1);
    a = measures.lmRight.Coefficients.Estimate(2);
    plot(x,a*x + b,'-','Color',sParams.clrRight)
end
if ~isempty(measures.lmRightFromOnset) %&& ~isnan(measures.lmLeft.Coefficients.pValue(1))
    % plot full regression line
    b = measures.lmRightFromOnset.Coefficients.Estimate(1);
    a = measures.lmRightFromOnset.Coefficients.Estimate(2);
    plot(x,a*x + b,'--','Color',sParams.clrRight)
end

legend off
yl = ylim;
ylim([0.001 yl(2)]);
if isfield(measures,'dblPeakTimeFromGratingPatches')
    plot(xlim,measures.dblPeakTimeFromGratingPatches*[1 1],'-','color',sParams.clrPatches)
end

xlabel('1/Speed (spp)');
ylabel('Peak time (s)');

subplot(3,2,3) % Corrected cumulative spikes
plot_rate_vs_time(measures,record,2); hold on
plot_corcumspikes_vs_time(2,measures,record,sParams)

subplot(3,2,4) % Corrected cumulative spikes
plot_corcumspikes_vs_time(3,measures,record,sParams)

subplot(3,2,5) % Corrected cumulative spikes
plot_corcumspikes_vs_time(4,measures,record,sParams)


% Panel Delta t
subplot(3,2,6);
hold on

if ~isempty(measures.lmLeft)
    plot(vecInvSpeed_pix(indLeft),measures.vecPeakTime(indLeft) - measures.lmLeft.Coefficients.Estimate(2).*vecInvSpeed_pix(indLeft),'o','Color',sParams.clrLeft,'MarkerFaceColor',sParams.clrLeft);
end
if ~isempty(measures.lmLeftFromOnset)
    plot(vecInvSpeed_pix(indLeft),measures.vecOnsetTime(indLeft) - measures.lmLeftFromOnset.Coefficients.Estimate(2).*vecInvSpeed_pix(indLeft),'o','Color',sParams.clrLeft);
end
if ~isempty(measures.lmRight)
    plot(vecInvSpeed_pix(indRight),measures.vecPeakTime(indRight) - measures.lmRight.Coefficients.Estimate(2).*vecInvSpeed_pix(indRight),'o','Color',sParams.clrRight,'MarkerFaceColor',sParams.clrRight);
end
if ~isempty(measures.lmRightFromOnset)
    plot(vecInvSpeed_pix(indRight),measures.vecOnsetTime(indRight) - measures.lmRightFromOnset.Coefficients.Estimate(2).*vecInvSpeed_pix(indRight),'o','Color',sParams.clrRight);
end
plot(xlim,measures.dblDeltaTLeft*[1 1],'-','Color',sParams.clrLeft);
plot(xlim,measures.dblDeltaTRight*[1 1],'-','Color',sParams.clrRight);
plot(xlim,measures.dblDeltaTLeftFromOnset*[1 1],'--','Color',sParams.clrLeft);
plot(xlim,measures.dblDeltaTRightFromOnset*[1 1],'-.','Color',sParams.clrRight);

if isfield(measures,'dblPeakTimeFromGratingPatches')
    plot(xlim,measures.dblPeakTimeFromGratingPatches*[1 1],'-','Color',sParams.clrPatches);
    plot(xlim,measures.dblOnsetTimeFromGratingPatches*[1 1],'--','Color',sParams.clrPatches);

end

xlabel('1/Speed (spp)');
ylabel('\Deltat (s)');
%     ylim([-2 2]);

end


function subplots_dot_diffhist( measures, record, sParams) %#ok<INUSD> 

figure('Name',['Dot diffhist ' num2str(measures.intIndex)],'NumberTitle','off');
indResponsive = (measures.vecZetaP < min(0.05,1/length(measures.vecZetaP)));

subplot(3,2,3);
hold on
for ind=1:length(record.sStimuli.vecStimStartX_pix)
    if isempty(measures.cellSpikeTimes{ind})
        continue
    end
    vecStimPos_pix = record.sStimuli.vecStimStartX_pix(ind) + measures.cellSpikeTimes{ind}*record.sStimuli.vecSpeed_pix(ind);
    vecCorCurSpike = (1:length(vecStimPos_pix))' - (vecStimPos_pix - vecStimPos_pix(1)) * length(vecStimPos_pix) / (vecStimPos_pix(end)-vecStimPos_pix(1));
    plot(vecStimPos_pix,vecCorCurSpike,'-');
end
xlabel('Stim position (pix)');
ylabel('CorCumSpikes (sp)');
xlim([-1200 1200])

subplot(3,2,4)
vecStimPosAtResponseOnset_pix = record.sStimuli.vecStimStartX_pix + measures.vecOnsetTime.*record.sStimuli.vecSpeed_pix;
plot(record.sStimuli.vecStimStartX_pix(indResponsive),vecStimPosAtResponseOnset_pix(indResponsive),'o')
xlabel('Stim position at start (pix)');
ylabel('Stim position at onset (pix)');
xlim([-1200 1200])
ylim([-1200 1200])
axis square
xyline

end

function subplots_flashing_dots( measures, record, sParams) %#ok<INUSD> 

figure('Name',['Flashing dots ' num2str(measures.intIndex)],'NumberTitle','off');

indResponsive = (measures.vecZetaP < min(0.05,1/length(measures.vecZetaP)));

subplot(3,2,1);
hold on
plot(record.sStimuli.vecStimStartX_pix,measures.vecPeakRate,'.-k');
plot(record.sStimuli.vecStimStartX_pix(indResponsive),measures.vecPeakRate(indResponsive),'ok');
xlabel('Stimulus X (pix)');
ylabel('Peak rate (sp/s)');

subplot(3,2,2);
hold on
plot(record.sStimuli.vecStimStartX_pix,measures.vecPeakTime,'--k');
plot(record.sStimuli.vecStimStartX_pix(indResponsive),measures.vecPeakTime(indResponsive),'ok');
xlabel('Stimulus X (pix)');
ylabel('Peak time (s)');


end



%% Helper functions

function    plot_rate_vs_time(measures,record,indLeft)
vecSpeed_pix = record.sStimuli.vecSpeed_pix;
vecStimPos_pix = record.sStimuli.vecStimStartX_pix(indLeft) + measures.cellEdges{indLeft}*vecSpeed_pix(indLeft);
histogram(measures.cellSpikeCounts{indLeft},vecStimPos_pix,'FaceColor',0.7*[1 1 1],'EdgeColor',0.7*[1 1 1]);
end

function plot_corcumspikes_vs_time(indLeft,measures,record,sParams)
vecSpeed_pix = record.sStimuli.vecSpeed_pix;

hold on
indRight = indLeft + 6;

% Left
if ~isempty(measures.cellSpikeTimes{indLeft})
    vecStimPos_pix = record.sStimuli.vecStimStartX_pix(indLeft) + measures.cellSpikeTimes{indLeft}*vecSpeed_pix(indLeft);
    vecCorCurSpike = (1:length(vecStimPos_pix))' - (vecStimPos_pix - vecStimPos_pix(1)) * length(vecStimPos_pix) / (vecStimPos_pix(end)-vecStimPos_pix(1));
    plot(vecStimPos_pix,vecCorCurSpike,'-','Color',sParams.clrLeft);
    indPeakTime = find(measures.cellSpikeTimes{indLeft}>=measures.vecPeakTime(indLeft),1);
    plot(vecStimPos_pix(indPeakTime),vecCorCurSpike(indPeakTime),'o','Color',sParams.clrLeft,'MarkerFaceColor',sParams.clrLeft);
    indOnsetTime = find(measures.cellSpikeTimes{indLeft}>=measures.vecOnsetTime(indLeft),1);
    plot(vecStimPos_pix(indOnsetTime),vecCorCurSpike(indOnsetTime),'o','Color',sParams.clrLeft);
end

% Right
if ~isempty(measures.cellSpikeTimes{indRight})
    vecStimPos_pix = record.sStimuli.vecStimStartX_pix(indRight) - measures.cellSpikeTimes{indRight}*vecSpeed_pix(indRight);
    vecCorCurSpike = (1:length(vecStimPos_pix))' - (vecStimPos_pix - vecStimPos_pix(1)) * length(vecStimPos_pix) / (vecStimPos_pix(end)-vecStimPos_pix(1));
    plot(vecStimPos_pix,(1:length(vecStimPos_pix))' - (vecStimPos_pix - vecStimPos_pix(1)) * length(vecStimPos_pix) / (vecStimPos_pix(end)-vecStimPos_pix(1)) ,'-','Color',sParams.clrRight);
    indPeakTime = find(measures.cellSpikeTimes{indRight}>=measures.vecPeakTime(indRight),1);
    plot(vecStimPos_pix(indPeakTime),vecCorCurSpike(indPeakTime),'o','Color',sParams.clrRight,'MarkerFaceColor',sParams.clrRight);
    indOnsetTime = find(measures.cellSpikeTimes{indRight}>=measures.vecOnsetTime(indRight),1);
    plot(vecStimPos_pix(indOnsetTime),vecCorCurSpike(indOnsetTime),'o','Color',sParams.clrRight);
end

drawnow
yl = ylim;

plot(measures.dblXRFLeft_pix*[1 1],ylim,'-','Color',sParams.clrLeft);
plot(measures.dblXRFLeftFromOnset_pix*[1 1],ylim,'--','Color',sParams.clrLeft);
plot(measures.dblXRFRight_pix*[1 1],ylim,'-','Color',sParams.clrRight);
plot(measures.dblXRFRightFromOnset_pix*[1 1],ylim,'--','Color',sParams.clrRight);

if isfield(measures,'dblXRFLeftFromGratingPatches_pix')
    plot(measures.dblXRFLeftFromGratingPatches_pix*[1 1],ylim,'-','Color',sParams.clrPatches);
    plot(measures.dblXRFRightFromGratingPatches_pix*[1 1],ylim,'-','Color',sParams.clrPatches);
end

ylim(yl);
xlabel('Stim position (pix)');
ylabel(['CorCumSpikes stim ' num2str(indLeft) ' (sp/s)']);
xlim([-1200 1200])

end

