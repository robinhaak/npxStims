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

%% Find related data
h_db = get_fighandle('Neuropixels database*');
if ~isempty(h_db)
    sUserData = get(h_db,'userdata');
    db = sUserData.db;
    strCrit = ['dataset=' record.dataset ',subject=' record.subject ...
        ',sessionid=' record.sessionid ];
    sRecordDotSpeeds = db(find_record(db,[strCrit ',stimulus=dot_speeds']));
    sRecordFlashingDots = db(find_record(db,[strCrit ',stimulus=flashing_dots']));
    sRecordDotSpeeds = db(find_record(db,[strCrit ',stimulus=dot_diffhist']));
else 
    sRecordDotSpeeds = [];
    sRecordFlashingDots = [];
    sRecordDotSpeeds = [];
end

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

%% Subplots_dot_speeds
function subplots_dot_speeds( measures, record, sParams)
indLeft = find([record.sStimuli.vecDirection]==0);
indRight = find([record.sStimuli.vecDirection]==180);

intNumStimuli = length(measures.vecPeakTime);
sParams.clrScheme = zeros(intNumStimuli,3);
sParams.clrScheme(indLeft,:) =   diag(linspace(0.4,1,length(indLeft))) * repmat(sParams.clrLeft,length(indLeft),1)  ;
sParams.clrScheme(indRight,:) =  diag(linspace(0.4,1,length(indRight))) * repmat(sParams.clrRight,length(indRight),1);
intMarkerSize = 20;


figure('Name',['Dots ' num2str(measures.intIndex)],'NumberTitle','off');

subplot(3,2,1)
hold on;
plot(record.sStimuli.vecSpeed_deg(indLeft),measures.vecPeakRate(indLeft),'-','Color',sParams.clrLeft,'MarkerFaceColor',sParams.clrLeft);
plot(record.sStimuli.vecSpeed_deg(indRight),measures.vecPeakRate(indRight),'-','Color',sParams.clrRight,'MarkerFaceColor',sParams.clrRight);
scatter(record.sStimuli.vecSpeed_deg,measures.vecPeakRate,intMarkerSize,sParams.clrScheme,'filled');

set(gca,'xscale','log')
xlabel('Speed (dps)');
ylabel('Peak rate (sp/s)');
legend('Left','Right','Location','Best');


subplot(3,2,2) % Peak time vs inverse speed
hold on;
vecInvSpeed_pix = 1./record.sStimuli.vecSpeed_pix;
scatter(vecInvSpeed_pix,measures.vecPeakTime,intMarkerSize,sParams.clrScheme,'filled');
scatter(vecInvSpeed_pix,measures.vecOnsetTime,intMarkerSize,sParams.clrScheme);

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


subplot(3,2,3) % Corrected cumulative spikes left
hold on
plot_corcumspikes_vs_time(2,measures,record,sParams)
plot_corcumspikes_vs_time(3,measures,record,sParams)
plot_corcumspikes_vs_time(4,measures,record,sParams)
plot(measures.dblXRFLeft_pix*[1 1],ylim,'-','Color',sParams.clrLeft);
plot(measures.dblXRFLeftFromOnset_pix*[1 1],ylim,'--','Color',sParams.clrLeft);
if isfield(measures,'dblXRFLeftFromGratingPatches_pix')
    plot(measures.dblXRFLeftFromGratingPatches_pix*[1 1],ylim,'-','Color',sParams.clrPatches);
end


subplot(3,2,4) % Corrected cumulative spikes right 
plot_corcumspikes_vs_time(8,measures,record,sParams)
plot_corcumspikes_vs_time(9,measures,record,sParams)
plot_corcumspikes_vs_time(10,measures,record,sParams)
if isfield(measures,'dblXRFRightFromGratingPatches_pix')
    plot(measures.dblXRFRightFromGratingPatches_pix*[1 1],ylim,'-','Color',sParams.clrPatches);
end
plot(measures.dblXRFRight_pix*[1 1],ylim,'-','Color',sParams.clrRight);
plot(measures.dblXRFRightFromOnset_pix*[1 1],ylim,'--','Color',sParams.clrRight);



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
%% Subplots_dot_diffhist
function subplots_dot_diffhist( measures, record, sParams) 

figure('Name',['Dot diffhist ' num2str(measures.intIndex)],'NumberTitle','off');
indResponsive = (measures.vecZetaP < min(0.05,1/length(measures.vecZetaP)));

sMeasuresDotSpeeds = get_related_measures( record, 'stimulus=dot_speeds', measures.intIndex );
sMeasuresFlashingDots = get_related_measures( record, 'stimulus=flashing_dots', measures.intIndex );
if isempty(sMeasuresFlashingDots)
    sMeasuresFlashingDots(1).dblXRFLeft_pix = NaN;
    sMeasuresFlashingDots(1).dblXRFRight_pix = NaN;
end


subplot(3,2,3); % CorCumSpikes
hold on
for ind=1:length(record.sStimuli.vecStimStartX_pix)
    if isempty(measures.cellSpikeTimes{ind})
        continue
    end
    vecStimPos_pix = record.sStimuli.vecStimStartX_pix(ind) + measures.cellSpikeTimes{ind}*record.sStimuli.vecSpeed_pix(ind);
    vecCorCumSpike = (1:length(vecStimPos_pix))' - (vecStimPos_pix - vecStimPos_pix(1)) * length(vecStimPos_pix) / (vecStimPos_pix(end)-vecStimPos_pix(1));
    plot(vecStimPos_pix,vecCorCumSpike,'-');
end
xlabel('Stim position (pix)');
ylabel('CorCumSpikes (sp)');
xlim([-1200 1200])
if ~isempty(sMeasuresDotSpeeds)
    if ~isempty(sMeasuresDotSpeeds.dblXRFLeftFromOnset_pix)
        plot(sMeasuresDotSpeeds.dblXRFLeftFromOnset_pix*[1 1],ylim,'--','Color',sParams.clrLeft);
    end
end
plot(sMeasuresFlashingDots.dblXRFLeft_pix*[1 1],ylim,'-','Color',sParams.clrFlashing);
plot(sMeasuresFlashingDots.dblXRFRight_pix*[1 1],ylim,'-','Color',sParams.clrFlashing);



subplot(3,2,4) % Start vs Onset position
hold on
vecStimPosAtResponseOnset_pix = record.sStimuli.vecStimStartX_pix + measures.vecOnsetTime.*record.sStimuli.vecSpeed_pix;
scatter(record.sStimuli.vecStimStartX_pix(indResponsive),vecStimPosAtResponseOnset_pix(indResponsive),20,'filled');
xlabel('Stim position at start (pix)');
ylabel('Stim position at onset (pix)');
xlim([-1200 1200])
ylim([-1200 1200])
axis square
xyline
if ~isempty(sMeasuresDotSpeeds)
    if ~isempty(sMeasuresDotSpeeds.dblXRFLeftFromOnset_pix)
        plot(sMeasuresDotSpeeds.dblXRFLeftFromOnset_pix*[1 1],ylim,'--','Color',sParams.clrLeft);
    end
end
plot(sMeasuresFlashingDots.dblXRFLeft_pix*[1 1],ylim,'-','Color',sParams.clrFlashing);
plot(sMeasuresFlashingDots.dblXRFRight_pix*[1 1],ylim,'-','Color',sParams.clrFlashing);



end
%% Subplots_flashing_dots
function subplots_flashing_dots( measures, record, sParams) %#ok<INUSD> 

figure('Name',['Flashing dots ' num2str(measures.intIndex)],'NumberTitle','off');

sMeasuresDotSpeeds = get_related_measures( record, 'stimulus=dot_speeds', measures.intIndex );
if isempty(sMeasuresDotSpeeds)
    sMeasuresDotSpeeds.dblXRFLeftFromOnset_pix = NaN;
    sMeasuresDotSpeeds.dblXRFRightFromOnset_pix = NaN;
end

intNumStimuli = length(measures.vecPeakTime);
sParams.clrScheme = parula(intNumStimuli);

subplot(3,2,1); % Stimulus X vs Peak rate
hold on
plot(record.sStimuli.vecStimStartX_pix,measures.vecPeakRate,'.-k');
scatter(record.sStimuli.vecStimStartX_pix(measures.vecResponsive),...
    measures.vecPeakRate(measures.vecResponsive),20,sParams.clrScheme(measures.vecResponsive,:),'filled');
xlabel('Stimulus X (pix)');
ylabel('Peak rate (sp/s)');
plot(measures.dblXRFLeft_pix*[1 1],ylim,'-','Color',sParams.clrFlashing);
plot(measures.dblXRFRight_pix*[1 1],ylim,'-','Color',sParams.clrFlashing);

plot(sMeasuresDotSpeeds.dblXRFLeftFromOnset_pix*[1 1],ylim,'--','Color',sParams.clrLeft);
plot(sMeasuresDotSpeeds.dblXRFRightFromOnset_pix*[1 1],ylim,'--','Color',sParams.clrRight);

subplot(3,2,2); % Stimulus X vs Peak time
hold on
plot(record.sStimuli.vecStimStartX_pix,measures.vecPeakTime,'--k');
% plot(record.sStimuli.vecStimStartX_pix(measures.vecResponsive),...
%     measures.vecPeakTime(measures.vecResponsive),'ok');

scatter(record.sStimuli.vecStimStartX_pix(measures.vecResponsive),...
    measures.vecPeakTime(measures.vecResponsive),20,sParams.clrScheme(measures.vecResponsive,:),'filled');

xlabel('Stimulus X (pix)');
ylabel('Peak time (s)');

subplot(3,2,3); % CorCumSpike vs time 
hold on

[~,ind] = max(measures.vecPeakRate);
vecSpikes = measures.cellSpikeTimes{ind};
vecCorCumSpike = corCumFun( vecSpikes );
plot(vecSpikes,vecCorCumSpike,'-','Color',sParams.clrScheme(ind,:));

% Left edge
ind = find(measures.vecResponsive,1,'first');
vecSpikes = measures.cellSpikeTimes{ind};
vecCorCumSpike = corCumFun( vecSpikes );
plot(vecSpikes,vecCorCumSpike,'-','Color',sParams.clrScheme(ind,:));

% Right edge
ind = find(measures.vecResponsive,1,'last');
vecSpikes = measures.cellSpikeTimes{ind};
vecCorCumSpike = corCumFun( vecSpikes );
plot(vecSpikes,vecCorCumSpike,'-','Color',sParams.clrScheme(ind,:));

xlabel('Time (s)');
ylabel('CorCum (sp)');


end



%% Helper functions


function    plot_rate_vs_time(measures,record,indLeft)
vecSpeed_pix = record.sStimuli.vecSpeed_pix;
vecStimPos_pix = record.sStimuli.vecStimStartX_pix(indLeft) + measures.cellEdges{indLeft}*vecSpeed_pix(indLeft);
histogram(measures.cellSpikeCounts{indLeft},vecStimPos_pix,'FaceColor',0.7*[1 1 1],'EdgeColor',0.7*[1 1 1]);
end

function plot_corcumspikes_vs_time(ind,measures,record,sParams)
vecSpeed_pix = record.sStimuli.vecSpeed_pix;

hold on

if ~isempty(measures.cellSpikeTimes{ind})
    vecStimPos_pix = record.sStimuli.vecStimStartX_pix(ind) + cos(record.sStimuli.vecDirection(ind)/180*pi) * measures.cellSpikeTimes{ind}*vecSpeed_pix(ind);
    vecCorCumSpike = (1:length(vecStimPos_pix))' - (vecStimPos_pix - vecStimPos_pix(1)) * length(vecStimPos_pix) / (vecStimPos_pix(end)-vecStimPos_pix(1));
    plot(vecStimPos_pix,vecCorCumSpike,'-','Color',sParams.clrScheme(ind,:));
    indPeakTime = find(measures.cellSpikeTimes{ind}>=measures.vecPeakTime(ind),1);
    plot(vecStimPos_pix(indPeakTime),vecCorCumSpike(indPeakTime),'o','Color',sParams.clrScheme(ind,:),'MarkerFaceColor',sParams.clrScheme(ind,:));
    indOnsetTime = find(measures.cellSpikeTimes{ind}>=measures.vecOnsetTime(ind),1);
    plot(vecStimPos_pix(indOnsetTime),vecCorCumSpike(indOnsetTime),'o','Color',sParams.clrScheme(ind,:));
end

xlabel('Stim position (pix)');
ylabel(['CorCumSpikes (sp)']);
xlim([-1200 1200])

end

function sRelatedRecord = get_related_record( record, strCrit )

h_db = get_fighandle('Neuropixels database*');
if ~isempty(h_db)
    sUserData = get(h_db,'userdata');
    db = sUserData.db;
    strCrit = ['dataset=' record.dataset ',subject=' record.subject ...
        ',date=' record.date ',' strCrit];
    sRelatedRecord = db(find_record(db,strCrit));
else
    sRelatedRecord = [];
end
end

function sRelatedMeasures = get_related_measures( record, strCrit, vecIndex)
if nargin<3
    vecIndex = [];
end
sRelatedRecord = get_related_record( record, strCrit );
if ~isempty(sRelatedRecord)
    sRelatedMeasures = sRelatedRecord.measures;
    if ~isempty(vecIndex)
        sRelatedMeasures = sRelatedMeasures(ismember([sRelatedMeasures.intIndex],vecIndex));
    end
else
    sRelatedMeasures = [];
end
end