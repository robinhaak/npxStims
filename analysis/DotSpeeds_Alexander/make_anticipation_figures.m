%MAKE_ANTICIPATION_FIGURES
%
% Makes some figures for the predictive motion project
%
% 2022, Alexander Heimel


%% Load database
sParams = RH_defaultParameters();
dbfilename = fullfile(sParams.strOutputPath,'Anticipation','Data_analysis','db_anticipation.mat');
load(dbfilename,'db')
logmsg('Loaded datebase as ''db''. Access with experiment_db(db).')
return

%% Combine moving dots and grating patches results
for i=2:length(db)
   db(i) = analyse_add_dots_to_patches( db(i), db, true);
end    


%% Record template
record.lab = 'Heimellab'; %#ok<UNRCH> 
record.project = 'Anticipation';
record.dataset = '22.35.03';
record.subject = '83554';
record.condition = '';
record.stimulus = 'MovingDots';
record.setup = 'Neuropixels';
record.investigator = 'RobinHaak';
record.date = '20221122';
record.logfile = '';
record.sessionid = '20221122_83554_15_55_45_RH_MovingDots';
record.sessnr = NaN;
record.version = '1.0';
record.path = '';
record.datatype = 'neuropixels'
record.analysisfunction = 'analyse_dotspeeds';
record.resultsfunction = 'results_dotspeeds'
record.measures = [];
db = record;

%% Reanalyse all records
for i=1:length(db)
    db(i) = analyse_dotspeeds(db(i)); 
end
strPath = fullfile(fileparts(which('make_anticipation_figures')),'..','..');
sParams.strOutputPath = strPath; % should be loaded from parameter file instead
dbfilename = fullfile(sParams.strOutputPath,'Anticipation','Data_analysis','db_anticipation.mat');
save(dbfilename,'db');

%% Make figure speed tuning
record = db(1);
figure
hold on

fld = 'vecPeakRate';
indLeft = find(record.sStimuli.vecDirection==0);
indRight = find(record.sStimuli.vecDirection==180);
matVal = cat(1,record.measures.(fld));
plot(record.sStimuli.vecSpeed_deg(indLeft),mean(matVal(:,indLeft),1,'OmitNan'),'.-r'     );
plot(record.sStimuli.vecSpeed_deg(indRight),mean(matVal(:,indRight),1,'OmitNan'),'.-b'     );
errorbar(record.sStimuli.vecSpeed_deg(indLeft),mean(matVal(:,indLeft),1,'OmitNan'),sem(matVal(:,indLeft),1),'.r'); % sem from invivotools
errorbar(record.sStimuli.vecSpeed_deg(indRight),mean(matVal(:,indRight),1,'OmitNan'),sem(matVal(:,indRight),1),'.b');
ylabel(fld)
legend('Left','Right')
xlabel('Speed (deg/s)')
title('Population')

%% Make figure deltaT versus speed
record = db(1);
figure
hold on

fld = 'vecPeakDeltaT';
vecSpeed_pix = record.sStimuli.vecSpeed_pix;
ind = find(record.sStimuli.vecDirection==0);
matVal = cat(1,record.measures.(fld));
plot(record.sStimuli.vecSpeed_deg(ind),mean(matVal(:,ind),1,'OmitNan'),'.-k'     );
errorbar(record.sStimuli.vecSpeed_deg(ind),mean(matVal(:,ind),1,'OmitNan'),sem(matVal(:,ind),1),'.k');

fld = 'vecMeanDeltaT';
ind = find(record.sStimuli.vecDirection==0);
matVal = cat(1,record.measures.(fld));
plot(record.sStimuli.vecSpeed_deg(ind),mean(matVal(:,ind),1,'OmitNan'),'.--k'     );
errorbar(record.sStimuli.vecSpeed_deg(ind),mean(matVal(:,ind),1,'OmitNan'),sem(matVal(:,ind),1),'.k');

legend('Peak','Mean','Location','Best');

ylabel(fld)
xlabel('Speed (deg/s)')
title('Population')

%% Get anticipatory indices

% some responsive channels
% 20221122_83554_15_55_45_RH_MovingDots: 142,145,147,148,149,150:155,157,167
% 20230117_85235_12_15_45_RH_MovingDots: 198,205,212-228,229,230,
% 20230117_singleUnits: 121, 154-156, 184-218
record = db(find_record(db,'subject=85235,sessionid=20230117_85235_12_15_45_RH_MovingDots'));

for i=1:length(db)
    disp(['Record ' recordfilter(db(i))]);
    ind=find([db(i).measures.dblDeltaTLeft]<0 | [db(i).measures.dblDeltaTRight]<0 );
    logmsg(['Possible anticipatory responses in : ' mat2str([db(i).measures(ind).intIndex])]);
end

%% Figure Delta T versus Depth
%record = db(find_record(db,'subject=85235,sessionid=20230117_85235_12_15_45_RH_MovingDots'));
%record = db(find_record(db,'subject=85235,sessionid=20230117_singleUnits_new'));
record = db(find_record(db,'subject=85235,sessionid=20230118_singleUnits,stimulus=MovingDots'));
%record = db(find_record(db,'subject=85235,sessionid=20230119_singleUnits,stimulus=MovingDots'));
%record = db(find_record(db,'subject=85234,sessionid=20230118_singleUnits,stimulus=MovingDots'));
%record = db(find_record(db,'subject=85234,sessionid=20230119_singleUnits,stimulus=MovingDots'));

vecDepth = [record.measures.dblDepth_um]; 
vecDeltaTLeft = [record.measures.dblDeltaTLeft];
vecDeltaTRight = [record.measures.dblDeltaTRight];
vecDeltaTLeftFromOnset = [record.measures.dblDeltaTLeftFromOnset];
vecDeltaTRightFromOnset = [record.measures.dblDeltaTRightFromOnset];
clrLeft = [1 0 0];
clrRight = [0 0 1];

figure;
hold on
plot(vecDeltaTLeft,vecDepth,'o','Color',clrLeft)
plot(vecDeltaTRight,vecDepth,'o','Color',clrRight)
plot(vecDeltaTLeftFromOnset,vecDepth,'x','Color',clrLeft)
plot(vecDeltaTRightFromOnset,vecDepth,'x','Color',clrRight)
plot([0 0],ylim,'-k');
%ylim([1400 2600])
xlim([-1.5 2]);
set(gca,'YDir','reverse');
legend('Left by peak time','Right by peak time','Left by onset time','Right by onset time','location','northeastoutside')
legend box off
xlabel('\DeltaT (s)')
ylabel('Depth (\mum)');
title(recordfilter(record))

%% Figure XRF versus Depth

record = db(find_record(db,'subject=85235,sessionid=20230117_85235_12_15_45_RH_MovingDots'));
record = db(find_record(db,'subject=85235,sessionid=20230117_singleUnits'));

vecDepth = [record.measures.dblDepth_um]; 
vecXRFLeft_pix = [record.measures.dblXRFLeft_pix];
vecXRFRight_pix = [record.measures.dblXRFRight_pix];
vecXRFLeftFromOnset_pix = [record.measures.dblXRFLeftFromOnset_pix];
vecXRFRightFromOnset_pix = [record.measures.dblXRFRightFromOnset_pix];
clrLeft = [1 0 0];
clrRight = [0 0 1];

figure;
hold on
plot(vecXRFLeft_pix,vecDepth,'o','Color',clrLeft)
plot(vecXRFRight_pix,vecDepth,'o','Color',clrRight)
plot(vecXRFLeftFromOnset_pix,vecDepth,'x','Color',clrLeft)
plot(vecXRFRightFromOnset_pix,vecDepth,'x','Color',clrRight)
ylim([1400 2600])
xlim([-1 1]*0.5*1920);
legend('Left by peak time','Right by peak time','Left by onset time','Right by onset time','location','northeastoutside')
legend box off
xlabel('X_{RF} (pix')
ylabel('Depth from Kilosort (\mum)');
set(gca,'YDir','reverse');

%% Figure Patch peak times versus Depth
record = db(find_record(db,'subject=85235,sessionid=20230117_85235_13_15_16_RH_GratingPatches,date=20230117'));
record = db(find_record(db,'subject=85235,sessionid=20230117_singleUnits,stimulus=GratingPatches'));
record = db(find_record(db,'subject=85235,sessionid=20230118_singleUnits,stimulus=GratingPatches'));
%record = db(find_record(db,'subject=85235,sessionid=20230119_singleUnits,stimulus=GratingPatches'));

indResponsive = [record.measures.boolResponsive]; 
vecDepth = [record.measures.dblDepth_um]; 
vecPeakTime = [record.measures.dblPeakTime]; 
vecOnsetTime = [record.measures.dblOnsetTime]; 

figure;
hold on

plot(vecPeakTime(indResponsive),vecDepth(indResponsive),'o')
plot(vecOnsetTime(indResponsive),vecDepth(indResponsive),'x')
xlabel('Peak time (s)')
ylabel('Depth (\mum)');
set(gca,'YDir','reverse');
axis square
title(recordfilter(record))

%% Figure onset patch versus onset moving dots
cellSessions = {...
    'subject=85234,sessionid=20230117_singleUnits',...
    'subject=85234,sessionid=20230118_singleUnits',...
    'subject=85234,sessionid=20230119_singleUnits',...
    'subject=85235,sessionid=20230117_singleUnits',...
    'subject=85235,sessionid=20230118_singleUnits',...
    'subject=85235,sessionid=20230119_singleUnits',...
    }

vecAllDeltaTLeft = [];
vecAllDeltaTRight = [];
vecAllDeltaTLeftFromOnset = [];
vecAllDeltaTRightFromOnset = [];
vecAllPeakTimeGratingPatches = [];
vecAllOnsetTimeGratingPatches = [];
vecAllSessionIndex = [];
vecAllIndex = [];

for s = 1:length(cellSessions)
    record = db(find_record(db,[cellSessions{s} ',stimulus=MovingDots']));
    indResponsive = [record.measures.boolResponsive];
    vecDeltaTLeft = [record.measures.dblDeltaTLeft];
    vecDeltaTRight = [record.measures.dblDeltaTRight];
    vecDeltaTLeftFromOnset = [record.measures.dblDeltaTLeftFromOnset];
    vecDeltaTRightFromOnset = [record.measures.dblDeltaTRightFromOnset];
    vecIndex = [record.measures.intIndex];
    vecDeltaTLeft(~indResponsive) = NaN;
    vecDeltaTRight(~indResponsive) = NaN;
    vecDeltaTLeftFromOnset(~indResponsive) = NaN;
    vecDeltaTRightFromOnset(~indResponsive) = NaN;
    
    vecAllDeltaTLeft =[vecAllDeltaTLeft vecDeltaTLeft];
    vecAllDeltaTRight =[vecAllDeltaTRight vecDeltaTRight];
    vecAllDeltaTLeftFromOnset =[vecAllDeltaTLeftFromOnset vecDeltaTLeftFromOnset];
    vecAllDeltaTRightFromOnset =[vecAllDeltaTRightFromOnset vecDeltaTRightFromOnset];
    vecAllSessionIndex = [vecAllSessionIndex s*ones(1,length(vecDeltaTRightFromOnset))];
    vecAllIndex = [vecAllIndex vecIndex];
    
    record = db(find_record(db,[cellSessions{s} ',stimulus=GratingPatches']));
    vecPeakTimeGratingPatches = NaN(size(vecDeltaTLeft));
    vecOnsetTimeGratingPatches = NaN(size(vecDeltaTLeft));
    for i = 1:length(indResponsive)
        ind = find([record.measures.intIndex]==vecIndex(i));
        if isempty(ind)
            continue
        end
        if ~record.measures(ind).boolResponsive
            continue
        end
        vecPeakTimeGratingPatches(i) = record.measures(ind).dblPeakTime;
        vecOnsetTimeGratingPatches(i) = record.measures(ind).dblOnsetTime;
    end % i
    
    vecAllPeakTimeGratingPatches = [vecAllPeakTimeGratingPatches vecPeakTimeGratingPatches ];
    vecAllOnsetTimeGratingPatches = [vecAllOnsetTimeGratingPatches vecOnsetTimeGratingPatches ];

    
end % s

figure
hold on
plot(vecAllPeakTimeGratingPatches,vecAllDeltaTRight,'ob' );
plot(vecAllPeakTimeGratingPatches,vecAllDeltaTLeft,'or' );
plot(vecAllOnsetTimeGratingPatches,vecAllDeltaTRightFromOnset,'xb' );
plot(vecAllOnsetTimeGratingPatches,vecAllDeltaTLeftFromOnset,'xr' );
xlim([-0.5 0.5]);
ylim([-0.5 0.5]);
xyline;
axis square
xlabel('Latency from patches (s)');
ylabel('Latency from moving dots (s)');

% putative anticipatory neurons
dblMinAnticipationTime = 0.100;

indPutativeAnticipatory = find( ...
    (vecAllDeltaTRightFromOnset < vecAllOnsetTimeGratingPatches - dblMinAnticipationTime & vecAllDeltaTRightFromOnset < 0.1) | ...
    (vecAllDeltaTLeftFromOnset < vecAllOnsetTimeGratingPatches - dblMinAnticipationTime & vecAllDeltaTLeftFromOnset < 0.1) );

%%
s = 2;
cellSessions{s}
ind = vecAllSessionIndex(indPutativeAnticipatory) == s;
vecAllIndex(indPutativeAnticipatory(ind))

