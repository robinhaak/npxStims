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
record = db(find_record(db,'subject=85235,sessionid=20230117_85235_12_15_45_RH_MovingDots'));
%record = db(find_record(db,'subject=85235,sessionid=20230117_singleUnits_new'));
%record = db(find_record(db,'subject=85235,sessionid=20230117_singleUnits'));

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
record = db(find_record(db,'subject=85235,stimulus=GratingPatches,date=20230117'));

indResponsive = [record.measures.dblResponseMax]>2; 
vecDepth = [record.measures.dblDepth_um]; 
vecPeakTime = [record.measures.dblPeakTime]; 

figure;
hold on

plot(vecPeakTime(indResponsive),vecDepth(indResponsive),'o')
xlabel('Peak time (s)')
ylabel('Depth (\mum)');
set(gca,'YDir','reverse');
axis square
title(recordfilter(record))


%% Compare times from multiunit and sorted
record = db(find_record(db,'subject=85235,sessionid=20230117_singleUnits_new'));

strSessionPath = fullfile(sParams.strOutputPath,record.project,'Data_collection',record.dataset,record.subject,record.sessionid);
strLog = fullfile(strSessionPath,[record.sessionid '.mat']);
sSorted = load(strLog);

record = db(find_record(db,'subject=85235,sessionid=20230117_85235_12_15_45_RH_MovingDots'));

strSessionPath = fullfile(sParams.strOutputPath,record.project,'Data_collection',record.dataset,record.subject,record.sessionid);
strLog = fullfile(strSessionPath,[record.sessionid '.mat']);
sUnsorted = load(strLog);

disp('WORKING HERE TO COMPARE SPIKE TIMES BETWEEN SORTED AND UNSORTED')

% vecSorte
% for i=1:length(sSorted.sCluster)




