%MAKE_ANTICIPATION_FIGURES
%
% Makes some figures for the predictive motion project
%
% 2022, Alexander Heimel

%% Analyse records
record.lab = 'Heimellab';
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
record.measures = [];

dbj = analyse_dotspeeds(record);
for i=2:length(dbj)
    dbj(i) = analyse_dotspeeds(dbj(i)); 
end
sParams.strOutputPath = getdesktopfolder(); % should be loaded from parameter file instead
dbfilename = fullfile(sParams.strOutputPath,'Anticipation','Data_analysis','dbj.mat');
save(dbfilename,'dbj');

%% Make figure speed tuning
record = dbj(1);
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
record = dbj(1);
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




