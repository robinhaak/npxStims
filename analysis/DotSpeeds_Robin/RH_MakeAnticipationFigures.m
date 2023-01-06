%RH_MAKEANTICIPATIONFIGURES
%create some figures for the motion anticipation project
%
%2022, Alexander Heimel, Robin Haak

%% analyse records 
record.lab = 'Heimellab';
record.project = 'Motion_anticipation';
record.dataset = '22.35.03';
record.subject = '83552';
record.condition = '';
record.stimulus = 'MovingDots';
record.setup = 'Neuropixels';
record.investigator = 'RobinHaak';
record.date = '20221201';
record.sessnr = 3; %dbl
record.sessionid = [record.subject '_' record.date '_' sprintf('%03d',record.sessnr)];
record.logfile = '';
record.version = '1.0';
record.path = '\\vs03\VS03-CSF-1\Haak\';
record.measures = [];

%loop through records
dbj = RH_AnalyseDotSpeeds(record);
for i=2:length(dbj)
    dbj(i) = RH_AnalyseDotSpeeds(dbj(i)); 
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




