function record = analyse_dotspeeds(record,verbose)
%ANALYSE_DOTSPEEDS analysis moving dots stimulus with different speeds
%
%  RECORD = ANALYSE_DOTSPEEDS( RECORD, VERBOSE )
%
%  Analysis responses to stimuli with different speeds
%
%  DeltaT>0 response to stimulus in the past.
%  DeltaT<0 response to stimulus in the future.
%  Calculate DeltaT based on PeakTime
%  Define
%     x0 as the starting position of the center of the stimulus
%     v the speed of the stimulus.
%     Then the center stimulus position at time t is x = x0 + v t
%     If xRF is the receptive field (RF) center, then the time tPass at
%     which the center passes the RF center is,
%         tPass = (xRF - x0)/v
%
%  Compute
%     tLeft = peakTime left stimulus
%     tRight = peakTime right stimulus
%  Assume tLeft = tPassLeft + DeltaT, and tRight = tPassRight + DeltaT
%  Then
%         tLeft = (xRF - x0Left) / v + DeltaT
%         tRight = - (xRF - X0Right) / v + DeltaT
%
%  from this follows:
%         DeltaT = (tLeft + tRight)/2 + (x0Left - x0Right)/(2v)
%         xRF = 1/2 v (tLeft - tRight)
%
% 2022-2023, Robin Haak, Alexander Heimel

if nargin<2 || isempty(verbose)
    verbose = true;
end

sParams.strOutputPath = getdesktopfolder(); % should be loaded from parameter file instead

% strPath = fullfile(fileparts(mfilename('fullpath')),'..','..');
% sParams.strOutputPath = strPath; % should be loaded from parameter file instead

sParams.separationFromPrevStimOff = 0.1; % s, time to stay clear of off-response for calculation of spontaneous rate


%% Load data
strSessionPath = fullfile(sParams.strOutputPath,record.project,'Data_collection',record.dataset,record.subject,record.sessionid);
if ~exist(strSessionPath,'dir')
    logmsg([strSessionPath ' does not exist']);
    return
end

strLog = [record.sessionid '.mat'];
load(fullfile(strSessionPath,strLog),'structEP','sParamsSGL');
load(fullfile(strSessionPath,'spikes.mat'),'vecSpikeCh','vecSpikeSecs');

%% Get stimulus onset times
%to be added: re-alignment of times based on diode data
vecStimOnSecs = structEP.ActOnNI; %NI stream times
vecStimOffSecs = structEP.ActOffNI;

%% Compute measures
vecStimID = structEP.vecStimID;
%vecIDs = unique(vecStimID);
% vecSpeeds = unique(structEP.sAllDots.vecSpeed_pix);
% intSpeed = length(vecSpeeds);
% vecDirs = unique(structEP.sAllDots.vecDirection);
% intDir = length(vecDirs);

%input channels to plot
vecChs  = 200:235; % 1 = tip of the probe  
% notable channels 147,150,151,153,157,167
% vecChs = [147 150 151 153 157 167];

% structEP.sAllDots contains all unique stimuli
record.sStimuli = structEP.sAllDots;
record.intScreenWidth_pix = structEP.sStimParams.intScreenWidth_pix;

vecSpeed_pix = record.sStimuli.vecSpeed_pix;

%vecStimRadius_pix = cellfun( @(x) abs(x(1,1)-x(3,1)), record.sStimuli.vecBoundingRect);
vecStimStartX_pix = cellfun( @(x) mean([x(1,1),x(3,1)]), record.sStimuli.vecBoundingRect) - record.intScreenWidth_pix/2 ; % center of screen is 0

disp('DELTAT CALCULATION IS DEPENDENT ON STIMULUS SPECIFICS! Assuming first 6 stimuli leftwards, and next 6 stimuli rightwards');
    indLeft = 1:6;
    indRight = 7:12;

measures = struct();
for c = 1:length(vecChs) % should be over units
    disp(['Computing channel ' num2str(vecChs(c)) ]) 
    measures(c).intCh = vecChs(c);

    vecSpikeTimesOnChannel = vecSpikeSecs(vecSpikeCh==vecChs(c));
    measures(c).dblRateSpontaneous = computeRateSpontaneous( vecSpikeTimesOnChannel, vecStimOnSecs,vecStimOffSecs, sParams);
    
    vecDuration = zeros(1,length(structEP.sAllDots.stimID));
    vecNRepeats = zeros(1,length(structEP.sAllDots.stimID));
    for i = 1:length(structEP.sAllDots.stimID)
        stimID = structEP.sAllDots.stimID(i);
        indStims = find(vecStimID == stimID);
        
        vecEventStarts = vecStimOnSecs(indStims);
        vecDuration(i) = mean( vecStimOffSecs(indStims) - vecStimOnSecs(indStims));
        vecNRepeats(i) = length(indStims);

        [dblZetaP,sZETA,sRate,~] = zetatest(vecSpikeTimesOnChannel,...
            vecEventStarts,vecDuration(i)); % ~ is essential for correct output
        
        measures(c).dblZetaP(i) = dblZetaP;
        measures(c).cellT{i} = sZETA.vecSpikeT;  % spiketimes with extra 0 and end point
        measures(c).cellSpikeT{i} = measures(c).cellT{i}(2:end-1);  % spiketimes 
        measures(c).vecNSpikes(i) = length( measures(c).cellSpikeT{i} );

        if ~isempty(sRate)
            measures(c).cellRate{i} = sRate.vecRate;
            measures(c).vecPeakTime(i) = sRate.dblPeakTime;
            measures(c).vecPeakRate(i) = sRate.dblPeakRate;
        else
            measures(c).cellRate{i} = [];
            measures(c).vecPeakTime(i) = NaN;
            measures(c).vecPeakRate(i) = NaN;
        end           
        measures(c).vecMeanRate(i) = (length(measures(c).cellSpikeT))/vecDuration(i)/vecNRepeats(i);
    end % stim i
    
    % DeltaT>0 response to stimulus in the past. 
    % DeltaT<0 response to stimulus in the future.
    % Calculate DeltaT based on PeakTime
    % Define 
    %     x0 as the starting position of the center of the stimulus
    %     v the speed of the stimulus.
    %     Then the center stimulus position at time t is x = x0 + v t
    %     If xRF is the receptive field (RF) center, then the time tPass at
    %     which the center passes the RF center is,
    %         tPass = (xRF - x0)/v
    %
    % Compute
    %     tLeft = peakTime left stimulus
    %     tRight = peakTime right stimulus
    % Assume tLeft = tPassLeft + DeltaT, and tRight = tPassRight + DeltaT 
    % Then
    %         tLeft = (xRF - x0Left) / v + DeltaT
    %         tRight = - (xRF - X0Right) / v + DeltaT
    %
    % from this follows:
    %         DeltaT = (tLeft + tRight)/2 + (x0Left - x0Right)/(2v)
    %         xRF = 1/2 v (tLeft - tRight)

    vecTLeft = measures(c).vecPeakTime(indLeft);
    vecTRight = measures(c).vecPeakTime(indRight);

    measures(c).vecPeakXRF_pix = vecSpeed_pix(indLeft) .* (vecTLeft - vecTRight) / 2; 
    measures(c).vecPeakDeltaT = (vecTLeft + vecTRight)/2  + (vecStimStartX_pix(indLeft)-vecStimStartX_pix(indRight)) ./ vecSpeed_pix(indLeft) / 2;

    % Calculate DeltaT based on all spikes
    % Assume rate(t) = rateSpontaneous + s(dir) g(t - DeltaT - tPass)
    %     with \int_-\inf^+\inf g(t) dt = 1, 
    %     support only between (-DeltaT-tPass) and (T-DeltaT-tPass)
    %     and  \int_-\inf^+\inf t g(t) dt = 0 (satisfied if g is symmetric around 0) 
    %     and tPass = (xRF - x0)/v, where xRF is RF position, x0 is start
    %     of stimulus and v is speed.
    % then from Expectation( nSpikes ) = \int_0^T rate(t) follows:
    %     s(dir) = nSpikes - T rateSpontaneous,
    %     with T time interval over which nSpikes is computed
    % and an integration of t * rate(t) over 0 to T, by substituting t'=t+DeltaT-tPass shows
    %     \int_0^T dt t rate(t) = 1/2 T^2 rateSpontaneous + s(dir)(DeltaT + tpass)
    % but the Expectation( spike time of spikes between 0 to T ) * Expectation(nSpikes) = \int_0^T dt t rate(t) 
    % thus 
    %    E(spike time of spikes between 0 to T) * E(nSpikes) = 1/2 T^2 rateSpontaneous + s(dir)(DeltaT + tpass)
    %
    %    Note that E(spike time of spikes between 0 to T) * E(nSpikes) = E(\sum_i t(i))
    % Now define:
    %   SLeft = nSpikesLeft - T rateSpontaneous
    %   SRight = nSpikesRight - T rateSpontaneous
    %   MLeft = E(\sum_i t(i))_left -1/2 T^2 rateSpontaneous 
    %   MRight = E(\sum_i t(i))_right -1/2 T^2 rateSpontaneous 
    % Then:
    %   MLeft / SLeft   = DeltaT + (xRF - x0Left)/v 
    %   MRight / SRight = DeltaT - (xRF - x0Right)/v 
    %
    % =>
    %   DeltaT = 1/2 (MLeft/SLeft + MRight/SRight + x0Left/v - x0Right/v)
    %   xRF = 1/2 v (MLeft/SLeft + MRight/SRight) + 1/2 (x0Left + x0Right)
    %
    %   when only spikes are occuring at tLeft and tRight, this simplifies
    %   to earlier expression.

    dblH = record.intScreenWidth_pix/2;
    vecSLeft = measures(c).vecNSpikes(indLeft) - measures(c).dblRateSpontaneous  * vecDuration(indLeft) .* vecNRepeats(indLeft);
    vecSRight = measures(c).vecNSpikes(indRight) - measures(c).dblRateSpontaneous * vecDuration(indRight) .* vecNRepeats(indRight);
    vecMLeft =   cellfun(@sum,measures(c).cellSpikeT(indLeft))  - 1/2 *  measures(c).dblRateSpontaneous * vecDuration(indLeft).^2 .* vecNRepeats(indLeft);
    vecMRight = cellfun(@sum,measures(c).cellSpikeT(indRight)) - 1/2 *  measures(c).dblRateSpontaneous * vecDuration(indRight).^2 .* vecNRepeats(indRight);

    measures(c).vecMeanDeltaT = 1/2 * (vecMLeft ./ vecSLeft + vecMRight ./ vecSRight + ...
        vecStimStartX_pix(indLeft)./ vecSpeed_pix(indLeft) - vecStimStartX_pix(indRight)./ vecSpeed_pix(indRight) )   ;

    measures(c).vecMeanXRF_pix = vecSpeed_pix(indLeft) .* ...
        (vecMLeft ./ vecSLeft - vecMRight ./ vecSRight) / 2  + ...
        1/2*(vecStimStartX_pix(indLeft) + vecStimStartX_pix(indRight));

    indNaN = abs(measures(c).vecMeanXRF_pix) > dblH;
    measures(c).vecMeanXRF_pix(indNaN) =  NaN;
    measures(c).vecMeanDeltaT(indNaN) = NaN;
end %c

record.measures = measures;

if verbose
    results_dotspeeds(record);
end



function dblRateSpontaneous = computeRateSpontaneous( vecSpikeTimes, vecStimOnSecs,vecStimOffSecs, sParams)
% 
intCount = 0;
dblPeriod = 0;
for i=1:length(vecStimOffSecs)-1
    intCount = intCount + ...
        length(find(vecSpikeTimes>(vecStimOffSecs(i)+sParams.separationFromPrevStimOff) & ...
       vecSpikeTimes<vecStimOnSecs(i+1)));
    dblPeriod = dblPeriod + vecStimOnSecs(i+1) - (vecStimOffSecs(i)+sParams.separationFromPrevStimOff);
end
dblRateSpontaneous = intCount / dblPeriod;
if dblPeriod<5 
    logmsg('Less than 5s to compute spontaneous rate.')
end
if intCount<10 
    logmsg('Less than 10 spikes to compute spontaneous rate.')
end


