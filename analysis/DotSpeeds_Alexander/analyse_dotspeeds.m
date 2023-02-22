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

sParams = RH_defaultParameters();

record = RH_analyse_synchronization(record,verbose);

%% Load data
% check which channels or clusters to analyse
vecClustersToAnalyze = get_channels2analyze( record );

strSessionPath = fullfile(sParams.strOutputPath,record.project,'Data_collection',record.dataset,record.subject,record.sessionid);
if ~exist(strSessionPath,'dir')
    logmsg([strSessionPath ' does not exist']);
    return
end

strLog = fullfile(strSessionPath,[record.sessionid '.mat']);
sVars = whos('-file',strLog);
if contains('sCluster',{sVars.name})
    % single unit data
    load(strLog,'sCluster','structEP','sSynthData'); %rm sSynthData later
    

    cellExpType = cellfun(@(x) x.structEP.strExpType, sSynthData.cellStim, 'UniformOutput', false);
    ind = find(contains(cellExpType,record.stimulus));
    if isempty(ind)
        logmsg(['Could not find stimulus ' record.stimulus ' in '  strLog]);
        return
    end
    if length(ind)>1
        logmsg(['More than one stimulus ' record.stimulus ' in '  strLog]);
        logmsg('Using first stimulus');
        ind = 1;
    end
    
    structEP = sSynthData.cellStim{ind}.structEP;
    sCluster = sSynthData.sCluster;

    
    if isempty(vecClustersToAnalyze)
        vecClustersToAnalyze = unique([sCluster.IdxClust]);
    end
    intNumClusters = length(vecClustersToAnalyze);
    vecSpikeTimesOfCluster = cell(intNumClusters,1);
    for c = 1:intNumClusters
        indClus = find([sCluster.IdxClust]==vecClustersToAnalyze(c),1);
        if ~isempty(indClus)
            vecSpikeTimesOfCluster{c} = sCluster(indClus).SpikeTimes;
        end
    end
    
        vecStimOnTime = structEP.vecStimOnTime;
        vecStimOffTime = structEP.vecStimOffTime;
%     sMeta = sSynthData.sMetaNI;
%     dblStartT = str2double(sMeta.firstSample) / str2double(sMeta.niSampRate);
%     vecStimOnTime = structEP.ActOnNI - dblStartT;
%     vecStimOffTime = structEP.ActOffNI - dblStartT;


else
    % unsorted channels
    load(strLog,'structEP');
    load(fullfile(strSessionPath,'spikes.mat'),'vecSpikeCh','vecSpikeSecs');
    
    if isempty(vecClustersToAnalyze)
        vecClustersToAnalyze = unique(vecSpikeCh);
    end
    intNumClusters = length(vecClustersToAnalyze);
    vecSpikeTimesOfCluster = cell(intNumClusters,1);
    for c = 1:intNumClusters
        vecSpikeTimesOfCluster{c} = vecSpikeSecs(vecSpikeCh==vecClustersToAnalyze(c));
    end
    
    %to be added: re-alignment of times based on diode data
    vecStimOnTime = structEP.ActOnNI; %NI stream times
    vecStimOffTime = structEP.ActOffNI;
end

%% Compute measures
vecStimID = structEP.vecStimID;
%vecIDs = unique(vecStimID);
% vecSpeeds = unique(structEP.sAllDots.vecSpeed_pix);
% intSpeed = length(vecSpeeds);
% vecDirs = unique(structEP.sAllDots.vecDirection);
% intDir = length(vecDirs);

% structEP.sAllDots contains all unique stimuli
record.sStimuli = structEP.sAllDots;
record.intScreenWidth_pix = structEP.sStimParams.intScreenWidth_pix;

vecSpeed_pix = record.sStimuli.vecSpeed_pix;

vecStimStartX_pix = cellfun( @(x) mean([x(1,1),x(3,1)]), record.sStimuli.vecBoundingRect) - record.intScreenWidth_pix/2 ; % center of screen is 0
record.sStimuli.vecStimStartX_pix = vecStimStartX_pix;

logmsg('Assuming Left and right stimuli are in the same order for delta T calculation');
indLeft = find([structEP.sAllDots.vecDirection]==0);
indRight = find([structEP.sAllDots.vecDirection]==180);

warning('off','zetatest:InsufficientSamples');

measures = struct([]);
indMeasures = 0;

dblIntertrialInterval =  structEP.sStimParams.dblSecsPreBlank + structEP.sStimParams.vecSecsPostBlank(1);

for c = 1:length(vecClustersToAnalyze) % over clusters or channels
    
    
    logmsg(['Computing cluster/channel index ' num2str(vecClustersToAnalyze(c)) ...
        ' (' num2str(c) ' of '  num2str(intNumClusters) '), ' num2str(length(vecSpikeTimesOfCluster{c})) ' spikes'])
    if length(vecSpikeTimesOfCluster{c}) < 1
        continue
    end
    
    indMeasures = indMeasures + 1;
    measure.intIndex = vecClustersToAnalyze(c);
    
    if exist('sCluster','var')
        ind = find([sCluster.IdxClust]==measure.intIndex);
        measure.dblDepth_um = sCluster(ind).Depth; %#ok<FNDSB> 
    else
        measure.dblDepth_um = 1000 - c; %before: measure.depth
    end
    
    measure.dblRateSpontaneous = computeRateSpontaneous( vecSpikeTimesOfCluster{c}, vecStimOnTime,vecStimOffTime, sParams);
    
    vecDuration = zeros(1,length(structEP.sAllDots.stimID));
    vecNRepeats = zeros(1,length(structEP.sAllDots.stimID));
     for i = 1:length(structEP.sAllDots.stimID)
        stimID = structEP.sAllDots.stimID(i);
        indStims = find(vecStimID == stimID);
        
        vecEventStarts = vecStimOnTime(indStims);
        vecDuration(i) = mean(vecStimOffTime(indStims) - vecStimOnTime(indStims));
        vecNRepeats(i) = length(indStims);
        
        [dblZetaP,sZETA,sRate,~] = zetatest(vecSpikeTimesOfCluster{c},...
            vecEventStarts,vecDuration(i)+dblIntertrialInterval); % ~ is essential for correct output
       
        
        
        measure.vecZetaP(i) = dblZetaP;
        measure.cellSpikeT{i} =sZETA.vecSpikeT(2:end-1);  % spiketimes, ZETA has padded extra points at beginning and end
        measure.vecNSpikes(i) = length(measure.cellSpikeT{i} );
        
        measure.cellRate{i} = [];
        measure.cellTime{i} = []; 
        measure.vecPeakTime(i) = NaN;
        measure.vecPeakRate(i) = NaN;
        
        if sParams.boolSmooth && ~isempty(vecSpikeTimesOfCluster{c})
            
            
                        dblBase = 1.5;
                        intSmoothSd = 2;
                        dblMinScale = round(log(1/10) / log(dblBase));
                        [sRate.vecT,sRate.vecRate,~] = getIFR(vecSpikeTimesOfCluster{c},vecEventStarts,vecDuration(i)+dblIntertrialInterval,intSmoothSd,dblMinScale,dblBase,false);
                        if ~isempty(sRate.vecT)
                            [sRate.dblPeakRate,ind] = max(sRate.vecRate);
                            sRate.dblPeakTime = sRate.vecT(ind);
                        else
                            sRate = [];
                        end
        end
        
        if ~isempty(sRate)
            measure.cellRate{i} = sRate.vecRate;
            measure.cellTime{i} = sRate.vecT;
            measure.vecPeakTime(i) = sRate.dblPeakTime;
            measure.vecPeakRate(i) = sRate.dblPeakRate;
            
            if sParams.boolFitGaussian
                logmsg(['Fitting ' num2str(i) ' of ' num2str(length(structEP.sAllDots.stimID))]);
                
                % fit erf to cumulative response
                intNumSpikes = length(sZETA.vecSpikeT);
                dblMaxTime = sZETA.vecSpikeT(end);
                x = sZETA.vecSpikeT;
                y = (1:intNumSpikes)';
                
                dblSpont = measure.dblRateSpontaneous*vecNRepeats(i);
                dblPeakTime = sRate.dblPeakTime;
                if isnan(dblPeakTime)
                    dblPeakTime = dblMaxTime/2;
                end
                
                P0 = [ 0     100        dblPeakTime   1          dblSpont];
                lb = [-10       0          0                  0                 0];
                ub = [ 10  intNumSpikes   dblMaxTime      dblMaxTime         3*dblSpont+5*vecNRepeats(i)];
                model = @(P,x) P(1) + P(5)*x + P(2)*0.5*(1+erf( (x-P(3))/P(4))) ;
                try
                    opt = optimoptions('lsqcurvefit');
                    opt.Display = 'none';
                    Pfit = lsqcurvefit(model,P0,x,y,lb,ub,opt);
                catch me
                    me.message
                    keyboard
                end
                measure.vecPeakTime(i) = Pfit(3);
                measure.vecOnsetTime(i) = Pfit(3) - 2*Pfit(4);
                %                 logmsg(['Erf fit peak time = '  num2str(measure.vecPeakTime(i))]);
                %                 logmsg(['Erf fit onset time = '  num2str(measure.vecOnsetTime(i))]);
                if 0
                    figure;
                    plot(x,y,'.')
                    hold on;
                    modelpred = model(Pfit,x);
                    plot(x,modelpred,'r-');
                    plot(measure.vecPeakTime(i)*[1 1],ylim,'-b');
                    plot(measure.vecOnsetTime(i)*[1 1],ylim,'-g');
                end 
                
%                 % fit gaussian to response
%                 strGaussEqn = 'a*exp(-((x-b)/sigma)^2/2)';
%                 d = measure.dblRateSpontaneous;
%                 a = measure.vecPeakRate(i)-d;
%                 b = measure.vecPeakTime(i);
%                 sigma = 0.2;  % notice: half sigma
%                 startPoints = [a b sigma];
%                 tic
%                 f = fit(measure.cellTime{i},measure.cellRate{i}-d,strGaussEqn,'Start', startPoints);
%                 toc
%                 ind = find(measure.cellTime{i}>measure.vecPeakTime(i)+f.sigma,1);
%                 if ~isempty(ind) && ind>2
%                     % now refit around the lefthand side of the peak
%                     startPoints = [f.a f.b f.sigma*1.1];
%                     f = fit(measure.cellTime{i}(1:ind),measure.cellRate{i}(1:ind)-d,strGaussEqn,'Start', startPoints);
%                 end
%                 measure.vecOnsetTime(i) = f.b - f.sigma*2;
%                 measure.vecPeakTime(i) = f.b;
%                 logmsg(['Gauss fit peak time = '  num2str(measure.vecPeakTime(i))]);
%                 logmsg(['Gauss fit onset time = '  num2str(measure.vecOnsetTime(i))]);
%                 if 1
%                     figure;
%                     plot(measure.cellTime{i},measure.cellRate{i},'k-');
%                     hold on
%                     plot(measure.cellTime{i}, d+f.a*exp(-((measure.cellTime{i}-f.b)/f.sigma).^2))
%                     plot(measure.vecOnsetTime(i),measure.dblRateSpontaneous,'o');
%                     plot(measure.vecPeakTime(i)*[1 1],ylim,'-b');
%                 end
            else
                sParams.dblOnsetResponseThreshold = 0.5;
                dblThreshold = measure.dblRateSpontaneous + sParams.dblOnsetResponseThreshold*(measure.vecPeakRate(i)-measure.dblRateSpontaneous);
                ind = find(measure.cellRate{i}>dblThreshold,1);
                if isempty(ind) || ind==1 % probably spontaneous rate is wrong
                    dblRateSpontaneous = median(measure.cellRate{i});
                    dblThreshold = dblRateSpontaneous + sParams.dblOnsetResponseThreshold*(measure.vecPeakRate(i)-dblRateSpontaneous);
                    ind = find(measure.cellRate{i}>dblThreshold,1);
                end
                if ~isempty(ind)
                    measure.vecOnsetTime(i) = measure.cellTime{i}(ind);
                else
                    measure.vecOnsetTime(i) = NaN;
                end
            end
            
        end
        measure.vecMeanRate(i) = (length(measure.cellSpikeT{i}))/vecDuration(i)/vecNRepeats(i);
        
        
        
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
    
    vecTLeft = measure.vecPeakTime(indLeft);
    vecTRight = measure.vecPeakTime(indRight);
    
    measure.vecPeakXRF_pix = vecSpeed_pix(indLeft) .* (vecTLeft - vecTRight) / 2;
    measure.vecPeakDeltaT = (vecTLeft + vecTRight)/2  + (vecStimStartX_pix(indLeft)-vecStimStartX_pix(indRight)) ./ vecSpeed_pix(indLeft) / 2;
    
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
    vecSLeft = measure.vecNSpikes(indLeft) -measure.dblRateSpontaneous  * vecDuration(indLeft) .* vecNRepeats(indLeft);
    vecSRight = measure.vecNSpikes(indRight) -measure.dblRateSpontaneous * vecDuration(indRight) .* vecNRepeats(indRight);
    vecMLeft =   cellfun(@sum,measure.cellSpikeT(indLeft))  - 1/2 * measure.dblRateSpontaneous * vecDuration(indLeft).^2 .* vecNRepeats(indLeft);
    vecMRight = cellfun(@sum,measure.cellSpikeT(indRight)) - 1/2 * measure.dblRateSpontaneous * vecDuration(indRight).^2 .* vecNRepeats(indRight);
    
    measure.vecMeanDeltaT = 1/2 * (vecMLeft ./ vecSLeft + vecMRight ./ vecSRight + ...
        vecStimStartX_pix(indLeft)./ vecSpeed_pix(indLeft) - vecStimStartX_pix(indRight)./ vecSpeed_pix(indRight) )   ;
    
    measure.vecMeanXRF_pix = vecSpeed_pix(indLeft) .* ...
        (vecMLeft ./ vecSLeft - vecMRight ./ vecSRight) / 2  + ...
        1/2*(vecStimStartX_pix(indLeft) + vecStimStartX_pix(indRight));
    
    indNaN = abs(measure.vecMeanXRF_pix) > dblH;
    measure.vecMeanXRF_pix(indNaN) =  NaN;
    measure.vecMeanDeltaT(indNaN) = NaN;
    
    % Comput fits on based of peak times for different velocities
    % peakTime =  1/v (xRF - x0) + DeltaT
    % we can fit xRF and DeltaT independently for each direction
    
    vecLeft = false(1,length(vecSpeed_pix));
    vecLeft(indLeft) = true;
    vecRight = false(1,length(vecSpeed_pix));
    vecRight(indRight) = true;
    vecResponsive = (measure.vecZetaP<sParams.dblThresholdResponsiveZetaP );
%     vecResponsive = vecResponsive &  (measure.vecPeakRate>measure.dblRateSpontaneous+2); % remove suppressed responses
    
    
    if sParams.boolOnlyUseMiddleRangeSpeeds
        vecMiddleRange = record.sStimuli.vecSpeed_deg>5 & record.sStimuli.vecSpeed_deg<100;
        vecResponsive = vecResponsive & vecMiddleRange;
    end
    
    if sum(vecLeft & vecResponsive)>1
        vecTime = measure.vecPeakTime;
        measure.lmLeft = compact(fitlm(1./vecSpeed_pix(vecLeft & vecResponsive),vecTime(vecLeft & vecResponsive)));
        measure.dblDeltaTLeft = measure.lmLeft.Coefficients.Estimate(1);
        measure.dblXRFLeft_pix = measure.lmLeft.Coefficients.Estimate(2) + mean(vecStimStartX_pix(indLeft));
        
        vecTime = measure.vecOnsetTime;
        measure.lmLeftFromOnset = compact(fitlm(1./vecSpeed_pix(vecLeft & vecResponsive),vecTime(vecLeft & vecResponsive)));
        measure.dblDeltaTLeftFromOnset = measure.lmLeftFromOnset.Coefficients.Estimate(1);
        measure.dblXRFLeftFromOnset_pix = measure.lmLeftFromOnset.Coefficients.Estimate(2) + mean(vecStimStartX_pix(indLeft));
    else
        measure.lmLeft = [];
        measure.dblDeltaTLeft = NaN;
        measure.dblXRFLeft_pix = NaN;
        measure.lmLeftFromOnset = [];
        measure.dblDeltaTLeftFromOnset = NaN;
        measure.dblXRFLeftFromOnset_pix = NaN;
    end
    if sum(vecRight & vecResponsive)>1
        vecTime = measure.vecPeakTime;
        measure.lmRight = compact(fitlm(1./vecSpeed_pix(vecRight & vecResponsive),vecTime(vecRight & vecResponsive)));
        measure.dblDeltaTRight = measure.lmRight.Coefficients.Estimate(1);
        measure.dblXRFRight_pix = -measure.lmRight.Coefficients.Estimate(2) + mean(vecStimStartX_pix(indRight));
        
        vecTime = measure.vecOnsetTime;
        measure.lmRightFromOnset = compact(fitlm(1./vecSpeed_pix(vecRight & vecResponsive),vecTime(vecRight & vecResponsive)));
        measure.dblDeltaTRightFromOnset = measure.lmRightFromOnset.Coefficients.Estimate(1);
        measure.dblXRFRightFromOnset_pix = -measure.lmRightFromOnset.Coefficients.Estimate(2) + mean(vecStimStartX_pix(indRight));
        
    else
        measure.lmRight = [];
        measure.dblDeltaTRight = NaN;
        measure.dblXRFRight_pix = NaN;
        measure.lmRightFromOnset = [];
        measure.dblDeltaTRightFromOnset = NaN;
        measure.dblXRFRightFromOnset_pix = NaN;
        
    end
    
    if sum(measure.vecZetaP<0.05)>2 
        measure.boolResponsive = true;
    else
        measure.boolResponsive = false;
    end
    if ~measure.boolResponsive
        % reduce database size
        logmsg('Cluster not responsive. Removing PSTHs');
        measure.cellTime = {};
        measure.cellRate = {};
    end
    
    measure = rmfield(measure,'cellSpikeT'); % to reduce database size

    if indMeasures == 1
        measures = measure;
    else
        measures(indMeasures) = measure;
    end
end % cluster c


if isfield(record,'measures') && ~isempty(record.measures) && length(measures)>=1
    % Merge with old measures
    sMeasuresOld = record.measures;
    vecIndicesOld = [sMeasuresOld.intIndex];
    vecIndicesNew = [measures.intIndex];
    vecOverlap = intersect(vecIndicesOld,vecIndicesNew);
    sMeasuresOld(ismember(vecIndicesOld,vecOverlap)) = [];
    try
        measures = [sMeasuresOld measures];
        vecIndices = [measures.intIndex];
        [~,ind] = sort(vecIndices);
        measures = measures(ind);
    catch me
        logmsg('Cannot merge with old measures. Removing old measures.')
    end
end


record.measures = measures;

logmsg(['Analyzed ' recordfilter(record)]);
end



