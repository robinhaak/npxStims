function sParams = RH_defaultParameters( )
%RH_defaultParameters
%
%  sStimParamSettings = rh_defaultParameters()
%
% 2023, Alexander Heimel, Robin Haak

sParams.strOutputPath = 'C:\Users\haak\Desktop';

sParams.separationFromPrevStimOff = 0.1; % s, time to stay clear of off-response for calculation of spontaneous rate
        sParams.dblBinWidth = 0.01; % Binwidth for PSTH for peak rate calculation

sParams.boolSmooth = false;
sParams.boolOnlyUseMiddleRangeSpeeds = true;
sParams.boolUseResponseOnset = true;
sParams.boolFitGaussian = false; % compute peak and onset times based on Gaussian fit
sParams.dblThresholdResponsiveZetaP = 0.05;
sParams.dblOnsetResponseThreshold = 0.5;

% Color scheme
sParams.clrLeft = [1 0 0];
sParams.clrRight = [0 0 1];
sParams.clrPatches = [0 0.6 0];



if exist('processparams_local.m','file')
    sOldParams = sParams;
    sParams = processparams_local( sParams );
    changed_process_parameters(sParams,sOldParams);
end
