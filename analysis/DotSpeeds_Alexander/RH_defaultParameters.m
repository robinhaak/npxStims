function sParams = rh_defaultParameters()
%rh_defaultParameters
%
%  sStimParamSettings = rh_defaultParameters()
%
% 2023, Alexander Heimel, Robin Haak

sParams.strOutputPath = 'Set true path here or set strOutputPath in processparams_local';

sParams.separationFromPrevStimOff = 0.1; % s, time to stay clear of off-response for calculation of spontaneous rate
sParams.boolSmooth = false;
sParams.boolOnlyUseMiddleRangeSpeeds = true;
sParams.boolUseResponseOnset = true;
sParams.boolFitGaussian = true; % compute peak and onset times based on Gaussian fit
sParams.dblThresholdResponsiveZetaP = 0.05;

if exist('processparams_local.m','file')
    sOldParams = sParams;
    sParams = processparams_local( sParams );
    changed_process_parameters(sParams,sOldParams);
end
