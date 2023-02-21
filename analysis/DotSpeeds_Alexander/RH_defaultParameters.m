function sParams = rh_defaultParameters()
%rh_defaultParameters
%
%  sStimParamSettings = rh_defaultParameters()
%
% 2023, Alexander Heimel, Robin Haak

strPath = fullfile(fileparts(mfilename('fullpath')),'..','..');
sParams.strOutputPath = strPath; % should be loaded from parameter file instead

sParams.separationFromPrevStimOff = 0.1; % s, time to stay clear of off-response for calculation of spontaneous rate
sParams.boolSmooth = true;
sParams.boolOnlyUseMiddleRangeSpeeds = true;
sParams.boolUseResponseOnset = true;
sParams.boolFitGaussian = false; % compute onset times based on Gaussian fit
sParams.dblThresholdResponsiveZetaP = 1;%0.05;

