function [dblFirstPeakTime,dblOnsetTime] = fitDoubleGaussianPSTH( vecSpikeT, dblSpont, dblPeakTime, verbose )
%fitDoubleGaussianPSTH fits erf to cumulative response
%
% [dblPeakTime,dblOnsetTime] = fitDoubleGaussianPSTH( vecSpikeT, dblSpont, dblPeakTime, verbose)
%    Input:
%      vecSpikeT are spike times relative to stimulus onset
%      dblSpont is spontaneous rate (sp/s), used as begin point in fitting
%      dblPeakTime is peak time, used as begin point in fitting
%
% 2023, Alexander Heimel

dblMaxTime = vecSpikeT(end);

if nargin<3 || isempty(dblPeakTime) || isnan(dblPeakTime)
    dblPeakTime = dblMaxTime/2;
end
if nargin<4 || isempty(verbose)
    verbose = false;
end

intNumSpikes = length(vecSpikeT);
x = vecSpikeT;
y = (1:intNumSpikes)';

P0 = [];
lb = []; % Lower bounds
ub = []; % Upper bounds

P0(1) = 0; % Offset number of spikes
lb(1) = -10; 
ub(1) = 10; 
P0(2) = 100; % Number of spikes for first Gaussian
lb(2) = 0;
ub(2) = intNumSpikes;
P0(3) = 0.10; % s, Time of first peak
lb(3) = 0;
ub(3) = dblMaxTime;
P0(4) = 1; % s, Width of first peak
lb(4) = 0;
ub(4) = dblMaxTime;
P0(5) = dblSpont; % sp/s, spontaneous rate
lb(5) = 0;
ub(5) = 3 * dblSpont + 5;
P0(6) = 100; % Number of spikes for second Gaussian
lb(6) = 0;
ub(6) = intNumSpikes;
P0(7) = P0(3) + 0.2; % s, Time of second peak
lb(7) = 0;
ub(7) = dblMaxTime;
P0(8) = 1; % s, Width of second peak
lb(8) = 0;
ub(8) = dblMaxTime;

model = @(P,x) P(1) + P(5)*x + P(2)*0.5*(1+erf( (x-P(3))/P(4))) + P(6)*0.5*(1+erf( (x-P(7))/P(8)));
try
    opt = optimoptions('lsqcurvefit');
    opt.Display = 'none';
    Pfit = lsqcurvefit(model,P0,x,y,lb,ub,opt);
catch me
    me.message
    keyboard
end

dblFirstPeakTime = Pfit(3);
dblOnsetTime = Pfit(3) - 2*Pfit(4);

if verbose
    figure;
    plot(x,y,'.')
    hold on;
    modelpred = model(Pfit,x);
    plot(x,modelpred,'r-');
    plot(dblFirstPeakTime*[1 1],ylim,'-b');
    plot(dblOnsetTime*[1 1],ylim,'-g');
end
