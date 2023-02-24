function [dblPeakTime,dblOnsetTime] = fitGaussianPSTH( vecSpikeT, dblSpont, dblPeakTime, verbose )
%fitGaussianPSTH fits erf to cumulative response
%
% [dblPeakTime,dblOnsetTime] = fitGaussianPSTH( vecSpikeT, dblSpont, dblPeakTime, verbose)
%
%    vecSpikeT are spike times relative to stimulus onset
%    dblSpont is spontaneous rate (sp/s), used as begin point in fitting
%    dblPeakTime is peak time, used as begin point in fitting
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

P0 = [ 0     100        dblPeakTime   1          dblSpont];
lb = [-10       0          0                  0                 0];
ub = [ 10  intNumSpikes   dblMaxTime      dblMaxTime         3*dblSpont+5];
model = @(P,x) P(1) + P(5)*x + P(2)*0.5*(1+erf( (x-P(3))/P(4))) ;
try
    opt = optimoptions('lsqcurvefit');
    opt.Display = 'none';
    Pfit = lsqcurvefit(model,P0,x,y,lb,ub,opt);
catch me
    me.message
    keyboard
end

dblPeakTime = Pfit(3);
dblOnsetTime = Pfit(3) - 2*Pfit(4);

if verbose
    figure;
    plot(x,y,'.')
    hold on;
    modelpred = model(Pfit,x);
    plot(x,modelpred,'r-');
    plot(dblPeakTime*[1 1],ylim,'-b');
    plot(dblOnsetTime*[1 1],ylim,'-g');
end
