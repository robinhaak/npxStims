function [dblPeakTime,dblOnsetTime] = fitGaussianPSTH( vecSpikeT, dblSpont, dblPeakTime, boolFitSpont, boolFitUpToPeak, verbose )
%fitGaussianPSTH fits erf to cumulative response
%
% [dblPeakTime,dblOnsetTime] = fitGaussianPSTH( vecSpikeT, dblSpont, dblPeakTime, boolFitSpont, boolFitUpToPeak, verbose )
%
%    vecSpikeT are spike times relative to stimulus onset
%    dblSpont is spontaneous rate (sp/s), used as begin point in fitting
%    dblPeakTime is peak time, used as begin point in fitting
%
% 2023, Alexander Heimel

dblOnsetTime = NaN;
dblMaxTime = vecSpikeT(end);

if nargin<3 || isempty(dblPeakTime) || isnan(dblPeakTime)
    dblPeakTime = dblMaxTime/2;
end
if nargin<4 || isempty(boolFitSpont)
    boolFitSpont = false;
end
if nargin<5 || isempty(boolFitUpToPeak)
    boolFitUpToPeak = false;
end
if nargin<6 || isempty(verbose)
    verbose = false;
end


if ~boolFitUpToPeak
    x = vecSpikeT;
    intNumSpikes = length(x);
    
    if intNumSpikes < 4
        return;
    end
    
    y = (1:intNumSpikes)';

    P0(1) = 0; % Offset number of spikes
    lb(1) = -10;
    ub(1) = 10;
    P0(2) = 100; % Number of spikes for first Gaussian
    lb(2) = 0;
    ub(2) = intNumSpikes;
    P0(3) = dblPeakTime; % s, Time of first peak
    lb(3) = 0;
    ub(3) = dblMaxTime;
    P0(4) = 1; % s, Width of first peak
    lb(4) = 0;
    ub(4) = dblMaxTime;
    P0(5) = dblSpont; % sp/s, spontaneous rate
    lb(5) = 0;
    ub(5) = 3 * dblSpont + 5;
else % boolFitUpToPeak
    ind = find(vecSpikeT>dblPeakTime + 0.05,1); % 50 ms slack
    x = vecSpikeT(1:ind);
    intNumSpikes = length(x);
    if intNumSpikes < 4
        return;
    end
    
    y = (1:intNumSpikes)';

    
    
    dblMaxTime = x(end);

    
    P0(1) = 0; % Offset number of spikes
    lb(1) = -10;
    ub(1) = 10;
    P0(2) = 100; % Number of spikes for first Gaussian
    lb(2) = 0;
    ub(2) = intNumSpikes;
    P0(3) = dblPeakTime; % s, Time of first peak
    lb(3) = 0;
    ub(3) = dblMaxTime;
    P0(4) = 1; % s, Width of first peak
    lb(4) = 0;
    ub(4) = dblMaxTime;
    P0(5) = dblSpont; % sp/s, spontaneous rate
    lb(5) = 0;
    ub(5) = 3 * dblSpont + 5;
end

opt = optimoptions('lsqcurvefit');
opt.Display = 'none';

if boolFitSpont
    model = @(P,x) P(1) + P(5)*x + P(2)*0.5*(1+erf( (x-P(3))/P(4))) ;
else
    model = @(P,x) P(1) + dblSpont*x + P(2)*0.5*(1+erf( (x-P(3))/P(4))) ;
    P0 = P0(1:4);
    lb = lb(1:4);
    ub = ub(1:4);
end


try
    Pfit = lsqcurvefit(model,P0,x,y,lb,ub,opt);
catch me
    logmsg(me.message)
    return
end

dblPeakTime = Pfit(3);
dblOnsetTime = Pfit(3) - 2*Pfit(4);


if verbose
    x = vecSpikeT;
    y = (1:length(x))';

    figure('Name','Timing','NumberTitle','off');
    subplot(2,1,1)
    histogram(x,'BinWidth',0.010);
    hold on;
    plot(dblPeakTime*[1 1],ylim,'-b');
    plot(dblOnsetTime*[1 1],ylim,'-g');
    
    subplot(2,1,2)
    plot(x,y,'.')
    hold on;
    modelpred = model(Pfit,x);
    plot(x,modelpred,'r-');
    plot(dblPeakTime*[1 1],ylim,'-b');
    plot(dblOnsetTime*[1 1],ylim,'-g');
end
