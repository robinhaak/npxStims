
function [centerX, centerY, sigmaX, sigmaY, amplitude] = fitGaussian2D_simple(rates, Xgrid, Ygrid)
    % Simple 2D Gaussian fitting (circular, no rotation)
    
    % Initial guess from data
    [maxVal, maxIdx] = max(rates(:));
    x0 = Xgrid(maxIdx);
    y0 = Ygrid(maxIdx);
    sigma0 = mean([range(Xgrid(:)), range(Ygrid(:))]) / 4;
    baseline = min(rates(:));
    
    % Define circular Gaussian model
    gauss2D = @(p, X, Y) p(4) + p(5) * exp(-((X - p(1)).^2 + (Y - p(2)).^2) / (2 * p(3)^2));
    
    % Initial parameters: [x0, y0, sigma, baseline, amplitude]
    params0 = [x0, y0, sigma0, baseline, maxVal - baseline];
    
    % Lower and upper bounds
    lb = [min(Xgrid(:)), min(Ygrid(:)), 0.5*sigma0, 0, 0];
    ub = [max(Xgrid(:)), max(Ygrid(:)), 2*sigma0, maxVal, 2*maxVal];
    
    % Fit using lsqcurvefit
    options = optimoptions('lsqcurvefit', 'Display', 'off', 'MaxIterations', 100);
    
    try
        % Create coordinate arrays
        coords = [Xgrid(:), Ygrid(:)];
        ratesVec = rates(:);
        
        % Fit
        params = lsqcurvefit(@(p, coords) gauss2D(p, coords(:,1), coords(:,2)), ...
            params0, coords, ratesVec, lb, ub, options);
        
        % Extract parameters
        centerX = params(1);
        centerY = params(2);
        sigmaX = abs(params(3));
        sigmaY = abs(params(3)); % circular
        amplitude = params(5);
        
        % Quality check: center should be within grid bounds
        if centerX < min(Xgrid(:)) || centerX > max(Xgrid(:)) || ...
           centerY < min(Ygrid(:)) || centerY > max(Ygrid(:))
            centerX = NaN;
            centerY = NaN;
        end
        
    catch
        % If fitting fails, return NaN
        centerX = NaN;
        centerY = NaN;
        sigmaX = NaN;
        sigmaY = NaN;
        amplitude = NaN;
    end
end