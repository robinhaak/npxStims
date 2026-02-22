% % =========================================================
% % Plot x z positions of all electrodes and saved channels
% %
% function plotSaved(xCoord, yCoord, shankInd, meta)
%     % get geometry
%     g = getGeomParams(meta);
% 
%     % calculate positions on one shank
%     nCol = g.elecPerShank/g.rowsPerShank;
%     rowInd = 0:g.elecPerShank-1;
%     rowInd = floor(rowInd./nCol);
%     oddRows = logical(mod(rowInd,2));
%     evenRows = ~oddRows;
% 
%     colInd = 0:g.elecPerShank-1;
%     colInd = mod(colInd,nCol);
% 
%     xall = colInd*g.horzPitch;
%     xall(evenRows) = xall(evenRows) + g.even_xOff;
%     xall(oddRows) = xall(oddRows) + g.odd_xOff;
% 
%     yall = rowInd*g.vertPitch;
% 
%     figure('Name','shank view','Units','Normalized', 'Position', [0.2,0.1,0.3,0.8])
%     for sI = 0:g.nShank-1
%         cc = find(shankInd == sI);
%         scatter( g.shankPitch*sI + xall, yall, 3, 'k', 'square' ); hold on;
%         scatter( g.shankPitch*sI + xCoord(cc), yCoord(cc), 30, 'green', 'square', 'filled' ); hold on; 
%     end
%     xlim([min(xall)-5, (g.nShank-1)*g.shankPitch + max(xall)+5]);
%     ylim([min(yall)-5, max(yall)+5]);
%     hold off;
% end % plotSaved

function plotSaved(xCoord, yCoord, shankInd, meta)

    % Get geometry parameters
    g = getGeomParams(meta);

    % ----- Compute electrode layout for ONE shank -----

    % Number of columns per shank
    nCol = g.elecPerShank / g.rowsPerShank;

    % Row indices
    rowInd = floor((0:g.elecPerShank-1) ./ nCol);

    % Odd/even row flags (for staggered layout)
    oddRows  = mod(rowInd, 2) == 1;
    evenRows = ~oddRows;

    % Column indices
    colInd = mod(0:g.elecPerShank-1, nCol);

    % X positions (with stagger offsets)
    xall = colInd * g.horzPitch;
    xall(evenRows) = xall(evenRows) + g.even_xOff;
    xall(oddRows)  = xall(oddRows)  + g.odd_xOff;

    % Y positions
    yall = rowInd * g.vertPitch;

    % ----- Plot -----

    figure('Name','shank view', ...
           'Units','Normalized', ...
           'Position',[0.2,0.1,0.3,0.8]);

    hold on

    % Loop over shanks (1-based indexing)
    for sI = 1:g.nShank

        % Indices of selected channels on this shank
        cc = (shankInd == sI);

        % Horizontal offset for this shank
        xOffset = g.shankPitch * (sI - 1);

        % Plot all electrodes (black)
        scatter(xOffset + xall, yall, ...
                3, 'k', 'square');

        % Overlay selected electrodes (green filled)
        scatter(xOffset + xCoord(cc), ...
                yCoord(cc), ...
                30, 'green', 'square', 'filled');
    end

    % Axis limits
    xlim([min(xall)-5, ...
          (g.nShank-1)*g.shankPitch + max(xall)+5]);

    ylim([min(yall)-5, max(yall)+5]);

    hold off
end