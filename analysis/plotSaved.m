% =========================================================
% Plot x z positions of all electrodes and saved channels
%
function plotSaved(xCoord, yCoord, shankInd, meta)
    % get geometry
    g = getGeomParams(meta);

    % calculate positions on one shank
    nCol = g.elecPerShank/g.rowsPerShank;
    rowInd = 0:g.elecPerShank-1;
    rowInd = floor(rowInd./nCol);
    oddRows = logical(mod(rowInd,2));
    evenRows = ~oddRows;

    colInd = 0:g.elecPerShank-1;
    colInd = mod(colInd,nCol);

    xall = colInd*g.horzPitch;
    xall(evenRows) = xall(evenRows) + g.even_xOff;
    xall(oddRows) = xall(oddRows) + g.odd_xOff;

    yall = rowInd*g.vertPitch;
    
    figure('Name','shank view','Units','Normalized', 'Position', [0.2,0.1,0.3,0.8])
    for sI = 0:g.nShank-1
        cc = find(shankInd == sI);
        scatter( g.shankPitch*sI + xall, yall, 3, 'k', 'square' ); hold on;
        scatter( g.shankPitch*sI + xCoord(cc), yCoord(cc), 30, 'green', 'square', 'filled' ); hold on; 
    end
    xlim([min(xall)-5, (g.nShank-1)*g.shankPitch + max(xall)+5]);
    ylim([min(yall)-5, max(yall)+5]);
    hold off;
end % plotSaved
