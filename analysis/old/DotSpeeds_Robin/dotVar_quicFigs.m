
indLeft = [1 3 5];
indRight = [2 4  6];

vecDepth = [];
vecDeltaPeak_l = [];
vecDeltaPeak_r = [];

for i = 1:numel(record.measures)
    if any(record.measures(i).dblZetaP([indLeft indRight])<0.1)
        figure;
        subplot(2,1,1); hold on
        for j = indLeft
            plot(record.measures(i).cellT{j},record.measures(i).cellRate{j});
        end
        scatter(record.measures(i).vecPeakTime(indLeft), record.measures(i).vecPeakRate(indLeft),'xk', 'HandleVisibility','off');
        legend('Continuous', 'Appear', 'Disappear');
        ylabel('Rate (spk/s)');
        fixfig;

        %
        vecDeltaPeak_l = [vecDeltaPeak_l record.measures(i).vecPeakTime(3)-record.measures(i).vecPeakTime(1)];

        subplot(2,1,2); hold on
        for j = indRight
            plot(record.measures(i).cellT{j},record.measures(i).cellRate{j});
        end
        scatter(record.measures(i).vecPeakTime(indRight), record.measures(i).vecPeakRate(indRight),'xk'); % 'HandleVisibility','off');
        ylabel('Rate (spk/s)');
        xlabel('Time(s)')
        fixfig;

        vecDeltaPeak_r = [vecDeltaPeak_r record.measures(i).vecPeakTime(4)-record.measures(i).vecPeakTime(2)];
        vecDepth = [vecDepth record.sSelNeuron(i).DepthBelowIntersect];


    else
        continue
    end
end %i


figure;
hold on;
scatter(vecDeltaPeak_l,vecDepth,'r')
scatter(vecDeltaPeak_r,vecDepth','b')
legend('Left', 'Right');

xline(0,'k--','HandleVisibility','off');
ylabel('Depth (um)');
xlabel('Delta peak (s)');
set(gca,'Ydir','reverse');
fixfig;
