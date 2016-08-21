%% rs1050225_MI_clean_LFP.mat

clear all;

% load this data set:
load rs1050225_MI_clean_LFP.mat
%%
% add the following to your path
% addpath(genpath('``server``:\Matt'))

X = zeros(length(lfp1), 96);
for j = 1:96
    try
        X(:,j) = eval(['lfp', num2str(j)]);
    catch
        disp(['no variable named lfp', num2str(j)]);
    end;
end;



reactionTimeUB = .75;
reactionTimeLB = .1;
movementTimeUB = 1.5;
movementTimeLB = .15;

beh = removeSlowReactionTimesFromBEH(beh, reactionTimeUB, reactionTimeLB);
beh = removeSlowMovementTimesFromBEH(beh, movementTimeUB, movementTimeLB);

md = zeros(8,1);
for i = 1:8
    md(i) = quantile(beh(beh(:,8) == i, 6) - beh(beh(:,8) == i, 5), .75);
end;

fastTrials = (beh(:,6) - beh(:,5)) < md(beh(:,8));
beh = beh(fastTrials,:);

numTrials = size(beh, 1);
eventTimes = round(beh(:, 5) * 1000);
%%
H = zeros(2001, numTrials, 96);
for j = 1:96
    Y = zeros(2001, numTrials);
    for i = 1:numTrials
        Y(:,i) = X(eventTimes(i)-1000:eventTimes(i)+1000,j);
    end;
    H(:, :, j) = hilbert(filterData(Y, 18, 3, 1000));
end;

absH = abs(H);
plotBetaAttenuation = squeeze(nanmean(absH,2));
for i = 1:96
    plot(plotBetaAttenuation(:,i))
    pause
end
clear lfp*
%%
[abat, am, aa] = findABAT_local_AT(H, 0.15, 0, 1000, 2);

%%
close all;
figure; hold on;
plotScalarOnUEA((abat(MIchans)-1000)./1000, MIchan2rc(MIchans,:), [], 1)
%plotChanonUEA(MIchan2rc(1:96,:));
[B, r2] = estimateResultant(abat(MIchans), MIchan2rc(MIchans, :));
disp('Rockstar');
mdl = LinearModel.fit(MIchan2rc(MIchans,:), abat(MIchans));
disp(mdl)
axis square; axis equal;
arrowhead([5.5 5.5], [5.5 5.5] + 5.5 * r2 * ([B(2), -B(1)])./norm([B(2), -B(1)]), 'lineWidth', 2, 'color', 'k');
disp(' ');
disp(' ');
disp(' ');

RS1050225MIABAT = abat(MIchans);
plotChanonUEA(MIchan2rc(1:96,:))