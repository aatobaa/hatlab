clear all;

load 'P-rs1050225_M1.mat'
load 'torqueRaw-rs1050225_M1.mat'
load 'rs1050225_MI_clean_LFP.mat'

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

behrs = beh;
behrs = removeSlowReactionTimesFromBEH(behrs, reactionTimeUB, reactionTimeLB);
behrs = removeSlowMovementTimesFromBEH(behrs, movementTimeUB, movementTimeLB);

md = zeros(8,1);
for i = 1:8
    md(i) = quantile(behrs(behrs(:,8) == i, 6) - behrs(behrs(:,8) == i, 5), .75);
end;

fastTrials = (behrs(:,6) - behrs(:,5)) < md(behrs(:,8));
behrs = behrs(fastTrials,:);

numTrials = size(behrs, 1);
eventTimes = round(behrs(:, 5) * 1000);

pjtrs = nan(numTrials, 1);
%time of peak joint torque; relative to movement onset
tpjtrs = nan(numTrials, 1);
torqueOnsetjtrs = nan(numTrials,1);
timeOfTorqueOnsetjtrs = nan(numTrials,1);

for i = 1:numTrials
    % relative to GO CUE
    idxs = find(kin.raw.stamps > behrs(i,4) & kin.raw.stamps < behrs(i,6));
    %shoulder torque
    st = real(torque(idxs,1));
    %elbow torque
    et = real(torque(idxs,2));
    %combined joint torque magnitude as an overall proxy for torque
    %activity
    cjt = sqrt(st.^2 + et.^2);
    
    % Use 0.1 for RS1050225 and 0.03 for V1050917
    [pks, locs] = findpeaks(cjt, 'MinPeakHeight', .1);
    if ~isempty(locs)
        pjtrs(i) = cjt(locs(1));
        % relative to MO
        tpjtrs(i) = kin.raw.stamps(idxs(locs(1))) - behrs(i,5);
        
        %Adil 7/26: torque Onset instead of peak torque
        tO = getTorqueOnset(cjt,2,.2);
        torqueOnsetMagjtrs(i) = cjt(tO);
        timeOfTorqueOnsetjtrs(i) = kin.raw.stamps(idxs(tO))-behrs(i,5);

    else
        disp(num2str(i))
    end;
%     subplot(1,3,1), plot(st)
%     ylabel('shoulder')
%     subplot(1,3,2), plot(et)
%     ylabel('elbow')
%     subplot(1,3,3), plot(cjt)
%     ylabel('combined')
    %saveas(gcf,strcat('example_',num2str(i)),'epsc')
    pause
end;
%%


clearvars -except tpjtrs behrs timeOfTorqueOnsetjtrs


%%
%%
load 'P-V1050913_M1Ab_PMvAc001.mat'
load 'torqueRaw-V1050913_M1Ab_PMvAc001.mat'
%%


reactionTimeUB = .75;
reactionTimeLB = .1;
movementTimeUB = 1.5;
movementTimeLB = .15;

behv13 = makeBEH(conditions, events);
behv13 = [behv13, (1:size(behv13, 1))'];
behv13 = removeSlowReactionTimesFromBEH(behv13, reactionTimeUB, reactionTimeLB);
behv13 = removeSlowMovementTimesFromBEH(behv13, movementTimeUB, movementTimeLB);

md = zeros(8,1);
for i = 1:8
    md(i) = quantile(behv13(behv13(:,8) == i, 6) - behv13(behv13(:,8) == i, 5), .75);
end;

fastTrials = (behv13(:,6) - behv13(:,5)) < md(behv13(:,8));
behv13 = behv13(fastTrials,:);

disp('done');


numTrials = size(behv13, 1);
pjtv13 = nan(numTrials, 1);
tpjtv13 = nan(numTrials, 1);
torqueOnsetjtv13 = nan(numTrials,1);
timeOfTorqueOnsetjtv13 = nan(numTrials,1);
for i = 1:numTrials
    % relative to GO CUE
    idxs = find(kin.raw.stamps > behv13(i,4) & kin.raw.stamps < behv13(i,6));
    st = real(torque(idxs,1));
    et = real(torque(idxs,2));
    
    cjt = sqrt(st.^2 + et.^2);
    
    % Use 0.1 for RS1050225 and 0.03 for V1050917
    [pks, locs] = findpeaks(cjt, 'MinPeakHeight', .03);
    if ~isempty(locs)
        pjtv13(i) = cjt(locs(1));
        % relative to MO
        tpjtv13(i) = kin.raw.stamps(idxs(locs(1))) - behv13(i,5);
        
        %Adil 7/26: torque Onset instead of peak torque
        tO = getTorqueOnset(cjt,2,.2);
        torqueOnsetMagjtv13(i) = cjt(tO);
        timeOfTorqueOnsetjtv13(i) = kin.raw.stamps(idxs(tO))-behv13(i,5);
    else
        disp(num2str(i))
    end;
end;

clearvars -except tpjtrs behrs timeOfTorqueOnsetjtrs behv13 tpjtv13...
    timeOfTorqueOnsetjtv13 



%%
%%
load 'P-V1050917_M1Ab_PMvAc002.mat'
load 'torqueRaw-V1050917_M1Ab_PMvAc002.mat'

%%
reactionTimeUB = .75;
reactionTimeLB = .1;
movementTimeUB = 1.5;
movementTimeLB = .15;

behv17 = makeBEH(conditions, events);
behv17 = [behv17, (1:size(behv17, 1))'];
behv17 = removeSlowReactionTimesFromBEH(behv17, reactionTimeUB, reactionTimeLB);
behv17 = removeSlowMovementTimesFromBEH(behv17, movementTimeUB, movementTimeLB);

md = zeros(8,1);
for i = 1:8
    md(i) = quantile(behv17(behv17(:,8) == i, 6) - behv17(behv17(:,8) == i, 5), .75);
end;

fastTrials = (behv17(:,6) - behv17(:,5)) < md(behv17(:,8));
behv17 = behv17(fastTrials,:);

disp('done');


numTrials = size(behv17, 1);


pjtv17 = nan(numTrials, 1);
tpjtv17 = nan(numTrials, 1);
torqueOnsetjtv17 = nan(numTrials,1);
timeOfTorqueOnsetjtv17 = nan(numTrials,1);

for i = 1:numTrials
    % relative to GO CUE
    idxs = find(kin.raw.stamps > behv17(i,4) & kin.raw.stamps < behv17(i,6));
    st = real(torque(idxs,1));
    et = real(torque(idxs,2));
    
    cjt = sqrt(st.^2 + et.^2);
    
    % Use 0.1 for RS1050225 and 0.03 for V1050917
    [pks, locs] = findpeaks(cjt, 'MinPeakHeight', .03);
    if ~isempty(locs)
        pjtv17(i) = cjt(locs(1));
        % relative to MO
        tpjtv17(i) = kin.raw.stamps(idxs(locs(1))) - behv17(i,5);

        
        %Adil 7/26: torque Onset instead of peak torque
        tO = getTorqueOnset(cjt,2,.2);
        torqueOnsetMagjtv17(i) = cjt(tO);
        timeOfTorqueOnsetjtv17(i) = kin.raw.stamps(idxs(tO))-behv17(i,5);
    else
        disp(num2str(i))
    end;
end;

clearvars -except tpjtrs behrs timeOfTorqueOnsetjtrs behv13 tpjtv13...
    timeOfTorqueOnsetjtv13 behv17 tpjtv17 timeOfTorqueOnsetjtv17
%%
%%
close all;
figure
set(gcf, 'position', [560 720 260 210])
hold on;
plot(reshape([(1:8)' (1:8)', nan(8,1)]' - .15, [], 1), reshape([grpstats(tpjtrs, behrs(:,8), 'meanci'), nan(8,1)]', [], 1), 'color', [1 .5 .6], 'lineWidth', 1.25)
plot(reshape([(1:8)' (1:8)', nan(8,1)]', [], 1), reshape([grpstats(tpjtv13, behv13(:,8), 'meanci'), nan(8,1)]', [], 1), 'color', [.25 .5 1], 'lineWidth', 1.25)
plot(reshape([(1:8)' (1:8)', nan(8,1)]' + .15, [], 1), reshape([grpstats(tpjtv17, behv17(:,8), 'meanci'), nan(8,1)]', [], 1), 'color', [.5 .75 1], 'lineWidth', 1.25)

plot((1:8)-.15, grpstats(tpjtrs, behrs(:,8), 'mean'), 'k.')
plot(1:8, grpstats(tpjtv13, behv13(:,8), 'mean'), 'k.')
plot((1:8)+.15, grpstats(tpjtv17, behv17(:,8), 'mean'), 'k.')

plot(1:8, grpstats([tpjtrs; tpjtv13; tpjtv17], [behrs(:,8); behv13(:,8); behv17(:,8)], 'mean'), 'color', 'k', 'lineWidth', 1.25)

xlim([0.5 8.5])
ylim([-.04 .12])
set(gca, 'xtick', 1:8)

ylabel('peak joint torque (s)')

%% Adil 7/26: Plots using Torque Onset 

close all;
figure
set(gcf, 'position', [560 720 260 210])
hold on;
plot(reshape([(1:8)' (1:8)', nan(8,1)]' - .15, [], 1), reshape([grpstats(timeOfTorqueOnsetjtrs, behrs(:,8), 'meanci'), nan(8,1)]', [], 1), 'color', [1 .5 .6], 'lineWidth', 1.25)
plot(reshape([(1:8)' (1:8)', nan(8,1)]', [], 1), reshape([grpstats(timeOfTorqueOnsetjtv13, behv13(:,8), 'meanci'), nan(8,1)]', [], 1), 'color', [.25 .5 1], 'lineWidth', 1.25)
plot(reshape([(1:8)' (1:8)', nan(8,1)]' + .15, [], 1), reshape([grpstats(timeOfTorqueOnsetjtv17, behv17(:,8), 'meanci'), nan(8,1)]', [], 1), 'color', [.5 .75 1], 'lineWidth', 1.25)

plot((1:8)-.15, grpstats(timeOfTorqueOnsetjtrs, behrs(:,8), 'mean'), 'k.')
plot(1:8, grpstats(timeOfTorqueOnsetjtv13, behv13(:,8), 'mean'), 'k.')
plot((1:8)+.15, grpstats(timeOfTorqueOnsetjtv17, behv17(:,8), 'mean'), 'k.')

plot(1:8, grpstats([timeOfTorqueOnsetjtrs; timeOfTorqueOnsetjtv13; timeOfTorqueOnsetjtv17], [behrs(:,8); behv13(:,8); behv17(:,8)], 'mean'), 'color', 'k', 'lineWidth', 1.25)

set(gca, 'xtick', 1:8)

xlim([0.5 8.5])
ylim([-.26 -.09])
ylabel('time of torque onset (s)')
xlabel('movement direction')
title('Torque Onset vs movement directions for Rockstar, V13 and V17')

%% Rockstar only
close all;
figure
set(gcf, 'position', [560 720 260 210])
hold on;
plot(reshape([(1:8)' (1:8)', nan(8,1)]' - .15, [], 1), reshape([grpstats(timeOfTorqueOnsetjtrs, behrs(:,8), 'meanci'), nan(8,1)]', [], 1), 'color', [1 .5 .6], 'lineWidth', 1.25)

plot((1:8)-.15, grpstats(timeOfTorqueOnsetjtrs, behrs(:,8), 'mean'), 'k.')

plot(1:8, grpstats([timeOfTorqueOnsetjtrs], [behrs(:,8)], 'mean'), 'color', 'k', 'lineWidth', 1.25)

set(gca, 'xtick', 1:8)

%xlim([0.5 8.5])
%ylim([-.26 -.09])
ylabel('time of torque onset (s)')

%% Fit line through peak torque vs torque onset
linScale = -.1:0.01:.35;
tpjtrs_nansremoved_Indx = ~isnan(tpjtrs)
timeOfTorqueOnsetjtrs_nansremoved_Indx = ~isnan(timeOfTorqueOnsetjtrs);

validIndx = tpjtrs_nansremoved_Indx & timeOfTorqueOnsetjtrs_nansremoved_Indx;

tpjtrs_nansremoved= tpjtrs(validIndx) 
timeOfTorqueOnsetjtrs_nansremoved = timeOfTorqueOnsetjtrs(validIndx) 

p = polyfit(tpjtrs_nansremoved, timeOfTorqueOnsetjtrs_nansremoved,1);
yfit = polyval(p,tpjtrs_nansremoved);
yresid = timeOfTorqueOnsetjtrs_nansremoved - yfit;
SSresid = sum(yresid.^2);

SStotal = sum((timeOfTorqueOnsetjtrs_nansremoved - mean(timeOfTorqueOnsetjtrs_nansremoved)).^2);

%SStotal = (length(timeOfTorqueOnsetjtrs_nansremoved)-1 * var(timeOfTorqueOnsetjtrs_nansremoved));
rsq = 1 - (SSresid/SStotal);
rsq
scatter(tpjtrs,timeOfTorqueOnsetjtrs)
hold on
plot(linScale, p(1) * linScale + p(2));
text(max(xlim)*.8,max(ylim)*.8,strcat('Rsq = ',num2str(rsq)),'FontSize',14)
xlabel('Time of Peak Joint Torque')
ylabel('Time of Torque Onset')
title('Time of Torque Onset vs Peak Torque, relative to Movement Onset (Rockstar)')