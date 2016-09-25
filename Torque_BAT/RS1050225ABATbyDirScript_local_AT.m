%% rs1050225_MI_clean_LFP.mat

%clear all;

% load this data set:
load rs1050225_MI_clean_LFP.mat

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
%changing from movement onset to torque onset
% eventTimes = round(xyz * 1000);
% aligned to movement onset 
eventTimes = round(beh(:, 5) * 1000);


H = zeros(2001, numTrials, 96);
for j = 1:96
    Y = zeros(2001, numTrials);
    for i = 1:numTrials
        if ~isnan(eventTimes(i))
            Y(:,i) = X(eventTimes(i)-1000:eventTimes(i)+1000,j);
        else
            Y(:,i) = nan;
        end
    end;
    H(:, :, j) = hilbert(filterData(Y, 18, 3, 1000));
end;

clear lfp*
%% EXPERIMENTAL ANALYSES
%%Looking at plots of beta averaged across electrodes 
hold on 
for nn = 1:96
    plot(abs(H(:,12,nn)))
    pause
end
    plot(nanmean(squeeze(abs(H(:,12,:))),2),'r','LineWidth',12)
%% Creating a beta_profile, downsampled to match the sampling frequency of 
% torque profiles. Then compute cross-correlation
beta_profiles = downsample(abs(nanmean(H,3)),2);
for trial = 1:numTrials
    crosscorr(beta_profiles(:,trial),torque_profiles(:,trial),999)
    pause
end
%% Single trials are now filtered 

torque_profiles = cjt_torque_profiles;
singleTrialBATs = zeros(numTrials,1);
for trial = 1:numTrials
%     trial = 1;
    % Normalize profiles
    minBeta = min(beta_profiles(:,trial));
    maxBeta= max(beta_profiles(:,trial));
    minTorque = min(torque_profiles(:,trial));
    maxTorque = max(torque_profiles(:,trial));

    a=-1; b=1;
    % Using formula from http://www.mathworks.com/matlabcentral/fileexchange/5103-toolbox-diffc/content/toolbox_diffc/toolbox/rescale.m
    normBeta = (b-a) .* (beta_profiles(:,trial) - minBeta)./(maxBeta-minBeta) + a;
    normTorque = (b-a) .* (torque_profiles(:,trial) - minTorque)./(maxTorque-minTorque) + a; 

        
    % try filtering
    fc = 5;
    fs = 10;
    wc = 0.005;
    [b,a] = butter(4,wc);
    filteredBeta = filtfilt(b,a,normBeta);
%     
%     plot(normBeta,'r','LineWidth',2)
%     hold on
%     plot(filteredBeta, 'b','LineWidth',2)
% %     plot(normTorque,'b','LineWidth',2)
%     hold off
%     legend('show','Beta','Torque')
%     
%     pause
    
    [singleTrialBAT, ~, ~, ~] = findABAT_local_AT(filteredBeta,.15,0,500,7);
    singleTrialBATs(trial) = singleTrialBAT;
%     pause
%     plot(beta_profiles(:,trial),'r')
%     hold on
%     plot(torque_profiles(:,trial),'b')
%     hold off
end
%% Expanding the above analysis with simulations 

torque_profiles = cjt_torque_profiles;
for trial = 1:1
    % Normalize profiles
    minBeta = min(beta_profiles(:,trial));
    maxBeta= max(beta_profiles(:,trial));
    minTorque = min(torque_profiles(:,trial));
    maxTorque = max(torque_profiles(:,trial));
    a=-1; b=1;
    % Using formula from http://www.mathworks.com/matlabcentral/fileexchange/5103-toolbox-diffc/content/toolbox_diffc/toolbox/rescale.m
    normBeta = (b-a) .* (beta_profiles(:,trial) - minBeta)./(maxBeta-minBeta) + a;
    normTorque = (b-a) .* (torque_profiles(:,trial) - minTorque)./(maxTorque-minTorque) + a; 
%     for X_OFFSET = 540:560
        for X_SCALE = 1:0.1:2
            xlinspace = X_SCALE*(1:1001);
            corrcoef(normTorque
            
            
            
            plot([normBeta(1:end-X_OFFSET); zeros(X_OFFSET,1)],'r','LineWidth',2)
            hold on
            plot(X_SCALE,[normTorque(X_OFFSET:end); zeros(X_OFFSET,1)],'b','LineWidth',2)
            hold off
            legend('show','Beta','Torque')
            pause
        %     plot(beta_profiles(:,trial),'r')
        %     hold on
        %     plot(torque_profiles(:,trial),'b')
        %     hold off
        end
%     end
end

%% 8 directions
fastTrials = (1:size(beh,1))';
method = 2;
[abatLeft, amLeft, aaLeft, adaLeft] = findABAT_local_AT(...
    H(:, fastTrials & (beh(:,8) == 5 | beh(:,8) == 5 | beh(:,8) == 5), :), ...
    .15, 0, 1000, method);

[abatRight, amRight, aaRight, adaRight] = findABAT_local_AT(...
    H(:, fastTrials & (beh(:,8) == 1 | beh(:,8) == 1 | beh(:,8) == 1), :), ...
    .15, 0, 1000, method);

[abatTop, amTop, aaTop, adaTop]    = findABAT_local_AT(...
    H(:, fastTrials & (beh(:,8) == 3 | beh(:,8) == 3 | beh(:,8) == 3), :), ...
    .15, 0, 1000, method);

[abatBottom, amBottom, aaBottom, adaBottom] = findABAT_local_AT(...
    H(:, fastTrials & (beh(:,8) == 7 | beh(:,8) == 7 | beh(:,8) == 7), :), ...
    .15, 0, 1000, method);

[abatNE, amNE, aaNE, adaNE] = findABAT_local_AT(...
    H(:, fastTrials & (beh(:,8) == 2 | beh(:,8) == 2 | beh(:,8) == 2), :), ...
    .15, 0, 1000, method);

[abatNW, amNW, aaNW, adaNW] = findABAT_local_AT(...
    H(:, fastTrials & (beh(:,8) == 4 | beh(:,8) == 4 | beh(:,8) == 4), :), ...
    .15, 0, 1000, method);

[abatSW, amSW, aaSW, adaSW] = findABAT_local_AT(...
    H(:, fastTrials & (beh(:,8) == 6 | beh(:,8) == 6 | beh(:,8) == 6), :), ...
    .15, 0, 1000, method);

[abatSE, amSE, aaSE, adaSE] = findABAT_local_AT(...
    H(:, fastTrials & (beh(:,8) == 8 | beh(:,8) == 8 | beh(:,8) == 8), :), ...
    .15, 0, 1000, method);
%}
%% blocked directions
% Movement directions
%      3
%    4   2  
%  5   x   1
%    6   8 
%      7

% -L-+-+-R-
%    |3|
%   4| |2  
% 5  |x|  1
%   6| |8 
%    |7|

% T    3
% |  4   2  
% +---------
% |5   x   1
% +---------
% |  6   8 
% B    7
%fastTrials = (beh(:,6) - beh(:,5)) < median(beh(:,6) - beh(:,5));
fastTrials = (1:size(beh,1))';

[abatLeft, amLeft, aaLeft, adaLeft] = findABAT_local_AT(...
    H(:, fastTrials & (beh(:,8) == 4 | beh(:,8) == 5 | beh(:,8) == 6), :), ...
    .15, 0, 1000, 2);

[abatRight, amRight, aaRight, adaRight] = findABAT_local_AT(...
    H(:, fastTrials & (beh(:,8) == 2 | beh(:,8) == 1 | beh(:,8) == 8), :), ...
    .15, 0, 1000, 2);

[abatTop, amTop, aaTop, adaTop]    = findABAT_local_AT(...
    H(:, fastTrials & (beh(:,8) == 4 | beh(:,8) == 3 | beh(:,8) == 2), :), ...
    .15, 0, 1000, 2);

[abatBottom, amBottom, aaBottom, adaBottom] = findABAT_local_AT(...
    H(:, fastTrials & (beh(:,8) == 8 | beh(:,8) == 7 | beh(:,8) == 6), :), ...
    .15, 0, 1000, 2);

[abatNE, amNE, aaNE, adaNE] = findABAT_local_AT(...
    H(:, fastTrials & (beh(:,8) == 3 | beh(:,8) == 2 | beh(:,8) == 1), :), ...
    .15, 0, 1000, 2);

[abatNW, amNW, aaNW, adaNW] = findABAT_local_AT(...
    H(:, fastTrials & (beh(:,8) == 3 | beh(:,8) == 4 | beh(:,8) == 5), :), ...
    .15, 0, 1000, 2);

[abatSW, amSW, aaSW, adaSW] = findABAT_local_AT(...
    H(:, fastTrials & (beh(:,8) == 5 | beh(:,8) == 7 | beh(:,8) == 6), :), ...
    .15, 0, 1000, 2);

[abatSE, amSE, aaSE, adaSE] = findABAT_local_AT(...
    H(:, fastTrials & (beh(:,8) == 8 | beh(:,8) == 7 | beh(:,8) == 1), :), ...
    .15, 0, 1000, 2);

%%
%close all;

minr = min([abatTop(MIchans); abatBottom(MIchans); abatLeft(MIchans); abatRight(MIchans); ...
            abatNE(MIchans);  abatNW(MIchans);     abatSW(MIchans);   abatSE(MIchans)]);
maxr = max([abatTop(MIchans); abatBottom(MIchans); abatLeft(MIchans); abatRight(MIchans); ...
            abatNE(MIchans);  abatNW(MIchans);     abatSW(MIchans);   abatSE(MIchans)]);
figure;
set(gcf, 'position', [360 280 760 650])
% top
subplot(3, 3, 2);
hold on;
plotScalarOnUEA((abatTop(MIchans)-1000)./1000, MIchan2rc(MIchans,:), ([minr maxr]-1000)./1000, 0, 'markerSize', 16);
[topB r2] = estimateResultant(abatTop(MIchans), MIchan2rc(MIchans, :));
%arrowhead([5.5 5.5], [topB(2) + 5.5, -topB(1) + 5.5], 'lineWidth', 2, 'color', 'k');
arrowhead([5.5 5.5], [5.5 5.5] + 5.5*r2*[topB(2) -topB(1)] ./ norm([topB(2) -topB(1)]), 'lineWidth', 2, 'color', 'k');

% left
subplot(3, 3, 4)
plotScalarOnUEA((abatLeft(MIchans)-1000)./1000, MIchan2rc(MIchans,:), ([minr maxr]-1000)./1000, 0, 'markerSize', 16);
[leftB r2] = estimateResultant(abatLeft(MIchans), MIchan2rc(MIchans, :));
%arrowhead([5.5 5.5], [leftB(2) + 5.5, -leftB(1) + 5.5], 'lineWidth', 2, 'color', 'k');
arrowhead([5.5 5.5], [5.5 5.5] + 5.5*r2*[leftB(2) -leftB(1)] ./ norm([leftB(2) -leftB(1)]), 'lineWidth', 2, 'color', 'k')

% bottom
subplot(3, 3, 8)
plotScalarOnUEA((abatBottom(MIchans)-1000)./1000, MIchan2rc(MIchans,:), ([minr maxr]-1000)./1000, 0, 'markerSize', 16);
[bottomB r2] = estimateResultant(abatBottom(MIchans), MIchan2rc(MIchans, :));
%arrowhead([5.5 5.5], [bottomB(2) + 5.5, -bottomB(1) + 5.5], 'lineWidth', 2, 'color', 'k');
arrowhead([5.5 5.5], [5.5 5.5] + 5.5*r2*[bottomB(2) -bottomB(1)] ./ norm([bottomB(2) -bottomB(1)]), 'lineWidth', 2, 'color', 'k');

% right
subplot(3, 3, 6)
plotScalarOnUEA((abatRight(MIchans)-1000)./1000, MIchan2rc(MIchans,:), ([minr maxr]-1000)./1000, 0, 'markerSize', 16);
[rightB r2] = estimateResultant(abatRight(MIchans), MIchan2rc(MIchans, :));
%arrowhead([5.5 5.5], [rightB(2) + 5.5, -rightB(1) + 5.5], 'lineWidth', 2, 'color', 'k');
arrowhead([5.5 5.5], [5.5 5.5] + 5.5*r2*[rightB(2) -rightB(1)] ./ norm([rightB(2) -rightB(1)]), 'lineWidth', 2, 'color', 'k');

% NE
subplot(3, 3, 3);
hold on;
plotScalarOnUEA((abatNE(MIchans)-1000)./1000, MIchan2rc(MIchans,:), ([minr maxr]-1000)./1000, 0, 'markerSize', 16);
[NEB r2] = estimateResultant(abatNE(MIchans), MIchan2rc(MIchans, :));
%arrowhead([5.5 5.5], [NEB(2) + 5.5, -NEB(1) + 5.5], 'lineWidth', 2, 'color', 'k');
arrowhead([5.5 5.5], [5.5 5.5] + 5.5*r2*[NEB(2) -NEB(1)] ./ norm([NEB(2) -NEB(1)]), 'lineWidth', 2, 'color', 'k');

% NW
subplot(3, 3, 1);
hold on;
plotScalarOnUEA((abatNW(MIchans)-1000)./1000, MIchan2rc(MIchans,:), ([minr maxr]-1000)./1000, 0, 'markerSize', 16);
[NWB r2] = estimateResultant(abatNW(MIchans), MIchan2rc(MIchans, :));
%arrowhead([5.5 5.5], [NWB(2) + 5.5, -NWB(1) + 5.5], 'lineWidth', 2, 'color', 'k');
arrowhead([5.5 5.5], [5.5 5.5] + 5.5*r2*[NWB(2) -NWB(1)] ./ norm([NWB(2) -NWB(1)]), 'lineWidth', 2, 'color', 'k');

% SW
subplot(3, 3, 7);
hold on;
plotScalarOnUEA((abatSW(MIchans)-1000)./1000, MIchan2rc(MIchans,:), ([minr maxr]-1000)./1000, 0, 'markerSize', 16);
[SWB r2] = estimateResultant(abatSW(MIchans), MIchan2rc(MIchans, :));
%arrowhead([5.5 5.5], [SWB(2) + 5.5, -SWB(1) + 5.5], 'lineWidth', 2, 'color', 'k');
arrowhead([5.5 5.5], [5.5 5.5] + 5.5*r2*[SWB(2) -SWB(1)] ./ norm([SWB(2) -SWB(1)]), 'lineWidth', 2, 'color', 'k');

% SE
subplot(3, 3, 9);
hold on;
plotScalarOnUEA((abatSE(MIchans)-1000)./1000, MIchan2rc(MIchans,:), ([minr maxr]-1000)./1000, 0, 'markerSize', 16);
[SEB r2] = estimateResultant(abatSE(MIchans), MIchan2rc(MIchans, :));
%arrowhead([5.5 5.5], [SEB(2) + 5.5, -SEB(1) + 5.5], 'lineWidth', 2, 'color', 'k');
arrowhead([5.5 5.5], [5.5 5.5] + 5.5*r2*[SEB(2) -SEB(1)] ./ norm([SEB(2) -SEB(1)]), 'lineWidth', 2, 'color', 'k');


% scale bar
subplot(3, 3, 5)
plotScalarOnUEA([], [], ([minr, maxr]-1000)./1000)

subplotHandles = get(gcf, 'children');
subplotHandles = subplotHandles(2:end); % first entry is the scale bar
%
%R = eye(2);
r = pi/2;
%R = [cos(r) sin(r); -sin(r) cos(r)]; % RS
R = [cosd(55) sind(55); -sind(55) cosd(55)]; %RX
%R = [cosd(3) sind(3); -sind(3) cosd(3)]; % Boo
%R = [cosd(55) sind(55); -sind(55) cosd(55)]; % V
for f = 1:length(subplotHandles)
    subplotChildren = get(subplotHandles(f), 'children');
    % apply rotation to axes
    xl = get(subplotHandles(f), 'XLim');
    yl = get(subplotHandles(f), 'YLim');
    xl = xl(:);
    yl = yl(:);
    % Typicaly, we specify the axis as [xMin, xMax, yMin, yMax].  This is
    % because x components and y components are independent.  If we want to
    % rotate our axis limits, we want to rotate points.  Specifically, 
    %  
    %   (xMin, yMax)                         (xMax, yMax)
    %              +-------------------------+
    %              |                         |
    %              |                         |
    %              +-------------------------+
    %   (xMin, yMin)                         (xMax, yMin) 
    % 
    % We then want to find our new axis limits, but still in the same xy
    % space.
    
    M = [xl, yl; xl(2:-1:1), yl] * R;
    
    set(subplotHandles(f), 'xLim', [min(M(:,1)), max(M(:,1))], ...
                           'yLim', [min(M(:,2)), max(M(:,2))]);
    
    
    % apply rotation to data
    for s = 1:length(subplotChildren)
        xData = get(subplotChildren(s), 'XData');
        yData = get(subplotChildren(s), 'YData');
        
        M = [xData(:), yData(:)] * R;
        
        set(subplotChildren(s), 'xData', M(:,1), 'yData', M(:,2));
        
    end;
    
end;
set(gcf, 'paperPosition', [0 0 10 10]);
%%

% test the hypothesis that blocked ABAT vectors are different by comparing
% the magnitude of the resultant difference to a null distribution.  



% Group A: left
% Group B: right

numTrialsA = sum(fastTrials & (beh(:,8) == 4 | beh(:,8) == 4 | beh(:,8) == 5));
numTrialsB = sum(fastTrials & (beh(:,8) == 2 | beh(:,8) == 8 | beh(:,8) == 1));
trialIdxs = find(fastTrials & beh(:,8) ~= 3 & beh(:,8) ~= 7);

if numTrialsA + numTrialsB ~= length(trialIdxs)
    disp('inconsistency between numTrials and trialIDXs');
end;

numBoot = 1000;
bsResultantMags = zeros(numBoot, 1);
bsBA = zeros(numBoot, 2);
bsBB = zeros(numBoot, 2);
tic
for b = 1:numBoot
    rIdxs = trialIdxs(randsample(1:(numTrialsA + numTrialsB), numTrialsA + numTrialsB));
    AIdxs = rIdxs(1:numTrialsA);
    BIdxs = setdiff(rIdxs, AIdxs);
    [abatA, ~, ~] = findABAT(H(:, AIdxs, :), .15, 0, 1000, 2);
    [abatB, ~, ~] = findABAT(H(:, BIdxs, :), .15, 0, 1000, 2);
    
    BA = estimateResultant(abatA(MIchans), MIchan2rc(MIchans, :));
    BB = estimateResultant(abatB(MIchans), MIchan2rc(MIchans, :));
    
    bsBA(b,:) = BA;
    bsBB(b,:) = BB;
    
    bsResultantMags(b) = norm(BA - BB);
    if ~mod(b, 100)
        toc
    end;
end;

%%
figure; hold on;
for b = 1:numBoot
    plot([0, (bsBA(b,2) - bsBB(b,2))], [0, -(bsBA(b,1) - bsBB(b,1))], 'color', [.5 .5 .5], 'lineWidth', 2)
end;

% actual resultant
plot([0, (leftB(2) - rightB(2))], [0, -(leftB(1) - rightB(1))], 'color', [.5 0 0], 'lineWidth', 2)
axis equal;
axis square;
set(gca, 'fontSize', 14)
%title('top bottom')
title('left right')

%%
% test the hypothesis that the two mean vectors are different by
% bootstrapping distribution of means from each population
% Movement directions
%      3
%    4   2  
%  5   x   1
%    6   8 
%      7

trialIdxsA = find(fastTrials & (beh(:,8) == 3 | beh(:,8) == 2 | beh(:,8) == 4));
trialIdxsB = find(fastTrials & (beh(:,8) == 8 | beh(:,8) == 6 | beh(:,8) == 7));
trialIdxs = find(fastTrials & beh(:,8) ~= 5 & beh(:,8) ~= 1);

numTrialsA = numel(trialIdxsA);
numTrialsB = numel(trialIdxsB);

if numTrialsA + numTrialsB ~= length(trialIdxs)
    disp('inconsistency between numTrials and trialIDXs');
end;

numBoot = 1000;

bsBA = zeros(numBoot, 2);
bsBB = zeros(numBoot, 2);
tic
for b = 1:numBoot
    AIdxs = trialIdxsA(randsample(1:numTrialsA, 50));
    BIdxs = trialIdxsB(randsample(1:numTrialsB, 50));
    
    [abatA, ~, ~] = findABAT(H(:, AIdxs, :), .15, 0, 1000, 2);
    [abatB, ~, ~] = findABAT(H(:, BIdxs, :), .15, 0, 1000, 2);
    
    BA = estimateResultant(abatA(MIchans), MIchan2rc(MIchans, :));
    BB = estimateResultant(abatB(MIchans), MIchan2rc(MIchans, :));
    
    bsBA(b,:) = BA;
    bsBB(b,:) = BB;
    
    bsResultantMags(b) = norm(BA - BB);
    if ~mod(b, 100)
        toc
    end;
end;


%%
%close all;
figure; hold on;
for b = 1:numBoot
    plot([0, bsBA(b,2)], [0, -bsBA(b,1)], 'color', 'r', 'lineWidth', 2)
end;
axis equal;

figure; hold on;
for b = 1:numBoot
    plot([0, bsBB(b,2)], [0, -bsBB(b,1)], 'color', 'g', 'lineWidth', 2);
end;
axis equal;
    
anglesBetweenConditions = acos( (dot(bsBA, bsBB, 2)) ./ ...
    (sqrt(dot(bsBA, bsBA, 2)) .* sqrt(dot(bsBB, bsBB, 2)))) * 180/pi;

observedDifference = acos((topB'*bottomB) ./ (norm(topB) * norm(bottomB))) * 180/pi;

%%
% Generate many bootstrapped samples from each blocked directional
% distribution and plot the result (relies on variables created in previous
% cells)

% Movement directions
%      3
%    4   2  
%  5   x   1
%    6   8 
%      7

close all;
figure; hold on;

topIdxs = find(fastTrials & (beh(:,8) == 2 | beh(:,8) == 3 | beh(:,8) == 4));
bottomIdxs = find(fastTrials & (beh(:,8) == 6 | beh(:,8) == 7 | beh(:,8) == 8));
leftIdxs = find(fastTrials & (beh(:,8) == 4 | beh(:,8) == 5 | beh(:,8) == 6));
rightIdxs = find(fastTrials & (beh(:,8) == 8 | beh(:,8) == 1 | beh(:,8) == 2));
NEIdxs = find(fastTrials & (beh(:,8) == 1 | beh(:,8) == 2 | beh(:,8) == 3));
NWIdxs = find(fastTrials & (beh(:,8) == 3 | beh(:,8) == 4 | beh(:,8) == 5));
SWIdxs = find(fastTrials & (beh(:,8) == 5 | beh(:,8) == 6 | beh(:,8) == 7));
SEIdxs = find(fastTrials & (beh(:,8) == 7 | beh(:,8) == 8 | beh(:,8) == 1));


numBoot = 1000;


bsTopB    = zeros(numBoot, 2);
bsBottomB = zeros(numBoot, 2);
bsLeftB   = zeros(numBoot, 2);
bsRightB  = zeros(numBoot, 2);
bsNEB     = zeros(numBoot, 2);
bsNWB     = zeros(numBoot, 2);
bsSEB     = zeros(numBoot, 2);
bsSWB     = zeros(numBoot, 2);

bsTopr2    = zeros(numBoot, 1);
bsBottomr2 = zeros(numBoot, 1);
bsLeftr2   = zeros(numBoot, 1);
bsRightr2  = zeros(numBoot, 1);
bsNEr2     = zeros(numBoot, 1);
bsNWr2     = zeros(numBoot, 1);
bsSEr2     = zeros(numBoot, 1);
bsSWr2     = zeros(numBoot, 1);

bsTopT    = zeros(numBoot, 1);
bsBottomT = zeros(numBoot, 1);
bsLeftT   = zeros(numBoot, 1);
bsRightT  = zeros(numBoot, 1);
bsNET     = zeros(numBoot, 1);
bsNWT     = zeros(numBoot, 1);
bsSET     = zeros(numBoot, 1);
bsSWT     = zeros(numBoot, 1);

absH = abs(H);
for b = 1:numBoot
    tic
    % create random indices
    ai = randsample(1:length(topIdxs), 50, 1);
    bi = randsample(1:length(bottomIdxs), 50, 1);
    ci = randsample(1:length(leftIdxs), 50, 1);
    di = randsample(1:length(rightIdxs), 50, 1);
    ei = randsample(1:length(NEIdxs), 50, 1);
    fi = randsample(1:length(NWIdxs), 50, 1);
    gi = randsample(1:length(SEIdxs), 50, 1);
    hi = randsample(1:length(SWIdxs), 50, 1);
    
    
    % find ABAT for each sample
    bsABATTop    = findABATbs(absH(:, topIdxs(ai), :), .15, 0, 1000, 2);
    bsABATBottom = findABATbs(absH(:, bottomIdxs(bi), :), .15, 0, 1000, 2);
    bsABATLeft   = findABATbs(absH(:, leftIdxs(ci), :), .15, 0, 1000, 2);
    bsABATRight  = findABATbs(absH(:, rightIdxs(di), :), .15, 0, 1000, 2);
    bsABATNE     = findABATbs(absH(:, NEIdxs(ei), :), .15, 0, 1000, 2);
    bsABATNW     = findABATbs(absH(:, NWIdxs(fi), :), .15, 0, 1000, 2);
    bsABATSE     = findABATbs(absH(:, SEIdxs(gi), :), .15, 0, 1000, 2);
    bsABATSW     = findABATbs(absH(:, SWIdxs(hi), :), .15, 0, 1000, 2);
    
    % estimate resultants
    [bsTopB(b,:) bsTopr2(b)]        = estimateResultant(bsABATTop(MIchans), MIchan2rc(MIchans, :));
    [bsBottomB(b,:) bsBottomr2(b)]  = estimateResultant(bsABATBottom(MIchans), MIchan2rc(MIchans, :));
    [bsLeftB(b,:) bsLeftr2(b)]      = estimateResultant(bsABATLeft(MIchans), MIchan2rc(MIchans, :));
    [bsRightB(b,:) bsRightr2(b)]    = estimateResultant(bsABATRight(MIchans), MIchan2rc(MIchans, :));
    [bsNEB(b,:) bsNEr2(b)]          = estimateResultant(bsABATNE(MIchans), MIchan2rc(MIchans, :));
    [bsNWB(b,:) bsNWr2(b)]          = estimateResultant(bsABATNW(MIchans), MIchan2rc(MIchans, :));
    [bsSEB(b,:) bsSEr2(b)]          = estimateResultant(bsABATSE(MIchans), MIchan2rc(MIchans, :));
    [bsSWB(b,:) bsSWr2(b)]          = estimateResultant(bsABATSW(MIchans), MIchan2rc(MIchans, :));
    
    % timing differences
    bsTopT(b)    = median(bsABATTop(MIchans));
    bsBottomT(b) = median(bsABATBottom(MIchans));
    bsLeftT(b)   = median(bsABATLeft(MIchans));
    bsRightT(b)  = median(bsABATRight(MIchans));
    bsNET(b)     = median(bsABATNE(MIchans));
    bsNWT(b)     = median(bsABATNW(MIchans));
    bsSET(b)     = median(bsABATSE(MIchans));
    bsSWT(b)     = median(bsABATSW(MIchans));
    toc
end;
disp('done');
%%
% raw data
close all;
subplot(3, 3, 1)
hold on;
for b = 1:numBoot
    plot([0, bsNWB(b,2)], [0, -bsNWB(b,1)], 'color', 'k', 'lineWidth', 1);
end;
axis(20 * [-1 1 -1 1]);
axis equal;

subplot(3, 3, 2)
hold on;
for b = 1:numBoot
    plot([0, bsTopB(b,2)], [0, -bsTopB(b,1)], 'color', 'k', 'lineWidth', 1);
end;
axis(20 * [-1 1 -1 1]);
axis equal;

subplot(3, 3, 3)
hold on;
for b = 1:numBoot
    plot([0, bsNEB(b,2)], [0, -bsNEB(b,1)], 'color', 'k', 'lineWidth', 1);
end;
axis(20 * [-1 1 -1 1]);
axis equal;

subplot(3, 3, 4)
hold on;
for b = 1:numBoot
    plot([0, bsLeftB(b,2)], [0, -bsLeftB(b,1)], 'color', 'k', 'lineWidth', 1);
end;
axis(20 * [-1 1 -1 1]);
axis equal;

subplot(3, 3, 6)
hold on;
for b = 1:numBoot
    plot([0, bsRightB(b,2)], [0, -bsRightB(b,1)], 'color', 'k', 'lineWidth', 1);
end;
axis(20 * [-1 1 -1 1]);
axis equal;

subplot(3, 3, 7)
hold on;
for b = 1:numBoot
    plot([0, bsSWB(b,2)], [0, -bsSWB(b,1)], 'color', 'k', 'lineWidth', 1);
end;
axis(20 * [-1 1 -1 1]);
axis equal;

subplot(3, 3, 8)
hold on;
for b = 1:numBoot
    plot([0, bsBottomB(b,2)], [0, -bsBottomB(b,1)], 'color', 'k', 'lineWidth', 1);  
end;
axis(20 * [-1 1 -1 1]);
axis equal;

subplot(3, 3, 9)
hold on;
for b = 1:numBoot
    plot([0, bsSEB(b,2)], [0, -bsSEB(b,1)], 'color', 'k', 'lineWidth', 1);
end;
axis(20 * [-1 1 -1 1]);
axis equal;

%%
% plot the distributions of bootstrapped beta attenuation directions.  The
% area of each bin is proportional to the value of each bin; not radius.

close all;
figure; 
maxVal = sqrt(400);
numWheels = 2;
numSpokes = 0;
edges = linspace(0, 2*pi, 21);


subplot(3, 3, 1)
[bsNWR, ~, edges] = rosecw(cart2pol(bsNWB(:,2), -bsNWB(:,1)), edges);
rosePlot(sqrt(bsNWR), edges, 'lineWidth', 2);
roseAxis(numWheels, numSpokes, maxVal);

subplot(3, 3, 2)
[bsTopR, ~, edges] = rosecw(cart2pol(bsTopB(:,2), -bsTopB(:,1)), edges);
rosePlot(sqrt(bsTopR), edges, 'lineWidth', 2);
roseAxis(numWheels, numSpokes, maxVal);

subplot(3, 3, 3)
[bsNER, ~, edges] = rosecw(cart2pol(bsNEB(:,2), -bsNEB(:,1)), edges);
rosePlot(sqrt(bsNER), edges, 'lineWidth', 2);
roseAxis(numWheels, numSpokes, maxVal);

subplot(3, 3, 4)
[bsLeftR, ~, edges] = rosecw(cart2pol(bsLeftB(:,2), -bsLeftB(:,1)), edges);
rosePlot(sqrt(bsLeftR), edges, 'lineWidth', 2);
roseAxis(numWheels, numSpokes, maxVal);

subplot(3, 3, 6)
[bsRightR, ~, edges] = rosecw(cart2pol(bsRightB(:,2), -bsRightB(:,1)), edges);
rosePlot(sqrt(bsRightR), edges, 'lineWidth', 2);
roseAxis(numWheels, numSpokes, maxVal);

subplot(3, 3, 7)
[bsSWR, ~, edges] = rosecw(cart2pol(bsSWB(:,2), -bsSWB(:,1)), edges);
rosePlot(sqrt(bsSWR), edges, 'lineWidth', 2);
roseAxis(numWheels, numSpokes, maxVal);

subplot(3, 3, 8)
[bsBottomR, ~, edges] = rosecw(cart2pol(bsBottomB(:,2), -bsBottomB(:,1)), edges);
rosePlot(sqrt(bsBottomR), edges, 'lineWidth', 2);
roseAxis(numWheels, numSpokes, maxVal);

subplot(3, 3, 9)
[bsSER, ~, edges] = rosecw(cart2pol(bsSEB(:,2), -bsSEB(:,1)), edges);
rosePlot(sqrt(bsSER), edges, 'lineWidth', 2);
roseAxis(numWheels, numSpokes, maxVal);

%%
close all;
figure; hold on;

cmap = hsv(8);
h = .1;


[~, t] = circ_kde(cart2pol(bsRightB(:,2), -bsRightB(:,1)), h);
plot(.125*cos(t), .125*sin(t), 'color', 'k', 'lineWidth', 1);
plot(.25*cos(t), .25*sin(t), 'color', 'k', 'lineWidth', 1);
%plot(5*cos(t), 5*sin(t), 'color', 'k', 'lineWidth', 2);
%plot(7*cos(t), 7*sin(t), 'color', 'k', 'lineWidth', 2);
%plot(10*cos(t), 10*sin(t), 'color', 'k', 'lineWidth', 2);


%text(1.25*cos(3/2), 1.25*sin(3/2), '10%', 'color', 'k')
%text(3.25*cos(3/2), 3.25*sin(3/2), '30%', 'color', 'k')
%text(5.25*cos(3/2), 5.25*sin(3/2), '50%', 'color', 'k')
%text(7.25*cos(3/2), 7.25*sin(3/2), '70%', 'color', 'k')
%text(10.25*cos(3/2), 10.25*sin(3/2), '100%', 'color', 'k')


[y, t] = circ_kde(cart2pol(bsRightB(:,2), -bsRightB(:,1)), h);
y = [y y(1)]; t = [t t(1)];
plot(10*(cos(t) .* y), 10*(sin(t) .* y), 'lineWidth', 3, 'color', white)
plot(10*(cos(t) .* y), 10*(sin(t) .* y), 'lineWidth', 1.5, 'color', lighter(cmap(1,:)))
[an ul ll] = circ_mean(cart2pol(bsRightB(:,2), -bsRightB(:,1)), bsRightr2(:));
%cid = circ_std(cart2pol(bsNEB(:,2), -bsNEB(:,1)), bsNEr2(:));
%line(.2*[0 cos(an)], .2*[0 sin(an)], 'lineWidth', 2, 'color', cmap(1,:), 'lineStyle', '-.')
%an = circ_mean(cart2pol(bsRightB(:,2), -bsRightB(:,1)));
line([.245*cos(an) .25*cos(an)], [.245*sin(an) .25*sin(an)], 'lineWidth', 2, 'color', cmap(1,:))
plot(.245*cos(linspace(ll, ul, 50)), .245 * sin(linspace(ll, ul, 50)), 'color', cmap(1, :), 'lineWidth', 2);
%arrowhead(.19*[cos(an) sin(an)], .2*[cos(an) sin(an)], 'lineWidth', 2, 'color', cmap(1,:))


[y, t] = circ_kde(cart2pol(bsNEB(:,2), -bsNEB(:,1)), h);
y = [y y(1)]; t = [t t(1)];
plot(10*(cos(t) .* y), 10*(sin(t) .* y), 'lineWidth', 3, 'color', white)
plot(10*(cos(t) .* y), 10*(sin(t) .* y), 'lineWidth', 1.5, 'color', cmap(2,:))
[an ul ll] = circ_mean(cart2pol(bsNEB(:,2), -bsNEB(:,1)), bsNEr2(:));
%line(.2*[0 cos(an)], .2*[0 sin(an)], 'lineWidth', 2, 'color', cmap(2,:), 'lineStyle', '-.')
%an = circ_mean(cart2pol(bsNEB(:,2), -bsNEB(:,1)));
line([.235*cos(an) .24*cos(an)], [.235*sin(an) .24*sin(an)], 'lineWidth', 2, 'color', cmap(2,:))
plot(.235*cos(linspace(ll, ul, 50)), .235 * sin(linspace(ll, ul, 50)), 'color', cmap(2, :), 'lineWidth', 2);

%arrowhead([0 0], .2*[cos(an) sin(an)], 'lineWidth', 2, 'color', cmap(2,:))


[y, t] = circ_kde(cart2pol(bsTopB(:,2), -bsTopB(:,1)), h);
y = [y y(1)]; t = [t t(1)];
plot(10*(cos(t) .* y), 10*(sin(t) .* y), 'lineWidth', 3, 'color', white)
plot(10*(cos(t) .* y), 10*(sin(t) .* y), 'lineWidth', 1.5, 'color', cmap(3,:))
[an ul ll] = circ_mean(cart2pol(bsTopB(:,2), -bsTopB(:,1)), bsTopr2(:));
%line(.2*[0 cos(an)], .2*[0 sin(an)], 'lineWidth', 2, 'color', cmap(3,:), 'lineStyle', '-.')
%an = circ_mean(cart2pol(bsTopB(:,2), -bsTopB(:,1)));
line([.225*cos(an) .23*cos(an)], [.225*sin(an) .23*sin(an)], 'lineWidth', 2, 'color', cmap(3,:))
plot(.225*cos(linspace(ll, ul, 50)), .225 * sin(linspace(ll, ul, 50)), 'color', cmap(3, :), 'lineWidth', 2);
%arrowhead([0 0], .2*[cos(an) sin(an)], 'lineWidth', 2, 'color', cmap(3,:))

[y, t] = circ_kde(cart2pol(bsNWB(:,2), -bsNWB(:,1)), h);
y = [y y(1)]; t = [t t(1)];
plot(10*(cos(t) .* y), 10*(sin(t) .* y), 'lineWidth', 3, 'color', white)
plot(10*(cos(t) .* y), 10*(sin(t) .* y), 'lineWidth', 1.5, 'color', cmap(4,:))
[an ul ll] = circ_mean(cart2pol(bsNWB(:,2), -bsNWB(:,1)), bsNWr2(:));
%line(.2*[0 cos(an)], .2*[0 sin(an)], 'lineWidth', 2, 'color', cmap(4,:), 'lineStyle', '-.')
%an = circ_mean(cart2pol(bsNWB(:,2), -bsNWB(:,1)));
line([.215*cos(an) .22*cos(an)], [.215*sin(an) .22*sin(an)], 'lineWidth', 2, 'color', cmap(4,:))
plot(.215*cos(linspace(ll, ul, 50)), .215 * sin(linspace(ll, ul, 50)), 'color', cmap(4, :), 'lineWidth', 2);
%arrowhead([0 0], .2*[cos(an) sin(an)], 'lineWidth', 2, 'color', cmap(4,:))

[y, t] = circ_kde(cart2pol(bsLeftB(:,2), -bsLeftB(:,1)), h);
y = [y y(1)]; t = [t t(1)];
plot(10*(cos(t) .* y), 10*(sin(t) .* y), 'lineWidth', 3, 'color', white)
plot(10*(cos(t) .* y), 10*(sin(t) .* y), 'lineWidth', 1.5, 'color', cmap(5,:))
[an ul ll] = circ_mean(cart2pol(bsLeftB(:,2), -bsLeftB(:,1)), bsLeftr2(:));
%line(.2*[0 cos(an)], .2*[0 sin(an)], 'lineWidth', 2, 'color', cmap(5,:), 'lineStyle', '-.')
%an = circ_mean(cart2pol(bsLeftB(:,2), -bsLeftB(:,1)));
line([.205*cos(an) .21*cos(an)], [.205*sin(an) .21*sin(an)], 'lineWidth', 2, 'color', cmap(5,:))
plot(.205*cos(linspace(ll, ul, 50)), .205 * sin(linspace(ll, ul, 50)), 'color', cmap(5, :), 'lineWidth', 2);
%arrowhead([0 0], .2*[cos(an) sin(an)], 'lineWidth', 2, 'color', cmap(5,:))

[y, t] = circ_kde(cart2pol(bsSWB(:,2), -bsSWB(:,1)), h);
y = [y y(1)]; t = [t t(1)];
plot(10*(cos(t) .* y), 10*(sin(t) .* y), 'lineWidth', 3, 'color', white)
plot(10*(cos(t) .* y), 10*(sin(t) .* y), 'lineWidth', 1.5, 'color', cmap(6,:))
[an ul ll] = circ_mean(cart2pol(bsSWB(:,2), -bsSWB(:,1)), bsSWr2(:));
%line(.2*[0 cos(an)], .2*[0 sin(an)], 'lineWidth', 2, 'color', cmap(6,:), 'lineStyle', '-.')
%an = circ_mean(cart2pol(bsSWB(:,2), -bsSWB(:,1)));
line([.195*cos(an) .2*cos(an)], [.195*sin(an) .2*sin(an)], 'lineWidth', 2, 'color', cmap(6,:))
plot(.195*cos(linspace(ll, ul, 50)), .195 * sin(linspace(ll, ul, 50)), 'color', cmap(6, :), 'lineWidth', 2);
%arrowhead([0 0], .2*[cos(an) sin(an)], 'lineWidth', 2, 'color', cmap(6,:))

[y, t] = circ_kde(cart2pol(bsBottomB(:,2), -bsBottomB(:,1)), h);
y = [y y(1)]; t = [t t(1)];
plot(10*(cos(t) .* y), 10*(sin(t) .* y), 'lineWidth', 3, 'color', white)
plot(10*(cos(t) .* y), 10*(sin(t) .* y), 'lineWidth', 1.5, 'color', cmap(7,:))
[an ul ll] = circ_mean(cart2pol(bsBottomB(:,2), -bsBottomB(:,1)), bsBottomr2(:));
%line(.2*[0 cos(an)], .2*[0 sin(an)], 'lineWidth', 2, 'color', cmap(7,:), 'lineStyle', '-.')
%an = circ_mean(cart2pol(bsBottomB(:,2), -bsBottomB(:,1)));
line([.185*cos(an) .19*cos(an)], [.185*sin(an) .19*sin(an)], 'lineWidth', 2, 'color', cmap(7,:))
plot(.185*cos(linspace(ll, ul, 50)), .185 * sin(linspace(ll, ul, 50)), 'color', cmap(7, :), 'lineWidth', 2);
%arrowhead([0 0], .2*[cos(an) sin(an)], 'lineWidth', 2, 'color', cmap(7,:))

[y, t] = circ_kde(cart2pol(bsSEB(:,2), -bsSEB(:,1)), h);
y = [y y(1)]; t = [t t(1)];
plot(10*(cos(t) .* y), 10*(sin(t) .* y), 'lineWidth', 3, 'color', white)
plot(10*(cos(t) .* y), 10*(sin(t) .* y), 'lineWidth', 1.5, 'color', cmap(8,:))
[an ul ll] = circ_mean(cart2pol(bsSEB(:,2), -bsSEB(:,1)), bsSEr2(:));
%line(.2*[0 cos(an)], .2*[0 sin(an)], 'lineWidth', 2, 'color', cmap(8,:), 'lineStyle', '-.')
%an = circ_mean(cart2pol(bsSEB(:,2), -bsSEB(:,1)));
%an = circ_median(cart2pol(bsSEB(:,2), -bsSEB(:,1)));
line([.175*cos(an) .18*cos(an)], [.175*sin(an) .18*sin(an)], 'lineWidth', 2, 'color', cmap(8,:))
plot(.175*cos(linspace(ll, ul, 50)), .175 * sin(linspace(ll, ul, 50)), 'color', cmap(8, :), 'lineWidth', 2);
%arrowhead([0 0], .2*[cos(an) sin(an)], 'lineWidth', 2, 'color', cmap(8,:))
%}

axis equal; axis off;
%%

% Hypothesis testing
disp(['P(NE = SW) < ', num2str( ...
circ_kuipertest(mod(cart2pol(bsNEB(:,2), -bsNEB(:,1)),2*pi), ...
                mod(cart2pol(bsSWB(:,2), -bsSWB(:,1)), 2*pi), 100, 1), '%1.3f\n')]);
            
disp(['P(Top = Bottom) < ', num2str( ...
circ_kuipertest(mod(cart2pol(bsBottomB(:,2), -bsBottomB(:,1)),2*pi), ...
                mod(cart2pol(bsTopB(:,2), -bsTopB(:,1)), 2*pi), 100, 1), '%1.3f\n')]);
            
disp(['P(NW = SE) < ', num2str( ...
circ_kuipertest(mod(cart2pol(bsNWB(:,2), -bsNWB(:,1)),2*pi), ...
                mod(cart2pol(bsSEB(:,2), -bsSEB(:,1)), 2*pi), 100, 1), '%1.3f\n')]);
            
disp(['P(Right = Left) < ', num2str( ...
circ_kuipertest(mod(cart2pol(bsRightB(:,2), -bsRightB(:,1)),2*pi), ...
                mod(cart2pol(bsLeftB(:,2), -bsLeftB(:,1)), 2*pi), 100, 1), '%1.3f\n')]);
            
%%
bsB = [bsRightB; bsNEB; bsTopB; bsNWB; bsLeftB; bsSWB; bsBottomB; bsSEB];
% Hypothesis testing
disp(['P(unequal medians) < ', num2str( ...
    circ_cmtest(atan2(bsB(:,2), bsB(:,1)), ...
                reshape(ones(numBoot,1) * (1:8), [], 1)), '%1.8f\n')]);

%% 
% Mutual information between ABAT vector and target location

bsTopR = bsTopR ./ sum(bsTopR);
bsBottomR = bsBottomR ./ sum(bsBottomR);
bsLeftR = bsLeftR ./ sum(bsLeftR);
bsRightR = bsRightR ./ sum(bsRightR);
bsNWR = bsNWR ./ sum(bsNWR);
bsNER = bsNER ./ sum(bsNER);
bsSWR = bsSWR ./ sum(bsSWR);
bsSER = bsSER ./ sum(bsSER);

px = bsTopR + bsBottomR + bsLeftR + bsRightR + bsNWR + bsNER + bsSWR + bsSER;
px = px./sum(px);

% I(x;y) = sum_y p(y) [\sum_x p(x|y) log2( p(x|y) / p(x) )]

MI = .125 * nansum(bsTopR .* log2(bsTopR./(px+eps))) + ...
     .125 * nansum(bsBottomR .* log2(bsBottomR./(px+eps))) + ...
     .125 * nansum(bsLeftR .* log2(bsLeftR./(px+eps))) + ...
     .125 * nansum(bsRightR .* log2(bsRightR./(px+eps))) + ...
     .125 * nansum(bsNWR .* log2(bsNWR./(px+eps))) + ...
     .125 * nansum(bsNER .* log2(bsNER./(px+eps))) + ...
     .125 * nansum(bsSWR .* log2(bsSWR./(px+eps))) + ...
     .125 * nansum(bsSER .* log2(bsSER./(px+eps)));
%%
% shuffle test for significance 
d = cell2mat(arrayfun(@(x, y) repmat(x, y, 1), ...
    (edges(1:20) + diff(edges)/2)', round(px(:).*8000), 'un', 0));

numBoot = 1000;
bsMI = zeros(numBoot, 1);
for b = 1:numBoot
    z = reshape(d(randsample(1:numel(d), numel(d))), 8, 1000);
    
    Z = zeros(8, numel(edges));
    for i = 1:8
        Z(i,:) = histc(z(i,:), edges);
    end;
    Z = bsxfun(@rdivide, Z, sum(Z, 2));
    bsMI(b) = sum(.125 * nansum(Z(:, 1:end-1) .* log2(Z(:, 1:end-1) ./ repmat(px', 8, 1)), 2));
end;



%%

% Distribution of median ABAT times as a function of target location
bins = 870:10:1130;

subplot(3, 3, 1)
hist(bsNWT, bins);
set(get(gca, 'children'), 'facecolor', 'white')
xlim([min(bins), max(bins)])

subplot(3, 3, 2)
hist(bsTopT, bins);
set(get(gca, 'children'), 'facecolor', 'white')
xlim([min(bins), max(bins)])

subplot(3, 3, 3)
hist(bsNET, bins);
set(get(gca, 'children'), 'facecolor', 'white')
xlim([min(bins), max(bins)])

subplot(3, 3, 4)
hist(bsLeftT, bins);
set(get(gca, 'children'), 'facecolor', 'white')
xlim([min(bins), max(bins)])

subplot(3, 3, 6)
hist(bsRightT, bins);
set(get(gca, 'children'), 'facecolor', 'white')
xlim([min(bins), max(bins)])

subplot(3, 3, 7)
hist(bsSWT, bins);
set(get(gca, 'children'), 'facecolor', 'white')
xlim([min(bins), max(bins)])

subplot(3, 3, 8)
hist(bsBottomT, bins);
set(get(gca, 'children'), 'facecolor', 'white')
xlim([min(bins), max(bins)])

subplot(3, 3, 9)
hist(bsSET, bins);
set(get(gca, 'children'), 'facecolor', 'white')
xlim([min(bins), max(bins)])

%% 
% Mutual information between median ABAT time and target location

bsptopt = histc(bsTopT, bins) ./ numBoot;
bspbott = histc(bsBottomT, bins) ./ numBoot;
bspleftt = histc(bsLeftT, bins) ./ numBoot;
bsprightt = histc(bsRightT, bins) ./ numBoot;
bspnwt = histc(bsNWT, bins) ./ numBoot;
bspnet = histc(bsNET, bins) ./ numBoot;
bspswt = histc(bsSWT, bins) ./ numBoot;
bspset = histc(bsSET, bins) ./ numBoot;

px = bsptopt + bspbott + bspleftt + bsprightt + bspnwt + bspnet + bspswt + bspset;
px = px./sum(px);

% I(x;y) = sum_y p(y) [\sum_x p(x|y) log2( p(x|y) / p(x) )]

MI = .125 * nansum(bsptopt .* log2(bsptopt./(px+eps)) + ...
                   bspbott .* log2(bspbott./(px+eps)) + ...
                   bspleftt .* log2(bspleftt./(px+eps)) + ...
                   bsprightt .* log2(bsprightt./(px+eps)) + ...
                   bspnwt .* log2(bspnwt./(px+eps)) + ...
                   bspnet .* log2(bspnet./(px+eps)) + ...
                   bspswt .* log2(bspswt./(px+eps)) + ...
                   bspset .* log2(bspset./(px+eps)));

% shuffle test for significance 
d = cell2mat(arrayfun(@(x, y) repmat(x, y, 1), bins(:), round(px(:).*8000), 'un', 0));

numBoot = 1000;
bsMI = zeros(numBoot, 1);
for b = 1:numBoot
    z = reshape(d(randsample(1:numel(d), numel(d))), 8, 1000);
    
    Z = zeros(8, numel(bins));
    for i = 1:8
        Z(i,:) = histc(z(i,:), bins);
    end;
    Z = bsxfun(@rdivide, Z, sum(Z, 2));
    bsMI(b) = .125 * sum(nansum(Z .* log2(Z ./ repmat(px', 8, 1))));
end;
%%
%%
%%
%%
close all;
figure; 
%subplot(2, 1, 1)
set(gcf, 'position', [200 300 200 200])
hold on;

cmap = hsv(8);
h = .1;


TH = zeros(numBoot, 8);
TH(:,1) = cart2pol(bsRightB(:,2), -bsRightB(:,1));
TH(:,2) = cart2pol(bsNEB(:,2), -bsNEB(:,1));
TH(:,3) = cart2pol(bsTopB(:,2), -bsTopB(:,1));
TH(:,4) = cart2pol(bsNWB(:,2), -bsNWB(:,1));
TH(:,5) = cart2pol(bsLeftB(:,2), -bsLeftB(:,1));
TH(:,6) = cart2pol(bsSWB(:,2), -bsSWB(:,1));
TH(:,7) = cart2pol(bsBottomB(:,2), -bsBottomB(:,1));
TH(:,8) = cart2pol(bsSEB(:,2), -bsSEB(:,1));

R = [bsRightr2 bsNEr2 bsTopr2 bsNWr2 bsLeftr2 bsSWr2 bsBottomr2 bsSEr2];

[~, t] = circ_kde(cart2pol(bsRightB(:,2), -bsRightB(:,1)), h);
plot(.125*cos(t), .125*sin(t), 'color', 'k', 'lineWidth', 1);
plot(.25*cos(t), .25*sin(t), 'color', 'k', 'lineWidth', 1);

for i = 1:8
    [y, t] = circ_kde(TH(:,i), h);
    y = [y y(1)]; t = [t t(1)];
    plot(10*(cos(t) .* y), 10*(sin(t) .* y), 'lineWidth', 2.5, 'color', white)
    plot(10*(cos(t) .* y), 10*(sin(t) .* y), 'lineWidth', 1.25, 'color', lighter(cmap(i,:)))
    [an ul ll] = circ_mean(TH(:,i), R(:,i));
    
    line([(.245 - (.01 * (i-1))) * cos(an) (.25 - (.01 * (i-1))) * cos(an)], ...
         [(.245 - (.01 * (i-1))) * sin(an) (.25 - (.01 * (i-1))) * sin(an)], ...
         'lineWidth', 1.25, 'color', cmap(i,:))
    
     plot((.245 - (.01 * (i-1))) * cos(linspace(ll, ul, 50)), ...
          (.245 - (.01 * (i-1))) * sin(linspace(ll, ul, 50)), ...
          'color', cmap(i, :), 'lineWidth', 1.25);
end;

axis equal; axis off
%%
%close all;
figure;
%subplot(2, 1, 2)
set(gcf, 'position', [200 300 200 200])
hold on;
ann = circ_mean(TH(:), R(:));

for i = 1:8
    [an ul ll] = circ_mean(TH(:,i), R(:,i));
    
    plot([i i], ([ul ll] - ann) * 180/pi, 'color', cmap(i,:));
    plot(i, (an - ann) * 180/pi, 'k.')
    
end;
set(gca, 'xtick', 1:8)
xlim([0.5 8.5])
ylabel('relative angle')
xlabel('movement direction')

k = 22;
for i = 1:7
    for j = i+1:8
        [p, tab] = circ_wwtest([TH(:,i), TH(:,j)], [ones(1000,1); 2*ones(1000,1)], [R(:,i), R(:,j)]);
        if p < .05 / 28
            disp([num2str(i), ',', num2str(j)])
            line([i, j], [k k], 'color', 'k')
            k = k + 1;
        end;
    end;
end;

ylim([-20 25])

