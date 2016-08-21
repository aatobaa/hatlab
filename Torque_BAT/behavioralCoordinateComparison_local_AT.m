% Run animalStrMIABAT2 first, then load p-file and torque file if
% necessary.
wristSpeed = sqrt( kin.raw.xvel.^2 + kin.raw.yvel.^2);
numTrials = size(beh, 1);
wristSpeedProfiles = zeros(500, numTrials);
xVelProfiles = zeros(500, numTrials);
yVelProfiles = zeros(500, numTrials);

for i = 1:numTrials
    idxs = find(kin.raw.stamps > beh(i,5) - .5 & kin.raw.stamps < beh(i,5) + .5);
    wristSpeedProfiles(:,i) = wristSpeed(idxs);
    xVelProfiles(:,i) = kin.raw.xvel(idxs);
    yVelProfiles(:,i) = kin.raw.yvel(idxs);
end;


% Find x and y velocities at the moment of maximum wrist speed
[~, idx] = max(wristSpeedProfiles);

maxXVel = zeros(numTrials, 1);
maxYVel = zeros(numTrials, 1);

for i = 1:numTrials
    maxXVel(i) = xVelProfiles(idx(i), i);
    maxYVel(i) = yVelProfiles(idx(i), i);
end;



%%
%close all;

figure;
hold on;
r = mean(sqrt(maxXVel.^2 + maxYVel.^2));

plot(r * cos(0:.01:2*pi), r * sin(0:.01:2*pi), 'color', .5 * [1 1 1])

colormap(hsv(10))
scatter(maxXVel, maxYVel, 100 * ones(size(beh(:,8))), beh(:,8), '.')
axis square;

axis(.75 * [-1 1 -1 1])
xlabel('x velocity')
ylabel('y velocity')

%%
%close all;
% Arm Trajectory
figure
hold on;
cmap = hsv(8);
for i = 1:numTrials
    idxs = find(kin.raw.stamps > beh(i,5) & kin.raw.stamps < beh(i,6));
    plot(kin.raw.x(idxs), kin.raw.y(idxs), 'color', cmap(beh(i,8),:));
end;
xlabel('trajectories in xy');
axis equal;

%%
% mean cartesian velocity
figure;
mwv = zeros(numTrials, 2);
% THIS IS FOR VELMA
%stamps = kin.raw.stamps(1):.05:kin.raw.stamps(end); stamps = stamps(2:end) - kin.raw.stamps(1);

% THIS IS FOR ROCKSTAR
stamps = 0:.05:kin.raw.stamps(end); stamps = stamps(2:end);
for i = 1:numTrials
    idxs = find(stamps > beh(i,5) & stamps < beh(i,6));
    mwv(i,:) = [mean(kin.binned.xvel(idxs)), mean(kin.binned.yvel(idxs))];
end;

colormap(hsv(10));
scatter(mwv(:,1), mwv(:,2), [], beh(:,8), '.')
xlabel('mean x velocity')
ylabel('mean y velocity')
axis equal;

%%
% mean joint torque
figure;
mjt = zeros(numTrials, 2);
for i = 1:numTrials
    idxs = find(stamps > beh(i,5) & stamps < beh(i,6));
        mjt(i,:) = mean(torque(idxs,:));
        %disp(num2str(i));
end;

colormap(hsv(10));
scatter(mjt(:,1), mjt(:,2), [], beh(:,8), '.')
xlabel('mean shoulder torque')
ylabel('mean elbow torque');
axis equal;

%%
% mean joint velocity
k1=kinfo.upperarm;
k2=kinfo.lowerarm;

x = kin.binned.x;
y = kin.binned.y;

r = sqrt(x.^2 + y.^2);
elb = acos( (r.^2 - k1^2 - k2^2)./(2*k1*k2));

switch kinfo.arm
    case 'left'
        sho = pi - atan2(y,x) - acos( (r.^2 + k1^2 - k2^2)./(2*k1*r));
    case 'right'
        sho = atan2(y,x) - acos( (r.^2 + k1^2 - k2^2)./(2*k1*r));
end

elbv = diff(elb)/kinfo.binSize;
elbv(end+1)=elbv(end);
shov = diff(sho)/kinfo.binSize;
shov(end+1) = shov(end);

mse = zeros(numTrials, 2);

for i = 1:numTrials
    idxs = find(stamps > beh(i,5) & stamps < beh(i,6));
    mse(i,:) = [mean(shov(idxs)), mean(elbv(idxs))];
end;

figure; 
colormap(hsv(10));
scatter(mse(:,1), mse(:,2), [], beh(:,8), '.');
xlabel('mean shoulder velocity')
ylabel('mean elbow velocity')
axis equal;
%%
% Arm trajectory in shoulder elbow space
figure; hold on;
for i = 1:numTrials
    idxs = find(stamps > beh(i,5) & stamps < beh(i,6));
    plot(sho(idxs), elb(idxs), 'color', cmap(beh(i,8),:));
end;
xlabel('shoulder')
ylabel('elbow')
axis equal;

%%
% mean joint power
mjp = zeros(numTrials, 2);

for i = 1:numTrials
    idxs = find(stamps > beh(i,5) & stamps < beh(i,5) + .5);
    mjp(i,:) = [mean(torque(idxs,1)' * shov(idxs)), ...
                mean(torque(idxs,2)' * elbv(idxs))];
end;

figure;
colormap(hsv(10));
scatter(mjp(:,1), mjp(:,2), [], beh(:,8), '.')

%% 
% peak cartesian velocity
figure;
pwv = zeros(numTrials, 2);
% THIS IS FOR VELMA
%stamps = kin.raw.stamps(1):.05:kin.raw.stamps(end); stamps = stamps(2:end) - kin.raw.stamps(1);

% THIS IS FOR ROCKSTAR
stamps = 0:.05:kin.raw.stamps(end); stamps = stamps(2:end);
for i = 1:numTrials
    idxs = find(stamps > beh(i,5) & stamps < beh(i,6));
    xv = kin.binned.xvel(idxs);
    yv = kin.binned.yvel(idxs);
    idx = find(sqrt(xv.^2 + yv.^2) == max(sqrt(xv.^2 + yv.^2)), 1);
    pwv(i,:) = [xv(idx), yv(idx)];
end;

colormap(hsv(10));
scatter(pwv(:,1), pwv(:,2), [], beh(:,8), '.')
xlabel('x velocity at peak wrist speed')
ylabel('y velocity at peak wrist speed')
axis equal;
%%
% peak joint torque
figure;
pjt = zeros(numTrials, 2);
tpjt = zeros(numTrials, 1);
pjtt = zeros(numTrials, 1);
for i = 1:numTrials
    idxs = find(stamps > beh(i,5) & stamps < beh(i,6) + 0);
    st = torque(idxs,1);
    et = torque(idxs,2);
    idx = find(sqrt(st.^2 + et.^2) == max(sqrt(st.^2 + et.^2)), 1);
    pjt(i,:) = [st(idx), et(idx)];
    %pjt(i,:) = [max(st(abs(st) == max(abs(st)))), max(et(abs(et) == max(abs(et))))];
    cjt = sqrt(st.^2 + et.^2);
    tpjt(i) = find(abs(cjt) == max(abs(cjt)));
    pjtt(i) = max(cjt(abs(cjt) == max(abs(cjt))));
end;

colormap(hsv(10));
scatter(pjt(:,1), pjt(:,2), [], beh(:,8), '.')
xlabel('peak shoulder torque')
ylabel('peak elbow torque');
axis equal;

%%
% peak joint velocity
figure;
pjv = zeros(numTrials, 2);
for i = 1:numTrials
    idxs = find(stamps > beh(i,5) & stamps < beh(i,6));
    sv = shov(idxs);
    ev = elbv(idxs);
    %idx = find(sqrt(xv.^2 + yv.^2) == max(sqrt(xv.^2 + yv.^2)), 1);
    %pwv(i,:) = [xv(idx), yv(idx)];
    pjv(i,:) = [max(sv(abs(sv) == max(abs(sv)))), max(ev(abs(ev) == max(abs(ev))))];
end;

colormap(hsv(10));
scatter(pjv(:,1), pjv(:,2), [], beh(:,8), '.');
xlabel('peak shoulder velocity');
ylabel('peak elbow velocity');
axis equal;
%%
% peak joint power
pjp = zeros(numTrials, 1);
for i = 1:numTrials
    idxs = find(stamps > beh(i,5) & stamps < beh(i,6));
    pjp(i) = max(sqrt( (shov(idxs) .* torque(idxs,1) + elbv(idxs) .* torque(idxs, 2)).^2));
    %sjp = 
    %pjp(i,:) = 
end
%%
% time of peak cartesian velocity

figure;
tpwv = zeros(numTrials, 2);

for i = 1:numTrials
    idxs = find(kin.raw.stamps > beh(i,5) & kin.raw.stamps < beh(i,6));
    tpwv(i,:) = [find(abs(kin.raw.xvel(idxs)) == max(abs(kin.raw.xvel(idxs)))), find(abs(kin.raw.yvel(idxs)) == max(abs(kin.raw.yvel(idxs))))];
end;

colormap(hsv(10));
scatter(tpwv(:,1), tpwv(:,2), [], beh(:,8), '.');

%%
k1=kinfo.upperarm;
k2=kinfo.lowerarm;

x = kin.raw.x - mean(kin.raw.x) + mean(kin.binned.x);
y = kin.raw.y - mean(kin.raw.y) + mean(kin.binned.y);

r = sqrt(x.^2 + y.^2);
kin.raw.elb = acos( (r.^2 - k1^2 - k2^2)./(2*k1*k2));

switch kinfo.arm
    case 'left'
        kin.raw.sho = pi - atan2(y,x) - acos( (r.^2 + k1^2 - k2^2)./(2*k1*r));
    case 'right'
        kin.raw.sho = atan2(y,x) - acos( (r.^2 + k1^2 - k2^2)./(2*k1*r));
end

kin.raw.elbv = diff(smoothWithGaussianKernel(kin.raw.elb, .05, 500))/.002;
kin.raw.elbv(end+1)=kin.raw.elbv(end);
kin.raw.shov = diff(smoothWithGaussianKernel(kin.raw.sho, .05, 500))/.002;
kin.raw.shov(end+1) = kin.raw.shov(end);
%%
% time of peak joint velocity
figure;
tpjv = zeros(numTrials, 2);
for i = 1:numTrials
    idxs = find(kin.raw.stamps > beh(i,5) & kin.raw.stamps < beh(i,6));
    sv = kin.raw.shov(idxs);
    svb = shov(stamps > beh(i,5) & stamps < beh(i,6));
    ev = kin.raw.elbv(idxs);
    evb = elbv(stamps > beh(i,5) & stamps < beh(i,6));
    %idx = find(sqrt(xv.^2 + yv.^2) == max(sqrt(xv.^2 + yv.^2)), 1);
    %pwv(i,:) = [xv(idx), yv(idx)];
    tpjv(i,:) = [find(abs(sv) == max(abs(sv))), find(abs(ev) == max(abs(ev)))];
    
    %{
    subplot(2, 1, 1)
    cla;
    hold on;
    plot(sv);
    %plot(stamps(stamps > beh(i,5) & stamps < beh(i,6)), svb, 'r');
    line(tpjv(i,1) * [1 1], get(gca, 'ylim'))
    
    subplot(2, 1, 2)
    cla;
    hold on;
    plot(ev);
    %plot(evb, 'r');
    line(tpjv(i,2) * [1 1], get(gca, 'ylim'));
    
    title(num2str(i));
    pause;
    %}
end;

colormap(hsv(10));
scatter(tpjv(:,1), tpjv(:,2), [], beh(:,8), '.');
xlabel('peak shoulder velocity');
ylabel('peak elbow velocity');
axis equal;

%%

figure;
tpjt = zeros(numTrials, 2);
for i = 1:numTrials
    idxs = find(kin.raw.stamps > beh(i,5) & kin.raw.stamps < beh(i,6));
    tpjt(i,:) = [find(abs(torque(idxs,1)) == max(abs(torque(idxs,1)))), ...
                 find(abs(torque(idxs,2)) == max(abs(torque(idxs,2))))];
end;

colormap(hsv(10));
scatter(tpjt(:,1), tpjt(:,2), [], beh(:,8), '.');
xlabel('time of peak shoulder torque')
ylabel('time of peak elbow torque');

%%
% time of peak joint torque
figure;
tpjt = zeros(numTrials, 1);

for i = 1:numTrials
end;

%% 
% Steve scott 2001 nature paper figure
close all;

% THIS IS FOR VELMA
%stamps = kin.raw.stamps(1):.05:kin.raw.stamps(end); stamps = stamps(2:end) - kin.raw.stamps(1);

% THIS IS FOR ROCKSTAR
stamps = 0:.05:kin.raw.stamps(end); stamps = stamps(2:end);
figure;
set(gcf, 'position', [3 500 1270 430]);
subplot(1, 4, 1)
% PEAK WRIST SPEED
pws = zeros(numTrials, 1);
for i = 1:numTrials
    idxs = find(stamps > beh(i,5) & stamps < beh(i,6));
    xv = kin.binned.xvel(idxs);
    yv = kin.binned.yvel(idxs);
    pws(i) = max(sqrt(xv.^2 + yv.^2));
end;

hold on;
%plot(pws .* cos((beh(:,8)-1) * pi/4), pws .* sin((beh(:,8)-1) * pi/4), 'k.');
for i = 0:7
    line(cos(i * pi/4) * [quantile(pws(beh(:,8) == i+1), .25) quantile(pws(beh(:,8) == i + 1), .75)], ...
         sin(i * pi/4) * [quantile(pws(beh(:,8) == i+1), .25) quantile(pws(beh(:,8) == i + 1), .75)], ...
         'lineWidth', 2, 'color', [.5 .5 .5]);
    plot(median(pws(beh(:,8) == i+1)) * cos(i * pi/4), median(pws(beh(:,8) == i+1)) * sin(i * pi/4), 'k.', 'markerSize', 20);
end;
axis equal;
axis off;

subplot(1, 4, 2);
% PEAK SHOULDER ELBOW VELOCITY
pse = zeros(numTrials, 1);
for i = 1:numTrials
    idxs = find(stamps > beh(i,5) & stamps < beh(i,6));
    pse(i) = max(sqrt(shov(idxs).^2 + elbv(idxs).^2));
end;

hold on;
%plot(pse .* cos((beh(:,8)-1) * pi/4), pse .* sin((beh(:,8)-1) * pi/4), 'k.');
for i = 0:7
    line(cos(i * pi/4) * [quantile(pse(beh(:,8) == i+1), .25) quantile(pse(beh(:,8) == i + 1), .75)], ...
         sin(i * pi/4) * [quantile(pse(beh(:,8) == i+1), .25) quantile(pse(beh(:,8) == i + 1), .75)], ...
         'lineWidth', 2, 'color', [.5 .5 .5]);
    plot(median(pse(beh(:,8) == i+1)) * cos(i * pi/4), median(pse(beh(:,8) == i+1)) * sin(i * pi/4), 'k.', 'markerSize', 20);
end;
axis equal;
axis off;

subplot(1, 4, 3);
% PEAK JOINT TORQUE
pjt = zeros(numTrials, 1);
for i = 1:numTrials
    idxs = find(stamps > beh(i,5) & stamps < beh(i,6));
    pjt(i) = max(sqrt(sum(torque(idxs,:).^2, 2)));
end;

hold on;
%plot(pjt .* cos((beh(:,8)-1) * pi/4), pjt .* sin((beh(:,8)-1) * pi/4), 'k.');
for i = 0:7
    line(cos(i * pi/4) * [quantile(pjt(beh(:,8) == i+1), .25) quantile(pjt(beh(:,8) == i + 1), .75)], ...
         sin(i * pi/4) * [quantile(pjt(beh(:,8) == i+1), .25) quantile(pjt(beh(:,8) == i + 1), .75)], ...
         'lineWidth', 2, 'color', [.5 .5 .5]);
    plot(median(pjt(beh(:,8) == i+1)) * cos(i * pi/4), median(pjt(beh(:,8) == i+1)) * sin(i * pi/4), 'k.', 'markerSize', 20);
end;
axis equal;
axis off;

subplot(1, 4, 4);
pjp = zeros(numTrials, 1);
for i = 1:numTrials
    idxs = find(stamps > beh(i,5) & stamps < beh(i,6));
    pjp(i) = max(sqrt( (shov(idxs) .* torque(idxs,1) + elbv(idxs) .* torque(idxs, 2)).^2 ...
        ));
end

hold on;
%plot(pjp .* cos((beh(:,8) - 1) * pi/4), pjp .* sin((beh(:,8)-1) * pi/4), 'k.');
for i = 0:7
    line(cos(i * pi/4) * [quantile(pjp(beh(:,8) == i+1), .25) quantile(pjp(beh(:,8) == i + 1), .75)], ...
         sin(i * pi/4) * [quantile(pjp(beh(:,8) == i+1), .25) quantile(pjp(beh(:,8) == i + 1), .75)], ...
         'lineWidth', 2, 'color', [.5 .5 .5]);
    plot(median(pjp(beh(:,8) == i+1)) * cos(i * pi/4), median(pjp(beh(:,8) == i+1)) * sin(i * pi/4), 'k.', 'markerSize', 20);
end;
axis equal;
axis off;
%%
close all;
sffm = events(15).times(ia);
effm = events(16).times(ia);
for i = 1:numTrials
    %idxs = find(stamps > beh(i,5) & stamps < beh(i,6));
    %idxs = find(stamps > sffm(i) & stamps < effm(i));
    idxs = find(stamps > min(beh(i,5), sffm(i)) & stamps < max(beh(i,6), effm(i)));
    
    subplot(5, 1, 1)
    plot(sqrt(kin.binned.xvel(idxs,1).^2 + kin.binned.yvel(idxs,1).^2))
    xlabel('wrist speed')
    
    subplot(5, 1, 2)
    plot(kin.binned.xvel(idxs,1));
    xlabel('xvel')
    
    subplot(5, 1, 3)
    plot(kin.binned.yvel(idxs,1));
    xlabel('yvel')
    
    subplot(5, 1, 4)
    plot(torque(idxs,1))
    xlabel('shoulder torque')
    
    subplot(5, 1, 5)
    plot(torque(idxs,2))
    xlabel('elbow torque')
    
    title(num2str(i));
    pause;
end;

%%
%%
% Money shot!
close all;
cmap = hsv(8);
figure;
set(gcf, 'position', [360 100 240 690])
%set(gcf, 'position', [360 130 230 800])
subplot(3, 1, 1)
% y and x come from RS1050225ABATbyDirAmplitudeTest
hold on;
disp('wrist speed')
mdl = LinearModel.fit(grpstats(y, x, 'mean'), grpstats(sqrt(sum(pwv.^2, 2)), beh(:,8), 'mean'))
% horizontal error bars
line(grpstats(y, x, 'meanci')', repmat(grpstats(sqrt(sum(pwv.^2, 2)), beh(:,8), 'mean'), 1, 2)', ...
    'color', [.5 .5 .5], 'lineWidth', 2)
line(repmat(grpstats(y, x, 'mean'), 1, 2)', grpstats(sqrt(sum(pwv.^2, 2)), beh(:,8), 'meanci')', ...
    'color', [.5 .5 .5], 'lineWidth', 2)
plot(grpstats(y, x, 'mean'), grpstats(sqrt(sum(pwv.^2, 2)), beh(:,8), 'mean'), ...
    'k.', 'markerSize', 13)
colormap(hsv(10))
%scatter(grpstats(y, x, 'mean'), grpstats(sqrt(sum(pwv.^2, 2)), beh(:,8), 'mean'), 2000 * ones(8,1), 1:8, '.')
xlim([-.1 .05])
lsline
xlabel('bat')
ylabel('peak wrist speed')
title(num2str(mdl.Rsquared.Ordinary, '%1.2f'));



subplot(3, 1, 2)
hold on;
disp('joint speed')
mdl = LinearModel.fit(grpstats(y, x, 'mean'), grpstats(sqrt(sum(pjv.^2, 2)), beh(:,8), 'mean'))
line(grpstats(y, x, 'meanci')', repmat(grpstats(sqrt(sum(pjv.^2, 2)), beh(:,8), 'mean'), 1, 2)', ...
    'color', [.5 .5 .5], 'lineWidth', 2)
line(repmat(grpstats(y, x, 'mean'), 1, 2)', grpstats(sqrt(sum(pjv.^2, 2)), beh(:,8), 'meanci')', ...
    'color', [.5 .5 .5], 'lineWidth', 2)
plot(grpstats(y, x, 'mean'), grpstats(sqrt(sum(pjv.^2, 2)), beh(:,8), 'mean'), ...
    'k.', 'markersize', 13)


colormap(hsv(10))
%scatter(grpstats(y, x, 'mean'), grpstats(sqrt(sum(pjv.^2, 2)), beh(:,8), 'mean'), 2000 * ones(8,1), 1:8, '.')
%xlim([-.1 .05])
lsline
xlabel('bat')
ylabel('peak joint speed')
title(num2str(mdl.Rsquared.Ordinary, '%1.2f'));

subplot(3, 1, 3)
hold on;
disp('joint torque')
mdl = LinearModel.fit(grpstats(y, x, 'mean'), grpstats(sqrt(sum(pjt.^2, 2)), beh(:,8), 'mean'))
line(grpstats(y, x, 'meanci')', repmat(grpstats(sqrt(sum(pjt.^2, 2)), beh(:,8), 'mean'), 1, 2)', ...
    'color', [.5 .5 .5], 'lineWidth', 2)
line(repmat(grpstats(y, x, 'mean'), 1, 2)', grpstats(sqrt(sum(pjt.^2, 2)), beh(:,8), 'meanci')', ...
    'color', [.5 .5 .5], 'lineWidth', 2)
plot(grpstats(y, x, 'mean'), grpstats(sqrt(sum(pjt.^2, 2)), beh(:,8), 'mean'), ...
    'k.', 'markerSize', 13)


colormap(hsv(10))
%scatter(grpstats(y, x, 'mean'), grpstats(sqrt(sum(pjt.^2, 2)), beh(:,8), 'mean'), 2000 * ones(8,1), 1:8, '.')
%xlim([-.1 .05])
lsline
xlabel('bat')
ylabel('peak joint torque')
title(num2str(mdl.Rsquared.Ordinary, '%1.2f'));

%%
% Supplemental figure using medians instead of means
close all;
figure;
set(gcf, 'position', [360 100 240 690])
%set(gcf, 'position', [360 130 230 800])
subplot(3, 1, 1)
% y and x come from RS1050225ABATbyDirAmplitudeTest
hold on;
disp('wrist speed')
mdl = LinearModel.fit(grpstats(y, x, 'median'), grpstats(sqrt(sum(pwv.^2, 2)), beh(:,8), 'median'))
% horizontal error bars
%line(grpstats(y, x, 'meanci')', repmat(grpstats(sqrt(sum(pwv.^2, 2)), beh(:,8), 'mean'), 1, 2)', ...
 %   'color', [.5 .5 .5], 'lineWidth', 2)
%line(repmat(grpstats(y, x, 'mean'), 1, 2)', grpstats(sqrt(sum(pwv.^2, 2)), beh(:,8), 'meanci')', ...
 %   'color', [.5 .5 .5], 'lineWidth', 2)
plot(grpstats(y, x, 'median'), grpstats(sqrt(sum(pwv.^2, 2)), beh(:,8), 'median'), ...
    'k.', 'markerSize', 13)
xlim([-.1 .05])
lsline
xlabel('bat')
ylabel('peak wrist speed')
title(num2str(mdl.Rsquared.Ordinary, '%1.2f'));



subplot(3, 1, 2)
hold on;
disp('joint speed')
mdl = LinearModel.fit(grpstats(y, x, 'median'), grpstats(sqrt(sum(pjv.^2, 2)), beh(:,8), 'median'))
%line(grpstats(y, x, 'meanci')', repmat(grpstats(sqrt(sum(pjv.^2, 2)), beh(:,8), 'mean'), 1, 2)', ...
 %   'color', [.5 .5 .5], 'lineWidth', 2)
%line(repmat(grpstats(y, x, 'mean'), 1, 2)', grpstats(sqrt(sum(pjv.^2, 2)), beh(:,8), 'meanci')', ...
 %   'color', [.5 .5 .5], 'lineWidth', 2)
plot(grpstats(y, x, 'median'), grpstats(sqrt(sum(pjv.^2, 2)), beh(:,8), 'median'), ...
    'k.', 'markersize', 13)
xlim([-.1 .05])
lsline
xlabel('bat')
ylabel('peak joint speed')
title(num2str(mdl.Rsquared.Ordinary, '%1.2f'));

subplot(3, 1, 3)
hold on;
disp('joint torque')
mdl = LinearModel.fit(grpstats(y, x, 'median'), grpstats(sqrt(sum(pjt.^2, 2)), beh(:,8), 'median'))
%line(grpstats(y, x, 'meanci')', repmat(grpstats(sqrt(sum(pjt.^2, 2)), beh(:,8), 'mean'), 1, 2)', ...
 %   'color', [.5 .5 .5], 'lineWidth', 2)
%line(repmat(grpstats(y, x, 'mean'), 1, 2)', grpstats(sqrt(sum(pjt.^2, 2)), beh(:,8), 'meanci')', ...
 %   'color', [.5 .5 .5], 'lineWidth', 2)
plot(grpstats(y, x, 'median'), grpstats(sqrt(sum(pjt.^2, 2)), beh(:,8), 'median'), ...
    'k.', 'markerSize', 13)
xlim([-.1 .05])
lsline
xlabel('bat')
ylabel('peak joint torque')
title(num2str(mdl.Rsquared.Ordinary, '%1.2f'));

%%
% Money shot with joint power
close all;
figure;
%set(gcf, 'position', [360 130 230 800])
subplot(2, 2, 1)
% y and x come from RS1050225ABATbyDirAmplitudeTest
hold on;
mdl = LinearModel.fit(grpstats(y, x, 'mean'), grpstats(sqrt(sum(pwv.^2, 2)), beh(:,8), 'mean'));
% horizontal error bars
line(grpstats(y, x, 'meanci')', repmat(grpstats(sqrt(sum(pwv.^2, 2)), beh(:,8), 'mean'), 1, 2)', ...
    'color', [.5 .5 .5], 'lineWidth', 2)
line(repmat(grpstats(y, x, 'mean'), 1, 2)', grpstats(sqrt(sum(pwv.^2, 2)), beh(:,8), 'meanci')', ...
    'color', [.5 .5 .5], 'lineWidth', 2)
plot(grpstats(y, x, 'mean'), grpstats(sqrt(sum(pwv.^2, 2)), beh(:,8), 'mean'), ...
    'k.', 'markerSize', 13)
lsline
xlabel('bat')
ylabel('peak wrist speed')
title(num2str(mdl.Rsquared.Ordinary, '%1.2f'));

subplot(2, 2, 2)
hold on;
mdl = LinearModel.fit(grpstats(y, x, 'mean'), grpstats(sqrt(sum(pjv.^2, 2)), beh(:,8), 'mean'));
line(grpstats(y, x, 'meanci')', repmat(grpstats(sqrt(sum(pjv.^2, 2)), beh(:,8), 'mean'), 1, 2)', ...
    'color', [.5 .5 .5], 'lineWidth', 2)
line(repmat(grpstats(y, x, 'mean'), 1, 2)', grpstats(sqrt(sum(pjv.^2, 2)), beh(:,8), 'meanci')', ...
    'color', [.5 .5 .5], 'lineWidth', 2)
plot(grpstats(y, x, 'mean'), grpstats(sqrt(sum(pjv.^2, 2)), beh(:,8), 'mean'), ...
    'k.', 'markersize', 13)
lsline
xlabel('bat')
ylabel('peak joint speed')
title(num2str(mdl.Rsquared.Ordinary, '%1.2f'));

subplot(2, 2, 3)
hold on;
mdl = LinearModel.fit(grpstats(y, x, 'mean'), grpstats(sqrt(sum(pjt.^2, 2)), beh(:,8), 'mean'));
line(grpstats(y, x, 'meanci')', repmat(grpstats(sqrt(sum(pjt.^2, 2)), beh(:,8), 'mean'), 1, 2)', ...
    'color', [.5 .5 .5], 'lineWidth', 2)
line(repmat(grpstats(y, x, 'mean'), 1, 2)', grpstats(sqrt(sum(pjt.^2, 2)), beh(:,8), 'meanci')', ...
    'color', [.5 .5 .5], 'lineWidth', 2)
plot(grpstats(y, x, 'mean'), grpstats(sqrt(sum(pjt.^2, 2)), beh(:,8), 'mean'), ...
    'k.', 'markerSize', 13)
lsline
xlabel('bat')
ylabel('peak joint torque')
title(num2str(mdl.Rsquared.Ordinary, '%1.2f'));

subplot(2, 2, 4)
hold on;
mdl = LinearModel.fit(grpstats(y, x, 'mean'), grpstats(sqrt(sum(pjp.^2, 2)), beh(:,8), 'mean'));
line(grpstats(y, x, 'meanci')', repmat(grpstats(sqrt(sum(pjp.^2, 2)), beh(:,8), 'mean'), 1, 2)', ...
    'color', [.5 .5 .5], 'lineWidth', 2)
line(repmat(grpstats(y, x, 'mean'), 1, 2)', grpstats(sqrt(sum(pjp.^2, 2)), beh(:,8), 'meanci')', ...
    'color', [.5 .5 .5], 'lineWidth', 2)
plot(grpstats(y, x, 'mean'), grpstats(sqrt(sum(pjp.^2, 2)), beh(:,8), 'mean'), ...
    'k.', 'markerSize', 13)
lsline
xlabel('bat')
ylabel('peak joint power')
title(num2str(mdl.Rsquared.Ordinary, '%1.2f'));

%%
% Use these data to train a multinomial logistic regression model to
% predict movement direction from kinematic quantities
close all;
kinQuant = mjt;


figure; hold on;

colormap(hsv(10));
scatter(kinQuant(:,1), kinQuant(:,2), [], beh(:,8), '.')



k = convhull(kinQuant(:,1), kinQuant(:,2));
[mgX, mgY] = meshgrid(linspace(min(kinQuant(:,1)), max(kinQuant(:,1)), 70), linspace(min(kinQuant(:,2)), max(kinQuant(:,2)), 70));
inPoly = inpolygon(mgX(:), mgY(:), kinQuant(k, 1), kinQuant(k, 2));

mgX = mgX(inPoly);
mgY = mgY(inPoly);

plot(kinQuant(k,1), kinQuant(k,2), 'k');



B = mnrfit(kinQuant, beh(:,8), 'model', 'nominal', 'interactions', 'on');
YHATp = mnrval(B, [mgX(:) mgY(:)], 'model', 'nominal', 'interactions', 'on');

[~, YHAT] = max(YHATp, [], 2);

colormap(hsv(10));
scatter(mgX(:), mgY(:), [], YHAT, 'x')
axis equal;
%%
figure;
BATIMES = YHATp * modelFitBATbyDir';
%probs = YHATp * (1:8)';


colormap( flipud([colorSpace(rgbStr2Mat('b2182b'), rgbStr2Mat('d6604d'), 32); ...
           colorSpace(rgbStr2Mat('d6604d'), rgbStr2Mat('f4a582'), 32); ...
           colorSpace(rgbStr2Mat('f4a582'), rgbStr2Mat('fddbc7'), 32); ...
           colorSpace(rgbStr2Mat('fddbc7'), rgbStr2Mat('f7f7f7'), 32); ...
           colorSpace(rgbStr2Mat('f7f7f7'), rgbStr2Mat('d1e5f0'), 32); ...
           colorSpace(rgbStr2Mat('d1e5f0'), rgbStr2Mat('92c5de'), 32); ...
           colorSpace(rgbStr2Mat('92c5de'), rgbStr2Mat('4393c3'), 32); ...
           colorSpace(rgbStr2Mat('4393c3'), rgbStr2Mat('2166ac'), 32)]));
scatter(mgX(:), mgY(:), [], BATIMES, 's', 'fill');
%surf(mgX, mgY, reshape(probs, 100, 100), 'edgeAlpha', 0);
%view([0 90]);
%axis([min(kinQuant(:,1)), max(kinQuant(:,1)), min(kinQuant(:,2)), max(kinQuant(:,2)) ]);



%imagesc(linspace(min(kinQuant(:,1)), max(kinQuant(:,1)), 100), ...
%        linspace(min(kinQuant(:,2)), max(kinQuant(:,2)), 100), ...
%        reshape(probs, 100, 100));

%label('mean x velocity')
%ylabel('mean y velocity')


% find the point with max entropy, and then draw the line segments
entr = -sum(YHATp .*log2(YHATp), 2);
line([mgX(BATIMES == min(BATIMES), 1), mgX(entr == max(entr), 1)], ...
     [mgY(BATIMES == min(BATIMES), 1), mgY(entr == max(entr), 1)], ...
     'lineWidth', 2, 'color', 'k');
 
line([mgX(BATIMES == max(BATIMES), 1), mgX(entr == max(entr), 1)], ...
     [mgY(BATIMES == max(BATIMES), 1), mgY(entr == max(entr), 1)], ...
     'lineWidth', 2, 'color', 'k');
 
axis equal;


%%
figure;
colormap(hot)
%imagesc(linspace(min(kinQuant(:,1)), max(kinQuant(:,1)), 100), ...
%        linspace(min(kinQuant(:,2)), max(kinQuant(:,2)), 100), ...
%        reshape(-sum(YHATp .* log2(YHATp), 2), 100, 100));
scatter(mgX(:), mgY(:), [], -sum(YHATp .* log2(YHATp), 2), 's', 'fill')
axis equal;

%%
close all;
figure;

for i = 1:numTrials
    idxs = find(stamps > beh(i,5) & stamps < beh(i,6));
    st = torque(idxs,1);
    et = torque(idxs,2);
    %idx = find(sqrt(st.^2 + et.^2) == max(sqrt(st.^2 + et.^2)), 1);
    %pjt(i,:) = [st(idx), et(idx)];
    %pjt(i,:) = [max(st(abs(st) == max(abs(st)))), max(et(abs(et) == max(abs(et))))];
    if beh(i,8) == 1
        subplot(3, 3, 6)
        hold on;
        plot(((1:length(idxs))-1) * .05, sqrt(st.^2 + et.^2));
        ylim([0 .5])
    elseif beh(i,8) == 2
        subplot(3, 3, 3)
        hold on;
        plot(((1:length(idxs))-1) * .05, sqrt(st.^2 + et.^2));
        ylim([0 .5])
    elseif beh(i,8) == 3
        subplot(3, 3, 2)
        hold on;
        plot(((1:length(idxs))-1) * .05, sqrt(st.^2 + et.^2));
        ylim([0 .5])
    elseif beh(i,8) == 4
        subplot(3, 3, 1)
        hold on;
        plot(((1:length(idxs))-1) * .05, sqrt(st.^2 + et.^2));
        ylim([0 .5])
    elseif beh(i,8) == 5
        subplot(3, 3, 4)
        hold on;
        plot(((1:length(idxs))-1) * .05, sqrt(st.^2 + et.^2));
        ylim([0 .5])
        ylabel('torque magnitude')
    elseif beh(i,8) == 6
        subplot(3, 3, 7)
        hold on;
        plot(((1:length(idxs))-1) * .05, sqrt(st.^2 + et.^2));
        ylim([0 .5])
    elseif beh(i,8) == 7
        subplot(3, 3, 8)
        hold on;
        plot(((1:length(idxs))-1) * .05, sqrt(st.^2 + et.^2));
        ylim([0 .5])
        xlabel('time (s)')
    elseif beh(i,8) == 8
        subplot(3, 3, 9)
        hold on;
        plot(((1:length(idxs))-1) * .05, sqrt(st.^2 + et.^2));
        ylim([0 .5])
    end;
end;

%%
torqueProfiles = nan(numTrials, 17);

for i = 1:numTrials
    idxs = find(stamps > beh(i,5) & stamps < beh(i,6));
    st = torque(idxs,1);
    et = torque(idxs,2);
    
    torqueProfiles(i,1:numel(idxs)) = sqrt(st.^2 + et.^2);
end;

%%
close all;
figure;

for i = 1:8
    if i == 1
        subplot(3, 3, 6)
        hold on;
        plot(0:.05:.8, nanmedian(torqueProfiles(beh(:,8) == 1, :)))
        ylim([0 .3])
        xlim([0, .5])
    elseif i == 2
        subplot(3, 3, 3)
        hold on;
        plot(0:.05:.8, nanmedian(torqueProfiles(beh(:,8) == 2, :)))
        ylim([0 .3])
        xlim([0, .5])
    elseif i == 3
        subplot(3, 3, 2)
        hold on;
        plot(0:.05:.8, nanmedian(torqueProfiles(beh(:,8) == 3, :)))
        ylim([0 .3])
        xlim([0, .5])
    elseif i == 4
        subplot(3, 3, 1)
        hold on;
        plot(0:.05:.8, nanmedian(torqueProfiles(beh(:,8) == 4, :)))
        ylim([0 .3])
        xlim([0, .5])
    elseif i == 5
        subplot(3, 3, 4)
        hold on;
        plot(0:.05:.8, nanmedian(torqueProfiles(beh(:,8) == 5, :)))
        ylim([0 .3])
        xlim([0, .5])
        ylabel('torque magnitude')
    elseif i == 6
        subplot(3, 3, 7)
        hold on;
        plot(0:.05:.8, nanmedian(torqueProfiles(beh(:,8) == 6, :)))
        ylim([0 .3])
        xlim([0, .5])
    elseif i == 7
        subplot(3, 3, 8)
        hold on;
        plot(0:.05:.8, nanmedian(torqueProfiles(beh(:,8) == 7, :)))
        ylim([0 .3])
        xlim([0, .5])
        xlabel('time (s)')
    elseif i == 8
        subplot(3, 3, 9)
        hold on;
        plot(0:.05:.8, nanmedian(torqueProfiles(beh(:,8) == 8, :)))
        ylim([0 .3])
        xlim([0, .5])
    end;
end;






