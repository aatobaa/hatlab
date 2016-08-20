%%
%close all;
figure;
hold on;
cmap = hsv(8);

for d = 1:8
    plot(((1:2001)-1000)./1000, mean(mean(abs(H(:,beh(:,8) == d,MIchans)), 3), 2), 'color', cmap(d,:), 'lineWidth', 2);
    %plot(mean(mean(abs(H(:,beh(:,8) == d,:)), 3), 2) + sem(mean(abs(H(:,beh(:,8) == d,:)), 3), 2), 'color', cmap(d,:));
    %plot(mean(mean(abs(H(:,beh(:,8) == d,:)), 3), 2) - sem(mean(abs(H(:,beh(:,8) == d,:)), 3), 2), 'color', cmap(d,:));
    errorArea(((1:2001)-1000)./1000,mean(mean(abs(H(:,beh(:,8) == d,MIchans)), 3), 2), ...
        2*sem(mean(abs(H(:,beh(:,8) == d,:)), 3), 2), ...
        2*sem(mean(abs(H(:,beh(:,8) == d,:)), 3), 2), 'edgeColor', cmap(d,:), 'faceColor', cmap(d,:))
end;

xlim([-.5 .5])

%%
%%
close all;
% eidx = 35 52 66 are good exemplars for RS
for eidx = 1:96
figure;
hold on;
cmap = hsv(8);


for d = 1:8
    plot(((1:2001)-1000)./1000, mean(abs(H(:, beh(:,8) == d, eidx)), 2), 'color', cmap(d,:), 'lineWidth', 2);
    %plot(mean(mean(abs(H(:,beh(:,8) == d,:)), 3), 2) + sem(mean(abs(H(:,beh(:,8) == d,:)), 3), 2), 'color', cmap(d,:));
    %plot(mean(mean(abs(H(:,beh(:,8) == d,:)), 3), 2) - sem(mean(abs(H(:,beh(:,8) == d,:)), 3), 2), 'color', cmap(d,:));
    errorArea(((1:2001)-1000)./1000, mean(abs(H(:, beh(:,8) == d, eidx)), 2), ...
        2*sem(abs(H(:,:,eidx)), 2), ...
        2*sem(abs(H(:,:,eidx)), 2), 'edgeColor', cmap(d,:), 'faceColor', cmap(d,:))
end;

xlim([-.5 .5])
end;
%%
%%
%close all;
figure;
%set(gcf, 'position', [100 100 300 round(5/3 * 300)]);
set(gcf, 'position', [100 100 500 200])
subplot(1, 2, 2);

hold on; 
plot(beh(:,8), (beh(:,5) - beh(:,4)) * 1000, '.', 'color', [.6 .6 .7])
plot(1:8, grpstats((beh(:,5) - beh(:,4)) * 1000, beh(:,8), 'median'), 'color', 'k')
xlim([.5 8.5])
set(gca, 'xtick', 1:8);
ylabel('reaction time (ms)')

subplot(1, 2, 1);
hold on;
plot(reshape(ones(numel(MIchans), 1) * (1:8), [], 1), ...
    ([abatRight(MIchans); abatNE(MIchans); abatTop(MIchans); abatNW(MIchans); ...
     abatLeft(MIchans); abatSW(MIchans); abatBottom(MIchans); abatSE(MIchans)]-1000)./1000, ...
     '.', 'color', [.6 .6 .7])
 
plot(1:8, (median([abatRight(MIchans), abatNE(MIchans), abatTop(MIchans), abatNW(MIchans), ...
     abatLeft(MIchans), abatSW(MIchans), abatBottom(MIchans), abatSE(MIchans)])-1000)./1000, 'lineWidth', 1, 'color', [0 0 0])
xlim([.5 8.5])

%
x = reshape(ones(numel(MIchans), 1) * (1:8), [], 1);
y = ([abatRight(MIchans); abatNE(MIchans); abatTop(MIchans); abatNW(MIchans); ...
     abatLeft(MIchans); abatSW(MIchans); abatBottom(MIchans); abatSE(MIchans)]-1000)./1000;
 
myFun = @(x, z) z(1)*cos((x-1)*2*pi/8 + z(2)) + z(3);
[z, r] = lsqcurvefit(@(z, x) myFun(x, z), ones(3,1), x, y);

plot(1:.01:8, myFun(1:.01:8, z), 'lineWidth', 1, 'color', rgbStr2Mat('6b7ce0'))

disp(['R2 = ', num2str(1 - (r ./ sum((y- mean(y)).^2)))]);
%}
set(gca, 'xtick', 1:8);
xlabel('movement direction');
ylabel('beta attenuation time');