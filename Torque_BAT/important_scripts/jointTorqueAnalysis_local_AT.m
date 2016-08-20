clear all;
load 'P-rs1050225_M1.mat'
load 'torqueRaw-rs1050225_M1.mat'

%load 'P-V1050917_M1Ab_PMvAc002.mat'
%load 'torqueRaw-V1050917_M1Ab_PMvAc002.mat'

%load 'P-V1050913_M1Ab_PMvAc00  1.mat'
%load 'torqueRaw-V1050913_M1Ab_PMvAc001.mat'


%%
pjt = nan(numTrials, 1);
tpjt = nan(numTrials, 1);

for i = 1:numTrials
    % relative to GO CUE
    idxs = find(kin.raw.stamps > beh(i,4) & kin.raw.stamps < beh(i,6));
    st = real(torque(idxs,1));
    et = real(torque(idxs,2));
    
    cjt = sqrt(st.^2 + et.^2);
    
    % Use 0.1 for RS1050225 and 0.03 for V1050917
    [pks, locs] = findpeaks(cjt, 'MinPeakHeight', .03);
    if ~isempty(locs)
        pjt(i) = cjt(locs(1));
        % relative to MO
        tpjt(i) = kin.raw.stamps(idxs(locs(1))) - beh(i,5);
    else
        disp(num2str(i))
    end;
end;
%%
close all;
scatter(1:numTrials, tpjt, [], beh(:,8), 'fill')
%%
close all;
colormap(hsv(8));
scatter(beh(:,5) - beh(:,4), tpjt, [], beh(:,8), 'fill')


%%
%close all;
figure;
hold on;
plot(beh(:,8), tpjt, '.')
plot(grpstats(tpjt, beh(:,8), 'median'), 'k')
xlim([.5 8.5])
%ylim([0 .1])

%%
close all;

plot((nanmean([abatRight(MIchans), abatNE(MIchans), abatTop(MIchans), abatNW(MIchans), ...
     abatLeft(MIchans), abatSW(MIchans), abatBottom(MIchans), abatSE(MIchans)])-1000)./1000, ...
     grpstats(tpjt, beh(:,8), 'median'), '.','MarkerSize',24)
%xlim([.5 8.5])

%%

torqueProfiles = zeros(numel(find(kin.raw.stamps > beh(i,5) - 0.5 & kin.raw.stamps < beh(i,5) + .8)), numTrials);
for i = 1:numTrials
    idxs = find(kin.raw.stamps > beh(i,5) - 0.5 & kin.raw.stamps < beh(i,5) + .8);
    
    
    st = real(torque(idxs,1));
    et = real(torque(idxs,2));
    
    cjt = sqrt(st.^2 + et.^2);
    
    torqueProfiles(:,i) = cjt;
end;

%%
close all;
figure;
hold on;
plot(kin.raw.stamps(idxs) - beh(end,5), torqueProfiles(:,281), 'r', 'lineWidth', 2)
plot(kin.raw.stamps(idxs) - beh(end,5), torqueProfiles(:,135), 'color', [1 .4 .6], 'lineWidth', 2)
plot(kin.raw.stamps(idxs) - beh(end,5), torqueProfiles(:,255), 'color', [1 .5 .5], 'lineWidth', 2)

plot(kin.raw.stamps(idxs) - beh(end,5), torqueProfiles(:,45), 'b', 'lineWidth', 2)
plot(kin.raw.stamps(idxs) - beh(end,5), torqueProfiles(:,105), 'color', [0 .5 1], 'lineWidth', 2)
plot(kin.raw.stamps(idxs) - beh(end,5), torqueProfiles(:,119), 'color', [.5 .25 1], 'lineWidth', 2)

xlim([-.25, 0.5])
ylim([0 .28])
xlabel('time (s)')
ylabel('joint torque magnitude (au)')
%%
close all;
figure
for i = 1:numTrials
    if beh(i,8) == 2
        plot(kin.raw.stamps(idxs) - beh(end,5), torqueProfiles(:,i), 'r')
        title(num2str(i))
        pause;
    end;
end;

%%
close all;

plot(1:8, grpstats(tpjt, beh(:,8), 'mean') + 2 * grpstats(tpjt, beh(:,8), 'mean'))

%%
%%
%%
%%
%%Adil: 7/26 These are mean BAT by direction
%nanmean(([abatRight(MIchans) abatNE(MIchans) abatTop(MIchans) abatNW(MIchans) abatLeft(MIchans) abatSW(MIchans) abatBottom(MIchans) abatSE(MIchans)]-1000)./1000)
%the above script is taken from running RSABATbyDirScript
rsba = [-0.0141    0.0057   -0.0647   -0.0423   -0.0295    0.0336   -0.0056   -0.0112];
%These are from grpstats(tpjtrs, beh(:,8), 'mean')
rsjt = [ 0.0300    0.0264    0.0232    0.0325    0.0447    0.0578    0.0263    0.0263];
%These are time of torque onset: grpstats(timeOfTorqueOnsetjtrs, beh(:,8), 'mean')
rsto = [ -0.1432 -0.1616 -0.1762 -0.1680 -0.1570 -0.1653 -0.1599 -0.1403];
%method 2 of computing torque onset:
rsto = [ -0.0871 -0.1005 -0.1205 -0.0894 -0.0757 -0.0841 -0.0964 -0.0903];

vba17 = [-0.0720   -0.1062   -0.1180   -0.0758   -0.0470   -0.0516   -0.0418   -0.0623];
vjt17 = [0.0345    0.0147   -0.0200    0.0145    0.0370    0.0403    0.0686    0.0514];

vba13 = [-0.0828   -0.0850   -0.0683   -0.0057   -0.0375    0.0096   -0.0602   -0.0682];
vjt13 = [0.0340    0.0133    0.0057    0.0445    0.0511    0.0590    0.0468    0.0600];






%rsba = grpstats(xRS, behRS(:,8), 'mean');
%vba = grpstats(xV, behV(:,8), 'mean');


%vba17 = [-0.1840   -0.2113   -0.1557   -0.1416   -0.1465   -0.1867   -0.2067   -0.1660];

close all;

figure;
% subplot(2,1,1)
% hold on;
% plot(rsba, rsjt, 'r.','MarkerSize',24)
% plot(vba17, vjt17, 'b.','MarkerSize',24)
% plot(vba13, vjt13, '.','MarkerSize',24, 'color', [0 .5 1])
% subplot(2,1,2)
hold on
plot(rsba, rsto, 'b.','MarkerSize',24)
plot(rsba, rsjt, 'r.','MarkerSize',24)
legend('onset','peak')
xlabel('Beta Attenuation Time')
ylabel('Torque time; onset and peak')
title('Beta attenuation vs Torque, binned by movement direction')

linScale = -.08:0.01:.04;
pp = polyfit(rsba, rsjt,1);
po = polyfit(rsba, rsto,1);
yfitp = polyval(pp,rsba);
yfito = polyval(po,rsba);
yresidp = rsjt - yfitp;
yresido = rsto - yfito;

SSresidp = sum(yresidp.^2);
SStotalp = sum((rsjt - mean(rsjt)).^2);

SSresido = sum(yresido.^2);
SStotalo = sum((rsto - mean(rsto)).^2);
%SStotal = (length(timeOfTorqueOnsetjtrs_nansremoved)-1 * var(timeOfTorqueOnsetjtrs_nansremoved));
rsqp = 1 - (SSresidp/SStotalp);
rsqo = 1 - (SSresido/SStotalo);
plot(linScale, pp(1) * linScale + pp(2),'r');
plot(linScale, po(1) * linScale + po(2),'b');
text(max(xlim)*.1,max(ylim)*.6,strcat('Rsq Peak Torque = ',num2str(rsqp)),'FontSize',14)
text(max(xlim)*.1,max(ylim)*.4,strcat('Rsq Torque Onset = ',num2str(rsqo)),'FontSize',14)



%%
close all;
aoctool([rsba(:); vba13(:); vba17(:)], [rsjt(:); vjt13(:); vjt17(:)], [ones(8,1); 2*ones(8,1); 2*ones(8,1)])

%%
numBoot = 10000;
B = zeros(numBoot, 1);
for i = 1:numBoot
    B(i) = corr(rsws(:), reshape(rsba(randsample(1:8, 8)), [], 1));
end;

%%
ds = dataset([rsba(:); vba(:)], [rsjt(:); vjt(:)], nominal([ones(8,1)-1; 1*ones(8,1)]), 'varNames', {'ba', 'jt', 'm'});

fullmdl = LinearModel.fit(ds, 'jt ~ m * ba');
nestedmdl = LinearModel.fit(ds, 'jt ~ m + m*ba - ba');

%%
close all;
figure; hold on;
plot(behRS(:,8) - .1, tpjtRS, '.', 'color', [.6 .6 .7])
plot((1:8)-0, grpstats(tpjtRS, behRS(:,8), 'median'))

plot(behV(:,8) + .1, tpjtV, '.', 'color', [.7 .6 .6])
plot((1:8)+0, grpstats(tpjtV, behV(:,8), 'median'), 'r')