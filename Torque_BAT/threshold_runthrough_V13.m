clear all;

load 'P-V1050913_M1Ab_PMvAc001.mat'
load 'torqueRaw-V1050913_M1Ab_PMvAc001.mat'


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
thresholdRange = 0.001:0.01:0.999;
linearModelStatistics = struct;
allR2 = nan(length(thresholdRange),1);
allSlope = nan(length(thresholdRange),1);
allIntercept = nan(length(thresholdRange),1);


t = 1;
for THRESHOLD = .06
    THRESHOLD
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
            tO = getTorqueOnset(cjt,2,THRESHOLD,i,'v');
            torqueOnsetMagjtv13(i) = cjt(tO);
            timeOfTorqueOnsetjtv13(i) = kin.raw.stamps(idxs(tO))-behv13(i,5);
%             saveas(gcf,strcat('example_',num2str(i)),'epsc')

%             pause
        else
            disp(num2str(i))
        end;
    end;
    
    vba13 = [-0.0828   -0.0850   -0.0683   -0.0057   -0.0375    0.0096   -0.0602   -0.0682];
    vjt13 = [0.0340    0.0133    0.0057    0.0445    0.0511    0.0590    0.0468    0.0600];

    vto13 = grpstats(timeOfTorqueOnsetjtv13, behv13(:,8), 'mean')';

    close all
    figure
    hold on

    plot(vba13, vto13, 'b.','MarkerSize',24)
    plot(vba13, vjt13, 'r.','MarkerSize',24)
    legend('onset','peak')
    xlabel('Beta Attenuation Time')
    ylabel('Torque time; onset and peak')
    title('Beta attenuation vs Torque, binned by movement direction (V13)')

    linScale = -.1:0.01:.02;
    pp = polyfit(vba13, vjt13,1);
    po = polyfit(vba13, vto13,1);
    yfitp = polyval(pp,vba13);
    yfito = polyval(po,vba13);
    yresidp = vjt13 - yfitp;
    yresido = vto13 - yfito;

    SSresidp = sum(yresidp.^2);
    SStotalp = sum((vjt13- mean(vjt13)).^2);

    SSresido = sum(yresido.^2);
    SStotalo = sum((vto13 - mean(vto13)).^2);
    rsqp = 1 - (SSresidp/SStotalp);
    rsqo = 1 - (SSresido/SStotalo);
    plot(linScale, pp(1) * linScale + pp(2),'r');
    plot(linScale, po(1) * linScale + po(2),'b');
    ylim([-.2 .1])
    text(max(xlim)*.1,max(ylim)*.6,strcat('Rsq Peak Torque = ',num2str(rsqp)),'FontSize',14)
    text(max(xlim)*.1,max(ylim)*.4,strcat('Rsq T. Onset = ',num2str(rsqo)),'FontSize',14)
    text(max(xlim)*.1,max(ylim)*.2,strcat('Slope T. Onset = ',num2str(po(1))),'FontSize',14)
    text(max(xlim)*.1,max(ylim)*.01,strcat('Intercept T. Onset = ',num2str(po(2))),'FontSize',14)
%     saveas(gcf,strcat('threshold_',num2str(THRESHOLD*1000)),'epsc');
    [THRESHOLD,rsqp,rsqo]
    allR2(t) = rsqo;
    allSlope(t) = po(1);
    allIntercept(t) = po(2);
    t = t + 1;
    close all
        
    [B, BINT, ~, ~, STATS] = regress(vto13', [ones(size(vba13')) vba13'])
    disp('      R2, F, P-VALUE, EST. VARIANCE ERROR')
    
end

linearModelStatistics.threshold = thresholdRange;
linearModelStatistics.slope = allSlope;
linearModelStatistics.intercept = allIntercept;
linearModelStatistics.R2 = allR2;

plot(allR2)
    
