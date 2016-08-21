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



thresholdRange = 0.001:0.01:0.999;
linearModelStatistics = struct;
allR2 = nan(length(thresholdRange),1);
allSlope = nan(length(thresholdRange),1);
allIntercept = nan(length(thresholdRange),1);
noisyTrials = nan(numTrials,1);

t = 1;
for THRESHOLD = thresholdRange
    THRESHOLD
    
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
            [tO, SKIP] = getTorqueOnset(cjt,2,THRESHOLD,i,'v');
            noisyTrials(i) = SKIP; 
            torqueOnsetMagjtv17(i) = cjt(tO);
            timeOfTorqueOnsetjtv17(i) = kin.raw.stamps(idxs(tO))-behv17(i,5);
            SKIP;
%             pause
        else
            disp(num2str(i))
        end;
    end;
    
    %These are trial-aligned to movement onset...
    vba17 = [-0.0720   -0.1062   -0.1180   -0.0758   -0.0470   -0.0516   -0.0418   -0.0623];
    vjt17 = [0.0345    0.0147   -0.0200    0.0145    0.0370    0.0403    0.0686    0.0514];

    vto17 = grpstats(timeOfTorqueOnsetjtv17, behv17(:,8), 'mean')';
%     vto17 = grpstats(timeOfTorqueOnsetjtv17(~noisyTrials), behv17(~noisyTrials,8), 'mean')';

    close all
    figure
    hold on

    plot(vba17, vto17, 'b.','MarkerSize',24)
    plot(vba17, vjt17, 'r.','MarkerSize',24)
    legend('onset','peak')
    xlabel('Beta Attenuation Time')
    ylabel('Torque time; onset and peak')
    title('Beta attenuation vs Torque, binned by movement direction (V17)')

    linScale = -.1:0.01:.02;
    pp = polyfit(vba17, vjt17,1);
    po = polyfit(vba17, vto17,1);
    yfitp = polyval(pp,vba17);
    yfito = polyval(po,vba17);
    yresidp = vjt17 - yfitp;
    yresido = vto17 - yfito;

    SSresidp = sum(yresidp.^2);
    SStotalp = sum((vjt17- mean(vjt17)).^2);

    SSresido = sum(yresido.^2);
    SStotalo = sum((vto17 - mean(vto17)).^2);
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
    
end

linearModelStatistics.threshold = thresholdRange;
linearModelStatistics.slope = allSlope;
linearModelStatistics.intercept = allIntercept;
linearModelStatistics.R2 = allR2;

plot(allR2)
    
