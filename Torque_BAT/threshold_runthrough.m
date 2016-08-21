% clear all;

load 'P-rs1050225_M1.mat'
load 'torqueRaw-rs1050225_M1.mat'
load 'rs1050225_MI_clean_LFP.mat'

% add the following to your path
% addpath(genpath('``server``:\Matt'))



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

REMOVETRIALS = [1,6,42,44,55,60,63,81,90,104,124,128,130,135,143,151,163,166,179,182,186,195,204,213,223,224,228];
GOODTRIALS = ~ismember(1:numTrials,REMOVETRIALS);

thresholdRange = 0.001:0.3:0.999;
linearModelStatistics = struct;
allR2 = nan(length(thresholdRange),1);
allSlope = nan(length(thresholdRange),1);
allIntercept = nan(length(thresholdRange),1);

t = 1;
for THRESHOLD = .06
    THRESHOLD
    pjtrs = nan(numTrials, 1);
    %time of peak joint torque; relative to go cue
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
            % relative to Go Cue
            tpjtrs(i) = kin.raw.stamps(idxs(locs(1))) - behrs(i,4);

            %Adil 7/26: torque Onset instead of peak torque
            tO = getTorqueOnset(cjt,2,THRESHOLD,i,'rs');
            torqueOnsetMagjtrs(i) = cjt(tO);
            timeOfTorqueOnsetjtrs(i) = kin.raw.stamps(idxs(tO))-behrs(i,4);

        else
            disp(num2str(i))
        end;
    %     subplot(1,3,1), plot(st)
    %     ylabel('shoulder')
    %     subplot(1,3,2), plot(et)
    %     ylabel('elbow')
    %     subplot(1,3,3), plot(cjt)
    %     ylabel('combined')
%         saveas(gcf,strcat('example_',num2str(i)),'epsc')
%         i
%         pause
    end;


    clear lfp*



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    rsba = [-0.0141    0.0057   -0.0647   -0.0423   -0.0295    0.0336   -0.0056   -0.0112];
    %These are from grpstats(tpjtrs, beh(:,8), 'mean')
    rsjt = [ 0.0300    0.0264    0.0232    0.0325    0.0447    0.0578    0.0263    0.0263];
    %These are time of torque onset: grpstats(timeOfTorqueOnsetjtrs, beh(:,8), 'mean')
    %rsto = [ -0.1432 -0.1616 -0.1762 -0.1680 -0.1570 -0.1653 -0.1599 -0.1403];
    %method 2 of computing torque onset:
    %rsto = [ -0.0871 -0.1005 -0.1205 -0.0894 -0.0757 -0.0841 -0.0964 -0.0903];
    rsto = grpstats(timeOfTorqueOnsetjtrs, behrs(:,8), 'mean')';
    %rsto = grpstats(timeOfTorqueOnsetjtrs(GOODTRIALS), behrs(GOODTRIALS,8), 'mean')';
    

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
    title('Beta attenuation vs Torque, binned by movement direction (RS)')

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
%%
linearModelStatistics.threshold = thresholdRange;
linearModelStatistics.slope = allSlope;
linearModelStatistics.intercept = allIntercept;
linearModelStatistics.R2 = allR2;

plot(allR2)
