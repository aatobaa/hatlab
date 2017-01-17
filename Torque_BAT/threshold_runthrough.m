% clear all;

% If you get "Subscripted assignment dimension mismatch error, then this is
% probably due to storing cjt into cjt_torque_profiles without using
% constant dimensions of 1001 elements (+- 1 second) 

% Each trial takes roughly 10 seconds; beh is measured in seconds. The
% first movement onset occurs roughly 40 seconds after the session starts
% (i.e. at the first value of kin.raw.stamps). Kinematics are sampled at
% 500 Hz.

load 'P-rs1050225_M1.mat'
load 'torqueRaw-rs1050225_M1.mat'
load 'rs1050225_MI_clean_LFP.mat'

% add the following to your path
% addpath(genpath('``server``:\Matt'))


% "We excluded trials with reactions times less than 100 ms and greater
% than 750 ms. Movement durations had to be between 150 ms and 1500 ms.  
reactionTimeUB = .75;
reactionTimeLB = .1;
movementTimeUB = 1.5;
movementTimeLB = .15;

behrs = beh;
behrs = removeSlowReactionTimesFromBEH(behrs, reactionTimeUB, reactionTimeLB);
behrs = removeSlowMovementTimesFromBEH(behrs, movementTimeUB, movementTimeLB);

% "We excluded the slowest 25% of trials from each movement direction"
md = zeros(8,1);
for i = 1:8
    md(i) = quantile(behrs(behrs(:,8) == i, 6) - behrs(behrs(:,8) == i, 5), .75);
end;

fastTrials = (behrs(:,6) - behrs(:,5)) < md(behrs(:,8));
behrs = behrs(fastTrials,:);

numTrials = size(behrs, 1);
eventTimes = round(behrs(:, 5) * 1000);

thresholdRange = 0.001:0.01:0.999; % Threshold for determining torque onset
linearModelStatistics = struct;
allR2 = nan(length(thresholdRange),1);
allSlope = nan(length(thresholdRange),1);
allIntercept = nan(length(thresholdRange),1);

cjt_torque_profiles = zeros(1000,numTrials);
st_torque_profiles = zeros(1000,numTrials);
et_torque_profiles = zeros(1000,numTrials); 

peak_torque_times_st = zeros(numTrials,1);
peak_torque_times_et = zeros(numTrials,1);
peak_torque_times_cjt = zeros(numTrials,1);
onset_torque_times_st = zeros(numTrials,1);
onset_torque_times_et = zeros(numTrials,1);
onset_torque_times_cjt = zeros(numTrials,1);


t = 1;
for THRESHOLD = 0.06
    THRESHOLD
    pjtrs = nan(numTrials, 1);
    %time of peak joint torque; relative to go cue
    tpjtrs = nan(numTrials, 1);
    torqueOnsetjtrs = nan(numTrials,1);
    timeOfTorqueOnsetjtrs = nan(numTrials,1);

    for i = 1:numTrials
        
        %% Select torque indices +- 1 second relative to movement onset.
        % Used for Crosscorrelation test
        AlignPlusMinusOne = 1;
        if AlignPlusMinusOne
            idxs = find(kin.raw.stamps > behrs(i,5) - 1 & kin.raw.stamps < behrs(i,5) + 1);
        end
        
        %% Select torque indices from Go to End of Movement.
        % Used for all other analyses
        AlignToGo = 0;
        if AlignToGo
        % relative to GO CUE
            idxs = find(kin.raw.stamps > behrs(i,4) & kin.raw.stamps < behrs(i,6));
        end
        
        %% Select Torque Profile of interest
        % Selects combined joint torque. Other joints unimplemented. Can
        % explore the effect of normalizing elbow and shoulder before
        % combining, and also look at each joint individually. 
        st = real(torque(idxs,1)); % shoulder torque
        et = real(torque(idxs,2)); % elbow torque
        cjt = sqrt(st.^2 + et.^2); % combined torque
        
        normalizeSubScript = 0;
        if normalizeSubScript
        %% Normalize shoulder and elbow EXPERIMENT:
        % 9/4/16 Test: Defining combined torque as the NORMALIZED combination of
        % shoulder and elbow torque, so modulations in either joint affect the
        % combined equally, regardless of the relative magnitude of the individual
        % componenent. 
        % 
            % Normalize profiles
            min_st = min(st);
            max_st= max(st);
            min_et = min(et);
            max_et = max(et);

            a=-1; b=1;
            % Using formula from http://www.mathworks.com/matlabcentral/fileexchange/5103-toolbox-diffc/content/toolbox_diffc/toolbox/rescale.m
            norm_st = (b-a) .* (st - min_st)./(max_st-min_st) + a;
            norm_et = (b-a) .* (et - min_et)./(max_et-min_et) + a; 

            st = norm_st;
            et = norm_et;
            cjt = sqrt(st.^2 + et.^2);

            doPlot = 0;
            if doPlot                
                subplot(4,1,1)
                getTorqueOnset(sqrt(st.^2),500,2,THRESHOLD,i,'norm_rs');
                ylabel('shoulder')

    %             legend('show','torque','torque onset','movement onset','peak torque');

                subplot(4,1,2)
                getTorqueOnset(sqrt(et.^2),500,2,THRESHOLD,i,'norm_rs');
                ylabel('elbow')

                subplot(4,1,3)
                tO = getTorqueOnset(cjt(1:end),500,2,THRESHOLD,i,'norm_rs');
                ylabel('combined')

                subplot(4,1,4)
                plot(1:2:2000,beta_profiles(1:1000,i),'LineWidth',2)
                ylabel('beta amplitude')
            end
            tO = getTorqueOnset(cjt(1:end),500,2,THRESHOLD,i,'norm_rs');
        end
        
        moSubPlots = 0;
        if moSubPlots
        %% Show torque profiles split by shoulder, elbow, combined, + beta
            if i == 233
                continue
            end
            subplot(4,1,1)
            [sto, spt, ~] = getTorqueOnset(sqrt(st.^2),500,2,THRESHOLD,i,'v');
            ylabel('shoulder')

            legend('show','torque','torque onset','movement onset','peak torque');

            subplot(4,1,2)
            [eto, ept, ~] = getTorqueOnset(sqrt(et.^2),500,2,THRESHOLD,i,'v');
            ylabel('elbow')

            subplot(4,1,3)
            [tO, cpt, ~] = getTorqueOnset(cjt(1:end),500,2,THRESHOLD,i,'v');
            ylabel('combined')
            subplot(4,1,4)
            plot(1:2:2000,beta_profiles(1:1000,i),'LineWidth',2)
            ylabel('beta amplitude') 

            peak_torque_times_st(i) = spt;
            peak_torque_times_et(i) = ept;
            peak_torque_times_cjt(i) = cpt;
            onset_torque_times_st(i) = sto;
            onset_torque_times_et(i) = eto;
            onset_torque_times_cjt(i) = tO;
            
            doPlot = 0;
            if doPlot
            %% This section addresses: "Movement Onset is not always 
             % identified at peak torque / torque onset
                subplot(2,1,1)
                hold on
                plot(peak_torque_times_st-500,'.','Color',[0 0 1],'MarkerSize',6)
                %B(1) is the intercept
                %B(2) is the slope
                [B, BINT_pst, ~, ~, STATS_pst] = regress(peak_torque_times_st' - 500,[ones(281,1) (1:281)']);
    %                 linScale = 1:281;
    %                 plot(linScale, B(1),'Color',[0 0 1],'LineWidth',13);

                plot(peak_torque_times_et-500,'.','Color', [0 0.9 0.2],'MarkerSize',6)
                [B, BINT_pet, ~, ~, STATS_pet] = regress(peak_torque_times_et' - 500,[ones(281,1) (1:281)']);
    %                 linScale = 1:281;
    %                 plot(linScale, B(2) * linScale + B(1),'Color',[0 0.9 0.2],'LineWidth',3);

                plot(peak_torque_times_cjt-500,'.', 'Color', [1 .8 0],'MarkerSize',6)
                [B, BINT_pcjt, ~, ~, STATS_pcjt] = regress(peak_torque_times_cjt' - 500,[ones(281,1) (1:281)']);
    %                 linScale = 1:281;
    %                 plot(linScale, B(2) * linScale + B(1),'Color',[1 .8 0],'LineWidth',3);


                line([0 281],[0 0],'LineStyle','--','Color','r','LineWidth',1)
                xlim([0 281])
                ylabel('Peak Torque')
                legend('show','shoulder torque','elbow torque','combined torque')
                title('Torque Events relative to Movement Onset')
                subplot(2,1,2)
                hold
                plot(onset_torque_times_st-500,'.','Color',[0 0 1],'MarkerSize',10)
                [B, BINT_ost, ~, ~, STATS_ost] = regress(onset_torque_times_st' - 500,[ones(281,1) (1:281)']);

                plot(onset_torque_times_et-500,'.','Color', [0 0.9 0.2],'MarkerSize',10)
                [B, BINT_oet, ~, ~, STATS_oet] = regress(onset_torque_times_et' - 500,[ones(281,1) (1:281)']);

                plot(onset_torque_times_cjt-500,'.','Color', [1 .8 0],'MarkerSize',10)
                [B, BINT_ocjt, ~, ~, STATS_ocjt] = regress(onset_torque_times_cjt' - 500,[ones(281,1) (1:281)']);

                line([0 281],[0 0],'LineStyle','--','Color','r','LineWidth',1)
                xlim([0 281])
                ylabel('Torque Onset')

                xlabel('Trial Number')
            end
        end
        
        %% Find Peak Torque magnitude and time 
        % Use 0.1 for RS1050225 and 0.03 for V1050917
        % If you get an error here, make sure the right "findpeaks"
        % function is being called. Should be in the Matlab Signal Toolbox
        [pks, locs] = findpeaks(cjt, 'MinPeakHeight', .1);
        if isempty(locs)
            disp(strcat(num2str(i), ' doesnt have peak torque'));
            continue
        end

        pjtrs(i) = cjt(locs(1)); % peaks are output in order of occurrence
        tpjtrs(i) = kin.raw.stamps(idxs(locs(1))) - behrs(i,5); % relative 
                                                        % to movement onset
           
        %% Find Torque Onset index and time 
        % tO is the index into the torque profile (i.e. cjt) where
        % torque onset is found. Torque onset is currently relative to
        % movement onset. 
        tO = getTorqueOnset(cjt,0,2,THRESHOLD,i,'rs');        
%             torqueOnsetMagjtrs(i) = cjt(tO);
        if ~isnan(tO)
            timeOfTorqueOnsetjtrs(i) = kin.raw.stamps(idxs(tO))-behrs(i,5);
        end

        plotMovementOnset = 0;
        if plotMovementOnset 
            %% If you care to plot where movement onset is, you'll need to
            % implement a tiny bit; instead of using tO = getTorqueOnset as
            % above, compute movementOnset first, and then plug that in place
            % of the ~. 
            % Since kinematics are sampled at 500 Hz, movement onset is
            % given by:
            movementOnset = round((behrs(i,5) - behrs(i,4))*500);        
        end
    end;

    clear lfp*

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    % nanmean(([abatRight(MIchans) abatNE(MIchans) abatTop(MIchans) abatNW(MIchans) abatLeft(MIchans) abatSW(MIchans) abatBottom(MIchans) abatSE(MIchans)]-1000)./1000)
    rsba = [-0.0141    0.0057   -0.0647   -0.0423   -0.0295    0.0336   -0.0056   -0.0112];
    %These are from grpstats(tpjtrs, beh(:,8), 'mean')
    rsjt = [ 0.0300    0.0264    0.0232    0.0325    0.0447    0.0578    0.0263    0.0263];
    %These are time of torque onset: grpstats(timeOfTorqueOnsetjtrs, beh(:,8), 'mean')
%     rsto = [ -0.1432 -0.1616 -0.1762 -0.1680 -0.1570 -0.1653 -0.1599 -0.1403];
    %method 2 of computing torque onset:
    %rsto = [ -0.0871 -0.1005 -0.1205 -0.0894 -0.0757 -0.0841 -0.0964 -0.0903];
    rsto = grpstats(timeOfTorqueOnsetjtrs, behrs(:,8), 'mean')';
    %rsto = grpstats(timeOfTorqueOnsetjtrs(GOODTRIALS), behrs(GOODTRIALS,8), 'mean')';
    
    %TESTING USING ALIGNMENT TO TORQUE ONSET 
%     rsba = [0.0841 0.1074 0.0777 0.0826 0.0638 0.1210 0.1137 0.0985];
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
    
    [B, BINT, ~, ~, STATS] = regress(rsto', [ones(size(rsba')) rsba'])
    disp('      R2, F, P-VALUE, EST. VARIANCE ERROR')
    
    
    % ANCOVA TESTING FOR SIGNIFICANCE; these values were hardcoded from
    % Matt for V13 and V17
%     testx = [rsba'; -0.0828; -0.0850; -0.0683; -0.0057; -0.0375; 0.0096; -0.0602; -0.0682; -0.0720; -0.1602; -0.1180; -0.0758; -0.0470; -0.0516; -0.0418; -0.0623];
%     testy = [rsto'; 0.0340; 0.0133; 0.0057; 0.0445; 0.0511; 0.0590; 0.0468; 0.0600; 0.0345; 0.0147; -0.0200; 0.0145; 0.0370; 0.0403; 0.0686; 0.0514];
%     aoctool(testx,testy,[ones(8,1);ones(16,1)+ones(16,1)])

end
%%
linearModelStatistics.threshold = thresholdRange;
linearModelStatistics.slope = allSlope;
linearModelStatistics.intercept = allIntercept;
linearModelStatistics.R2 = allR2;

plot(allR2)

%% Plot torque onset time as a function of direction 

torqueRight = timeOfTorqueOnsetjtrs(behrs(:,8) == 1);
torqueNE = timeOfTorqueOnsetjtrs(behrs(:,8) == 2);
torqueTop = timeOfTorqueOnsetjtrs(behrs(:,8) == 3);
torqueNW = timeOfTorqueOnsetjtrs(behrs(:,8) == 4);
torqueLeft = timeOfTorqueOnsetjtrs(behrs(:,8) == 5);
torqueSW = timeOfTorqueOnsetjtrs(behrs(:,8) == 6);
torqueBottom = timeOfTorqueOnsetjtrs(behrs(:,8) == 7);
torqueSE = timeOfTorqueOnsetjtrs(behrs(:,8) == 8);


x = [1 * ones(size(torqueRight)); 2 * ones(size(torqueNE)); ...
    3 * ones(size(torqueTop)); 4 * ones(size(torqueNW)); ...
    5 * ones(size(torqueLeft)); 6 * ones(size(torqueSW)); ...
    7 * ones(size(torqueBottom)); 8 * ones(size(torqueSE))];

plot(x,[torqueRight; torqueNE; torqueTop; torqueNW; torqueLeft; torqueSW;...
    torqueBottom; torqueSE]','.')
hold on
plot(1:8, [nanmean(torqueRight); nanmean(torqueNE); nanmean(torqueTop); nanmean(torqueNW); nanmean(torqueLeft); nanmean(torqueSW);...
    nanmean(torqueBottom); nanmean(torqueSE)]','lineWidth',1)
xlabel('Movement Direction (E, NE, N, NW, W, SW, S, SE')
ylabel('Torque Onset Time')
title('Distribution of torque onset time by movement direction (Rockstar)')

%% ANCOVA - replicates figure 2 (compute PARALLEL lines)
    testx = [rsba'; -0.0828; -0.0850; -0.0683; -0.0057; -0.0375; 0.0096; -0.0602; -0.0682; -0.0720; -0.1602; -0.1180; -0.0758; -0.0470; -0.0516; -0.0418; -0.0623];
    testy = [rsjt'; 0.0340; 0.0133; 0.0057; 0.0445; 0.0511; 0.0590; 0.0468; 0.0600; 0.0345; 0.0147; -0.0200; 0.0145; 0.0370; 0.0403; 0.0686; 0.0514];
    aoctool(testx,testy,[ones(8,1);ones(16,1)+ones(16,1)])
