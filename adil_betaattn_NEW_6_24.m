clearvars -except model* BAO_coef* chan2rc
% clearvars -except BAT_R2 BAT_M11014 BAT_M11014_boostrap MAT_M11014 BAT_M11017_boostrap MAT_M11017 BAT_M11019_boostrap MAT_M11019 BAT_M11021_boostrap MAT_M11021 BAT_M11025_boostrap MAT_M11025 BAT_S11014 MAT_S11014;
%Change file names at line 9, 149, 204, 215, 218,219
%'M1TM_20111014.mat'
%'M1TM_20111017.mat'
%'M1TM_20111019.mat'
%'M1TM_20111021.mat'
%'M1TM_20111025.mat'
filename = 'M1TM_20111025';
NUM_ELEC = 96;
NUM_OBS = 4001;
PEAK_FRQ = 18;
BANDWIDTH = 3;
SAMPLE_FRQ = 1000;
THRESHOLD = 0.15;
START_TIME = 500;
END_TIME = 3500;
% This not needed because it is computed later. NUM_DATASETS = 5;

%% Build Data
all_electrodes = build_data(strcat(filename,'.mat'), NUM_ELEC, NUM_OBS, PEAK_FRQ, BANDWIDTH, SAMPLE_FRQ);
NUM_TRIALS = size(all_electrodes,2);
%% Compute Power Spectrum (empty)

%% Simulation: Simulate
NUM_SIMS = 20;
BIN = 30;
TRIAL_RANGE = 1:NUM_TRIALS;
bootstrapped_times = simulate_BAT(all_electrodes,NUM_ELEC, TRIAL_RANGE, NUM_SIMS, BIN, START_TIME, END_TIME,[82]);

% 
% 
% bootstrapped_times = zeros(NUM_ELEC,NUM_SIMS);
% if BIN > max(TRIAL_RANGE) - min(TRIAL_RANGE)
%     error('BIN must be smaller than the number of trials')
% end
% for j = 1:NUM_SIMS
%     i = round(rand(BIN,1)*(max(TRIAL_RANGE)-min(TRIAL_RANGE)) + min(TRIAL_RANGE));
%     BAM = compute_BAM(all_electrodes(:,i,:), NUM_ELEC, START_TIME, END_TIME, [82]);
%     bootstrapped_times(:,j) = BAM;    
%     j
% end        
% bootstrapped_times = simulate_BAT(all_electrodes,NUM_ELEC, 0:10, 1000,10,[82]);
%%

%% PLOT: To visualize (debugging purposes), plot BAT
for i = 1:96
    data = nanmean(all_electrodes(500:3500, :,i),2);
    block_test = round((3500 - 500) / 3);
    logitc = @(x, c) (c(1) - c(2)) ./ (1 + exp(-c(3) .* (x - c(4)))) + c(2);
    plot(data,'LineWidth',4) 
    LB = [min(data(1:block_test)) min(data(block_test*2:block_test*3)) -1, 1];
    UB = [max(data(1:block_test)) max(data(block_test*2:block_test*3)) -.01, 3000];
    x = 500:3500;
    clear cHat
    if(~isnan(sum(data)))
        cHat = fitSigmoid2(x, data, LB, UB);
        hold on
        plot(logitc(x, cHat), 'r')
        line([beta_attn_med(i) - START_TIME beta_attn_med(i) - START_TIME],get(gca,'ylim'),'Color','g');
        beta_attn_med(i)
%         line([cHat(4) - START_TIME cHat(4) - START_TIME],get(gca,'ylim'),'Color','g');
%         beta_attn_med(i)
%         cHat(4)
        %THE FOURTH ARGUMENT OF cHat IS MY NEW BAT
        hold off
    end
%     plot(nanmean(all_data(500:3500,:,i),2)) 
    %The above gathers the mean of the 2nd dimension of the matrix, the
    %trials, for each row, i.e. for each second of time. 
    line([1500 1500],get(gca,'ylim'));
    line([beta_attn_times(i)-START_TIME beta_attn_times(i)-START_TIME],get(gca,'ylim'),'Color','r');
    xlim([0,3000])
%     saveas(gcf,'beta_attn_sample_elec1','eps');
    pause;
end
%% Compute BAT for actual experiment
beta_attn_times = compute_BAT(all_electrodes, NUM_ELEC, START_TIME, END_TIME, THRESHOLD, [82]);
%% Compute BAM for actual experiment
beta_attn_med = compute_BAM(all_electrodes(:,:,:), NUM_ELEC, START_TIME, END_TIME, [82]);
%% Simulation: Average each electrode across all trials
simulated_fit = mean(bootstrapped_times,2); %USE NANMEAN?
%% Check if inside confidence interval
simulated_sd = std(bootstrapped_times);
is_fitted = sum(beta_attn_times > quantile(bootstrapped_times, 0.025, 2) & beta_attn_times < quantile(bootstrapped_times, .975, 2));
num_working_electrodes = sum(beta_attn_times > 0);
is_fitted/num_working_electrodes
%% Chan2RC preparation 1/2 (needs internet)
chan2rc = makechan2rc_mac('y','mio');
%% Chan2RC preparation 2/2 (doesn't need internet)
plot_this_matrix = zeros(10,10);
plot_this_matrix_bootstrap = zeros(10,10,NUM_SIMS);
%%
%currentBAT = BAT_M11014; %Use when you store Beta attenuation times so you
%don't have to produce the same matrix over and over.
%beta_attn_times = currentBAT;
%%
%Create matrix to plot with beta attenuation times
for i=1:128
    if i <= NUM_ELEC
        index = chan2rc(i,:);
        plot_this_matrix(index(2),-index(1)+11) = beta_attn_times(i);
    end
end
%%
%Create matrix to plot with bootstrapped beta attenuation times
for i=1:128
    for j=1:NUM_SIMS
        if i <= NUM_ELEC
            index = chan2rc(i,:);
            plot_this_matrix_bootstrap(index(2),-index(1)+11,j) = bootstrapped_times(i,j);
        end
    end
end
        
%% THE COMMAND HERE SHOULD ONLY BE RUN ONCE 
beta_attn_times(~beta_attn_times) = nan; %Turn the zeros into nans so they don't screw up the color plot
beta_attn_med(~beta_attn_med) = nan;
bootstrapped_times(~bootstrapped_times) = nan;
%% Fit Linear Model and Compute Statistics
%Do I need to do this for the bootstrapped times?
plotX = ones(96,3);
plotX(:,2) = chan2rc(1:96,1);
plotX(:,3) = chan2rc(1:96,2);
ds = dataset(beta_attn_times, plotX(:,2),plotX(:,3), 'VarNames', {'BAT','X','Y'});
model = strcat('model_',filename);
eval(strcat(model, ' = LinearModel.fit(ds, ''BAT~X+Y'');'));
R2 = eval(strcat(model,'.Rsquared.Ordinary'));
F = (eval(strcat(model,'.SSR / ',model,'.NumPredictors', '/ ',model,'.MSE')));
pValue = 1 - fcdf(F, eval(strcat(model,'.NumPredictors')), eval(strcat(model,'.DFE')));

%% Fit Linear Model for Bootstrapped Data
plotB = ones(96,3);
plotB(:,2) = chan2rc(1:96,1);
plotB(:,3) = chan2rc(1:96,2);
for i = 1:NUM_SIMS
    i
    ds_b = dataset(bootstrapped_times(:,i), plotB(:,2),plotB(:,3), 'VarNames', {'BAT','X','Y'});
    simulated_model = LinearModel.fit(ds_b, 'BAT~X+Y');
    R2_b = simulated_model.Rsquared.Ordinary;
    coefs_b = simulated_model.Coefficients.Estimate;
    coefs_b_norm = norm([(coefs_b(3)/2),(coefs_b(2)/2)]);
    if ~exist('BAO_b_coef','var')
        BAO_b_coef = {coefs_b filename R2_b 'angle'};
    else
        BAO_b_coef = [BAO_b_coef; {coefs_b filename R2_b 'angle'}];
    end
%     %I don't think I need F and p-value, so I'm leaving them out. 
%     
%     model_b = strcat('model_b_',filename);
%     eval(strcat(model_b, ' = LinearModel.fit(ds_b, ''BAT~X+Y'');'));
%     R2_b = eval(strcat(model_b,'.Rsquared.Ordinary'));
%     F = (eval(strcat(model_b,'.SSR / ',model_b,'.NumPredictors', '/ ',model_b,'.MSE')));
%     pValue = 1 - fcdf(F, eval(strcat(model_b,'.NumPredictors')), eval(strcat(model_b,'.DFE')));
end





%% PLOT: Compute BAO and plot 10x10 square color-plot of BAT
close all;
hold on;
axis off;
axis square;
axis([1 10 1 10]);

% caxis([2100 2600]);
for i=1:10
    for j=1:10
        c = plot_this_matrix(i,j);
        if ~isnan(c) && c~=0
            scatter(i,j,500,c,'filled');
            elec_num = find(chan2rc(:,2)==i & chan2rc(:,1) == -(j-11));
            %text(i,j,num2str(elec_num)); %Labels which electrode is each plotted circle
        end
    end
end
coefs = eval(strcat(model,'.Coefficients.Estimate;'));
coefs_norm = norm([(coefs(3)/2),(coefs(2)/2)]);
%5.5 is the center of a plot with axis from 1 to 10; 10-1 = 9 / 2 = 4.5 + 1
% WHY IS THE Y-COORDINATE BEING SUBTRACTED HERE? 
arrow([5.5,5.5],[5.5+((coefs(3)/2)/coefs_norm)*(R2*4.5),5.5-((coefs(2)/2)/coefs_norm)*(R2*4.5)],'Width',3);

if ~exist('BAO_coef','var')
    BAO_coef = {coefs filename R2 'angle'};
else
    BAO_coef = [BAO_coef; {coefs filename R2 'angle'}];
end

% set(get(colorbar),'Location','southoutside');\
c_bar = colorbar('southoutside');
set(get(c_bar, 'xlabel'), 'string', 'time (ms)')


saveas(gcf,strcat('BAT_PLOT_',filename),'epsc');
%% PLOT add-on: Plot all bootstrapped arrows
for j = 1:NUM_SIMS
    plotX = ones(96,3);
    plotX(:,2) = chan2rc(1:96,1);
    plotX(:,3) = chan2rc(1:96,2);
    coefs = regress(bootstrapped_times(:,j),plotX);
    coefs_norm = norm([(coefs(3)/2),(coefs(2)/2)]);
    arrow([5.5,5.5],[5.5+(coefs(3)/2)/coefs_norm*3,5.5-(coefs(2)/2)/coefs_norm*3],'Width',3);
end
% saveas(gcf,strcat('BAT_PLOT_BOOTSTRAP_',filename','eps')
%%
% BAT_M11025 = beta_attn_times;
% BAT_M11025_boostrap = bootstrapped_times;

%%
%scatter(chan2rc(1:96,2), -chan2rc(1:96,1) + 11, 50 * ones(96,1), beta_attn_times(:), 'filled')
%%The above line plots everything in the for loop in a better way.

%Use ctrl+d to see what's under the hood of matlab

%PLOTTING CRITERIA
%Always label axes
%Consider the aspect ratio (make it square)
%Create a legend.

%Do pairwise comparisons of beta attenuation times (not spacially) with
%every other data set

%[DONE] Replot maps so they all have same color scale
%2 beta attenuations times for each electrode; plot first component on x
%axis, second component on y axis
%Our figures are on the bottoms; 6 total scatter plots.
%IF NO ELECTRODE, DON'T PLOT IT. 
%close all;
%hold on;
%axis on;
%axis square;
%axis([1500 2500 1500 2500]);
%scatter(BAT_M11021,BAT_M11025);
%h = lsline;
%p2 = polyfit(get(h,'xdata'),get(h,'ydata'),1)
%text(2300,1600,['R=' num2str(p2(1))]);
%title('M11021 vs M11025');
%saveas(gcf,'M11021 M11025','eps')

%arrow([1500,1500],[1500+1*1000,1500+0.5733*1000],'Width',5,'Length',30);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%=======Due Feb 20=======%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Try replotting the arrow for every X trials (to test learning
%%hypothesis)
%Where X is determined by how small you can get while still seeing the
%sigmoidal shape in beta amplitude (so we can nicely filter where beta attn
%"happens"

%If we do need 50 trials, do a boostrap analysis, draw 1000 random samples
%of 50 trials. (or draw 50 random trials, go through the whole motion)

%To quantify circularity vs unimodality, use "vector strength"
%Grab Matt's code or code yourself. Google it, one line mathematical
%formula by Goldberg 1969ish. Or Matt/code folder
%Look at the power spectrum across all electrodes to find the peak in beta
%(make sure it's 18). Figure out how to compute a power spectrum. Vim's
%book. 

%R IS NOT SLOPE... IT'S SQRT OF R2
%Check correlations




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Feb 27
%Compute the arrow of BAT and test the unimodality 

%% PLOT: SIMULATION BAO. Compute Vector Strength. Run after SIMULATION
%NOTE: Should we weight these by R2 value, when computing vector strength?
%TOTAL_SIM is NUM_SIMS * NUM_DATASETS
% [NUM_TOTAL_SIM,~] = size(BAO_coef);

%Plot a circle
ang=0:0.01:2*pi; 
r = 1;
xp=r*cos(ang);
yp=r*sin(ang);
plot(0+xp,0+yp, 'Color', 'black');
axis([-1 1 -1 1]);
axis square;
axis off;
%end circle plot

for i = 1:NUM_SIMS
    coefs_b = cell2mat(BAO_b_coef(i,1));
    t_filename = char(BAO_b_coef(i,2));
    coefs_b_norm = norm([(coefs_b(3)/2),(coefs_b(2)/2)]);
    x = ((coefs_b(3)/2)/coefs_b_norm);
    y = ((coefs_b(2)/2)/coefs_b_norm);
    %index 3 is the R2 of the vector 
    r2_b = cell2mat(BAO_b_coef(i,3));
    %index 4 is the angle of the vector in radians
    BAO_b_coef(i,4) = {acos(x)};  
    arrow([0,0],[0+x*r2_b,0-y*r2_b],'Width',1, 'Length', 3);
    if(coefs_b(3)>0)
        text(0+x*r2_b+0.02,0-y*r2_b, t_filename(length(t_filename)-2:length(t_filename)));
    else
        text(0+x*r2_b-.12,0-y*r2_b, t_filename(length(t_filename)-2:length(t_filename)));
    end
    
    if i == NUM_SIMS
        vs = vectorStrength(cell2mat(BAO_b_coef(:,4)));
        wvs = weightedVectorStrength(cell2mat(BAO_b_coef(:,4)),cell2mat(BAO_b_coef(:,3)));
        text(0,-0.8,strcat('Vector Strength: ',num2str(vs)), 'FontName', 'Arial','FontSize',12, 'FontWeight','bold');
        text(0,-.9,strcat('Weighted (R2) Vector Strength: ',num2str(wvs)), 'FontName', 'Arial','FontSize',12, 'FontWeight','bold');
    end
end

%% PLOT: BAO for each dataset. Compute Vector Strength. Run after computing BAO for all days
%NOTE: Should we weight these by R2 value, when computing vector strength?
[NUM_DATASETS,~] = size(BAO_coef);
ang=0:0.01:2*pi; 
r = 1;
xp=r*cos(ang);
yp=r*sin(ang);
plot(0+xp,0+yp, 'Color', 'black');
axis([-1 1 -1 1]);
axis square;
axis off;
for i = 1:NUM_DATASETS
    coefs = cell2mat(BAO_coef(i,1));
    t_filename = char(BAO_coef(i,2));
    coefs_norm = norm([(coefs(3)/2),(coefs(2)/2)]);
    x = ((coefs(3)/2)/coefs_norm);
    y = ((coefs(2)/2)/coefs_norm);
    %index 3 is the R2 of the vector 
    r2 = cell2mat(BAO_coef(i,3));
    %index 4 is the angle of the vector in radians
    BAO_coef(i,4) = {acos(x)};  
    arrow([0,0],[0+x*r2,0-y*r2],'Width',1, 'Length', 3);
    if(coefs(3)>0)
        text(0+x*r2+0.02,0-y*r2, t_filename(length(t_filename)-2:length(t_filename)));
    else
        text(0+x*r2-.12,0-y*r2, t_filename(length(t_filename)-2:length(t_filename)));
    end
    
    if i == NUM_DATASETS
        vs = vectorStrength(cell2mat(BAO_coef(:,4)));
        wvs = weightedVectorStrength(cell2mat(BAO_coef(:,4)),cell2mat(BAO_coef(:,3)));
        text(0,-0.8,strcat('Vector Strength: ',num2str(vs)), 'FontName', 'Arial','FontSize',12, 'FontWeight','bold');
        text(0,-.9,strcat('Weighted (R2) Vector Strength: ',num2str(wvs)), 'FontName', 'Arial','FontSize',12, 'FontWeight','bold');
    end
end

%% Simulation: Compute p-value for Vector Strength
%Generate random set of NUM_DATASETS vectors and R2 values [theta r2] rand from 0 to 2pi
for j= 1:3
    NUM_SIMS = 100000;
    NUM_DATASETS = 20;
    simulated_wvs = zeros(NUM_SIMS,1);
    simulated_vs = zeros(NUM_SIMS,1);
    for i = 1:NUM_SIMS
        rand_vectors = [rand([NUM_DATASETS 1])*2*pi rand([NUM_DATASETS 1])];
        sim_vs = vectorStrength(rand_vectors(:,1));
        sim_wvs = weightedVectorStrength(rand_vectors(:,1),rand_vectors(:,2));
        simulated_vs(i,:) = sim_vs;
        simulated_wvs(i,:) = sim_wvs;
    end

    sum_vs = sum(simulated_vs >= vs);
    sum_wvs = sum(simulated_wvs >= wvs);
    p_valueVS = sum_vs / NUM_SIMS;
    p_valueWVS = sum_wvs / NUM_SIMS;

    [p_valueVS p_valueWVS]
end
%%





