
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%                       PLOTS                          %%%%%%
 %%%%%          Using the 'results' data structure          %%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CONSTANTS 

%Datasets for Youke:
DATASETS = ['M1TM_20111014'; 'M1TM_20111017'; 'M1TM_20111019'; 'M1TM_20111021'; 'M1TM_20111025'];
%Datasets for Big Papi:
% DATASETS = ['M1TM_20101020'; 'M1TM_20101021'; 'M1TM_20101022'; 'M1TM_20101025'; 'M1TM_20101026'; 'M1TM_20101028'];
        


NUM_ELEC = 96;
NUM_OBS = 4001;
%PEAK_FRQ Should be calculated again via power spectrum
PEAK_FRQ = 18;
BANDWIDTH = 3;
SAMPLE_FRQ = 1000;
%THRESHOLD for linear interpolation
THRESHOLD = 0.15;
%Segment of the data to use for computing BAT and BAM
START_TIME = 500;
END_TIME = 3500;   
%%

INDX_FILENAME = 1;
INDX_BAM = 2;
INDX_BAT = 3;
INDX_MODEL = 4;
INDX_R2 = 5;
INDX_PValue = 6;
INDX_VS = 7;
INDX_WVS = 8;
INDX_TRIAL_BIN = 9;

%% PLOT: BAT_PLOT_M1TM_20111014
%Produce a BAT 10x10 plot for all datasets 
for d = 1:size(DATASETS,1)

    %% Preallocate
    plot_this_matrix_BAM = zeros(10,10);
    plot_this_matrix_BAT = zeros(10,10);

    %% BAM Plot Prep
    for i=1:128
        if i <= NUM_ELEC
            index = chan2rc(i,:);
            plot_this_matrix_BAM(index(2),-index(1)+11) = results{d,INDX_BAM}(i);
        end
    end
    %% BAT Plot Prep
    for i=1:128
        if i <= NUM_ELEC
            index = chan2rc(i,:);
            plot_this_matrix_BAT(index(2),-index(1)+11) = results{d,INDX_BAT}(i);
        end
    end
    
    %% Convert 0 to NaN for color plot 
%     results{d,INDX_BAT}(~results{d,INDX_BAT}) = nan; %Turn the zeros into nans so they don't screw up the color plot
%     results{d,INDX_BAM}(~results{d,INDX_BAM}) = nan;
    %% Make plot
    close all;
    hold on;
    axis off;
    axis square;
    axis([1 10 1 10]);

    % caxis([2100 2600]);
    for i=1:10
        for j=1:10
            c = plot_this_matrix_BAM(i,j);
            if ~isnan(c) && c~=0
                scatter(i,j,500,c,'filled');
                elec_num = find(chan2rc(:,2)==i & chan2rc(:,1) == -(j-11));
                %text(i,j,num2str(elec_num)); %Labels which electrode is each plotted circle
            end
        end
    end
    %% Compute BAO and place arrow
    model = results{d, INDX_MODEL};
    coefs = model.Coefficients.Estimate;
    coefs_norm = norm([(coefs(3)/2),(coefs(2)/2)]);
    R2 = results{d,INDX_R2};
    %5.5 is the center of a plot with axis from 1 to 10; 10-1 = 9 / 2 = 4.5 + 1
    % WHY IS THE Y-COORDINATE BEING SUBTRACTED HERE? 
    arrow([5.5,5.5],[5.5+((coefs(3)/2)/coefs_norm)*(R2*4.5),5.5-((coefs(2)/2)/coefs_norm)*(R2*4.5)],'Width',3);

    % set(get(colorbar),'Location','southoutside');
    c_bar = colorbar('southoutside');
    set(get(c_bar, 'xlabel'), 'string', 'time (ms)')

    saveas(gcf,strcat('BAT_PLOT_',results{d,INDX_FILENAME}),'epsc');

    % NOTE: MODEL IS ONLY COMPUTED FOR BAT; NOT BAM. 
     
end

%% PLOT: Vector_strength_M1TM_201110
%BAO for each dataset. Compute Vector Strength. Run after computing BAO for all days
[NUM_DATASETS,~] = size(results);
ang=0:0.01:2*pi; 
r = 1;
xp=r*cos(ang);
yp=r*sin(ang);
plot(0+xp,0+yp, 'Color', 'black');
axis([-1 1 -1 1]);
axis square;
axis off;
angles = zeros(NUM_DATASETS,1);
for d = 1:NUM_DATASETS
    model = results{d, INDX_MODEL};
    coefs = model.Coefficients.Estimate;
    coefs_norm = norm([(coefs(3)/2),(coefs(2)/2)]);
    t_filename = results{d, INDX_FILENAME};
    
    x = ((coefs(3)/2)/coefs_norm);
    y = ((coefs(2)/2)/coefs_norm);
    
    r2 = results{d,INDX_R2};
%     coefs = cell2mat(BAO_coef(i,1));
%     t_filename = char(BAO_coef(i,2));
%     coefs_norm = norm([(coefs(3)/2),(coefs(2)/2)]);
    %index 3 is the R2 of the vector 
%     r2 = cell2mat(BAO_coef(i,3));
    %index 4 is the angle of the vector in radians
    
    angles(d) = atan2(y,x);  
    
    
    arrow([0,0],[0+x*r2,0-y*r2],'Width',1, 'Length', 3);
    if(coefs(3)>0)
        text(0+x*r2+0.02,0-y*r2, t_filename(length(t_filename)-2:length(t_filename)));
    else
        text(0+x*r2-.12,0-y*r2, t_filename(length(t_filename)-2:length(t_filename)));
    end
    
    if d == NUM_DATASETS
        vs = vectorStrength(angles);
        wvs = weightedVectorStrength(angles,cell2mat(results(:,INDX_R2)));
        
        %% Simulation: Compute p-value for Vector Strength
        %Generate random set of NUM_DATASETS vectors and R2 values [theta r2] rand from 0 to 2pi
        NUM_SIMS_p = 100000;
%         NUM_DATASETS = 20;
        simulated_wvs = zeros(NUM_SIMS_p,1);
        simulated_vs = zeros(NUM_SIMS_p,1);
        for i = 1:NUM_SIMS_p
            rand_vectors = [rand([NUM_DATASETS 1])*2*pi rand([NUM_DATASETS 1])];
            sim_vs = vectorStrength(rand_vectors(:,1));
            sim_wvs = weightedVectorStrength(rand_vectors(:,1),rand_vectors(:,2));
            simulated_vs(i,:) = sim_vs;
            simulated_wvs(i,:) = sim_wvs;
        end

        sum_vs = sum(simulated_vs >= vs);
        sum_wvs = sum(simulated_wvs >= wvs);
        p_valueVS = sum_vs / NUM_SIMS_p;
        p_valueWVS = sum_wvs / NUM_SIMS_p;
        %%
        
        text(0,-0.8,strcat('Vector Strength: ',num2str(vs), ' p = ', num2str(p_valueVS)), 'FontName', 'Arial','FontSize',12, 'FontWeight','bold');
        text(0,-.9,strcat('Weighted (R2) Vector Strength: ',num2str(wvs), ' p = ', num2str(p_valueWVS)), 'FontName', 'Arial','FontSize',12, 'FontWeight','bold');
    end
end
saveas(gcf,strcat('Vector_strength_',results{d,INDX_FILENAME}(1:11)),'epsc');

%% PLOT: BAT_avg_time_M1TM_201110
[NUM_DATASETS,~] = size(results);
BAT_avgs = zeros(NUM_DATASETS,1);
BAT_extrema = zeros(NUM_DATASETS,2);

for d = 1:NUM_DATASETS
    BAT = results{d, INDX_BAT};
    try
        BAT(~BAT) = nan;
    catch
        
    end    
    BAT_avgs(d) = nanmean(BAT);
    BAT_extrema(d,:) = [min(BAT) max(BAT)];
end
global_min = min(BAT_extrema(:,1));
dist_to_min = BAT_avgs - BAT_extrema(:,1);
dist_to_max = BAT_avgs - BAT_extrema(:,2);

figure
errorbar(1:NUM_DATASETS, BAT_avgs, dist_to_min, dist_to_max)
ax = gca; 
xlabels = results(:,INDX_FILENAME)';
set(ax,'XTick', 1:NUM_DATASETS)
set(ax,'XTickLabel', xlabels)

text(1,global_min,'min, avg, max BAT', 'FontName', 'Arial','FontSize',12, 'FontWeight','bold');
saveas(gcf,strcat('BAT_avg_time_',results{1,INDX_FILENAME}(1:11)),'epsc');

%% PLOT: BAM_avg_time_M1TM_201110
[NUM_DATASETS,~] = size(results);
BAM_avgs = zeros(NUM_DATASETS,1);
BAM_extrema = zeros(NUM_DATASETS,2);

for d = 1:NUM_DATASETS
    BAM = results{d, INDX_BAM};
    try
        BAM(~BAM) = nan;
    catch
        
    end    
    BAM_avgs(d) = nanmean(BAM);
    BAM_extrema(d,:) = [min(BAM) max(BAM)];
end
global_min = min(BAM_extrema(:,1));
dist_to_min = BAM_avgs - BAM_extrema(:,1);
dist_to_max = BAM_avgs - BAM_extrema(:,2);

figure
errorbar(1:NUM_DATASETS, BAM_avgs, dist_to_min, dist_to_max)
ax = gca; 
xlabels = results(:,INDX_FILENAME)';
set(ax,'XTick', 1:NUM_DATASETS)
set(ax,'XTickLabel', xlabels)

text(1,global_min,'min, avg, max BAM', 'FontName', 'Arial','FontSize',12, 'FontWeight','bold');
saveas(gcf,strcat('BAM_avg_time_',results{1,INDX_FILENAME}(1:11)),'epsc');
%%

%% PLOT: R2_M1TM_201110
figure
bar(cell2mat(results(:,INDX_R2)))
axis([0.5 5.5 0 1]);
ax = gca;
xlabels = results(:,INDX_FILENAME)';
set(ax,'XTick', 1:NUM_DATASETS)
set(ax,'XTickLabel', xlabels)
text(1,.8,'R2 values', 'FontName', 'Arial','FontSize',12, 'FontWeight','bold');
saveas(gcf,strcat('R2_',results{1,INDX_FILENAME}(1:11)),'epsc');

%% PLOT: VS_M1TM_201110
figure
bar(cell2mat(results(:,INDX_VS)))
axis([0.5 5.5 0 1]);
ax = gca;
xlabels = results(:,INDX_FILENAME)';
set(ax,'XTick', 1:NUM_DATASETS)
set(ax,'XTickLabel', xlabels)
text(1,.95,'Vector Strengths (from 1000x simulation)', 'FontName', 'Arial','FontSize',12, 'FontWeight','bold');
saveas(gcf,strcat('VS_',results{1,INDX_FILENAME}(1:11)),'epsc');

%% PLOT: WVS_M1TM_201110
figure
bar(cell2mat(results(:,INDX_WVS)))
axis([0.5 5.5 0 1]);
ax = gca;
xlabels = results(:,INDX_FILENAME)';
set(ax,'XTick', 1:NUM_DATASETS)
set(ax,'XTickLabel', xlabels)
text(1,.99,'Weighted Vector Strengths (from 1000x simulation)', 'FontName', 'Arial','FontSize',12, 'FontWeight','bold');
saveas(gcf,strcat('WVS_',results{1,INDX_FILENAME}(1:11)),'epsc');

%% LEARNING PLOTS: (change this title) R2 of trial bins
[NUM_DATASETS,~] = size(results);

figure
xlabel('Trial bin')
ylabel('R2')
for d = 1:NUM_DATASETS
    trials = results{d,INDX_TRIAL_BIN};
    [NUM_TRIALS,~] = size(trials);
    hold all
    r2 = cell2mat(trials(:,INDX_R2));
    r2(r2==1) = nan;
    plot(r2,'--.','MarkerSize',26)
    ax = gca;
    set(ax,'XTick', 1:NUM_TRIALS)
    
%     set(ax,'XLabel','hi')
    
%     text(1,.8,'R2 values', 'FontName', 'Arial','FontSize',12, 'FontWeight','bold');    
end
legend('show')

%% LEARNING PLOTS: (change this title) R2 of trial bins with DAYS on X
[NUM_DATASETS,~] = size(results);
trial_data = [];

figure
xlabel('Day')
ylabel('R2')
for d = 1:NUM_DATASETS
    trials = results{d,INDX_TRIAL_BIN};
    [NUM_TRIALS,~] = size(trials);
    for t = 1:NUM_TRIALS
        trial_data(t,d) = cell2mat(trials(t,INDX_R2));
    end
end
% trial_data(~trial_data) = nan;
% trial_data(trial_data==1) = nan;
for t = 1:NUM_TRIALS
    hold all
    plot(trial_data(t,:), '.','MarkerSize', 26)
end
ax = gca;
set(ax,'XTick', 1:NUM_DATASETS)
legend('show')

%% LEARNING PLOTS: (change this title) VS of trial bins
[NUM_DATASETS,~] = size(results);

figure
xlabel('Trial bin')
ylabel('VS')
for d = 1:NUM_DATASETS
    trials = results{d,INDX_TRIAL_BIN};
    [NUM_TRIALS,~] = size(trials);
    hold all
    vs = cell2mat(trials(:,INDX_VS));
%     vs(vs==1) = nan;
    plot(vs,'--.','MarkerSize',26)
    ax = gca;
    set(ax,'XTick', 1:NUM_TRIALS)
    
%     set(ax,'XLabel','hi')
    
%     text(1,.8,'R2 values', 'FontName', 'Arial','FontSize',12, 'FontWeight','bold');    
end
legend('show')

%% LEARNING PLOTS: (change this title) VS of trial bins with DAYS on X
[NUM_DATASETS,~] = size(results);
trial_data = [];

figure
xlabel('Day')
ylabel('VS')
for d = 1:NUM_DATASETS
    trials = results{d,INDX_TRIAL_BIN};
    [NUM_TRIALS,~] = size(trials);
    for t = 1:NUM_TRIALS
        trial_data(t,d) = cell2mat(trials(t,INDX_VS));
    end
end
trial_data(trial_data==0) = nan;
% trial_data(trial_data==1) = nan;
for t = 1:NUM_TRIALS
    hold all
    plot(trial_data(t,:), '.','MarkerSize', 26)
end
ax = gca;
set(ax,'XTick', 1:NUM_DATASETS)
legend('show')



%% PLOT add-on: Plot all bootstrapped arrows
for j = 1:NUM_SIMS
    plotX = ones(96,3);
    plotX(:,2) = chan2rc(1:96,1);
    plotX(:,3) = chan2rc(1:96,2);
    coefs = regress(bootstrapped_BAM(:,j),plotX);
    coefs_norm = norm([(coefs(3)/2),(coefs(2)/2)]);
    arrow([5.5,5.5],[5.5+(coefs(3)/2)/coefs_norm*3,5.5-(coefs(2)/2)/coefs_norm*3],'Width',3);
end
% saveas(gcf,strcat('BAT_PLOT_BOOTSTRAP_',filename','eps')

%%
%scatter(chan2rc(1:96,2), -chan2rc(1:96,1) + 11, 50 * ones(96,1), beta_attn_times(:), 'filled')
%%The above line plots everything in the for loop in a better way.

%PLOTTING CRITERIA
%Always label axes
%Consider the aspect ratio (make it square)
%Create a legend.

%Do pairwise comparisons of beta attenuation times (not spacially) with
%every other data set

%[DONE] Replot maps so they all have same color scale
%2 beta attenuations times for each electrode; plot first component on x
%axis, second component on y axis

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


%To quantify circularity vs unimodality, use "vector strength"
%formula by Goldberg 1969ish. 
%Look at the power spectrum across all electrodes to find the peak in beta
%(make sure it's 18). Vim's
%book. 



%% PLOT: BAT_PLOT_M1TM_20111014
%Produce a BAT 10x10 plot for all datasets 
for d = 1:size(DATASETS,1)

    %% Preallocate
    plot_this_matrix_BAM = zeros(10,10);
    plot_this_matrix_BAT = zeros(10,10);

    %% BAM Plot Prep
    for i=1:128
        if i <= NUM_ELEC
            index = chan2rc(i,:);
            plot_this_matrix_BAM(index(2),-index(1)+11) = results{d,INDX_BAM}(i);
        end
    end
    %% BAT Plot Prep
    for i=1:128
        if i <= NUM_ELEC
            index = chan2rc(i,:);
            plot_this_matrix_BAT(index(2),-index(1)+11) = results{d,INDX_BAT}(i);
        end
    end
    
    %% Convert 0 to NaN for color plot 
%     results{d,INDX_BAT}(~results{d,INDX_BAT}) = nan; %Turn the zeros into nans so they don't screw up the color plot
%     results{d,INDX_BAM}(~results{d,INDX_BAM}) = nan;
    %% Make plot
    close all;
    hold on;
    axis off;
    axis square;
    axis([1 10 1 10]);
    cmin = min(results{d,INDX_BAT});
    cmax = max(results{d,INDX_BAT});
    
    caxis([cmin cmax]);
    colormap(flipud(gray));
    set(gcf,'color','w');
    for i=1:10
        for j=1:10
            c = plot_this_matrix_BAT(i,j);
            if ~isnan(c) && c~=0
                scatter(i,j,1500,c,'filled');
                elec_num = find(chan2rc(:,2)==i & chan2rc(:,1) == -(j-11));
                %text(i,j,num2str(elec_num)); %Labels which electrode is each plotted circle
            end
        end
    end
    %% Compute BAO and place arrow
    model = results{d, INDX_MODEL};
    coefs = model.Coefficients.Estimate;
    coefs_norm = norm([(coefs(3)/2),(coefs(2)/2)]);
    R2 = results{d,INDX_R2};
    %5.5 is the center of a plot with axis from 1 to 10; 10-1 = 9 / 2 = 4.5 + 1
    % WHY IS THE Y-COORDINATE BEING SUBTRACTED HERE? 
    arrow([5.5,5.5],[5.5+((coefs(3)/2)/coefs_norm)*(R2*4.5),5.5-((coefs(2)/2)/coefs_norm)*(R2*4.5)],'Width',3);

    % set(get(colorbar),'Location','southoutside');\
    c_bar = colorbar('southoutside');
    set(get(c_bar, 'xlabel'), 'string', 'time (ms)')

    saveas(gcf,strcat('BAT_PLOT_',results{d,INDX_FILENAME}),'epsc');

    % NOTE: MODEL IS ONLY COMPUTED FOR BAT; NOT BAM. 
     
end

