%'results' is of the form:
% {'filename' BAM BAT model R2 PValue VS WVS {nx8 cell}}
%'{nx8 cell}' is of the form:
% {'set of trials' BAM BAT model R2 PValue VS WVS}

%List of Datasets:
%Datasets for Youke:
DATASETS = ['M1TM_20111014'; 'M1TM_20111017'; 'M1TM_20111019'; 'M1TM_20111021'; 'M1TM_20111025'];
%Datasets for Big Papi:
% DATASETS = [];

results = cell(size(DATASETS,1),9);

%Eventually want to add P-Values for weighted and unweighted vector
%strengths (currently simulating to calculate them every time).

%Indices for the cell of results.
INDX_FILENAME = 1;
INDX_BAM = 2;
INDX_BAT = 3;
INDX_MODEL = 4;
INDX_R2 = 5;
INDX_PValue = 6;
INDX_VS = 7;
INDX_WVS = 8;
INDX_TRIAL_BIN = 9;

%Remember to initialize results to have the # of indices created. 

%%
for d = 1:size(DATASETS,1)
    %% FOR TESTING PURPOSES
%     d = 1;
    %% Set Constants
    clearvars -except d results chan2rc DATASETS INDX*
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
    filename = DATASETS(d,:);
    results{d,INDX_FILENAME} = filename;
    
    %% Build Data
    all_electrodes = build_data(strcat(filename,'.mat'), NUM_ELEC, NUM_OBS, PEAK_FRQ, BANDWIDTH, SAMPLE_FRQ);
    NUM_TRIALS = size(all_electrodes,2);
    
    %% Compute empirical BAM
    'Computing BAM'
    beta_attn_med = compute_BAM(all_electrodes(:,:,:), NUM_ELEC, START_TIME, END_TIME, [82]);
    
    %% Compute empirical BAT
    'Computing BAT'
    beta_attn_times = compute_BAT(all_electrodes, NUM_ELEC, START_TIME, END_TIME, THRESHOLD, [82]);
    
    %% Add BAM and BAT to results
    
    results{d,INDX_BAM} = beta_attn_med;
    results{d,INDX_BAT} = beta_attn_times;
    
    %% Chan2RC preparation 1/2 (needs internet)
%     chan2rc = makechan2rc_mac('y','mio');
    
    %% Chan2RC preparation 2/2 (doesn't need internet)
    plot_this_matrix = zeros(10,10);
    
    %% BAT Plot Prep
    for i=1:128
        if i <= NUM_ELEC
            index = chan2rc(i,:);
            plot_this_matrix(index(2),-index(1)+11) = beta_attn_times(i);
        end
    end

    %% Convert 0 to NaN for color plot 
    beta_attn_times(~beta_attn_times) = nan; %Turn the zeros into nans so they don't screw up the color plot
    beta_attn_med(~beta_attn_med) = nan;
    
    %% Fit Linear Model and Compute Statistics
    'Fitting Linear Model'
    plotX = ones(96,3);
    plotX(:,2) = chan2rc(1:96,1);
    plotX(:,3) = chan2rc(1:96,2);
    ds = dataset(beta_attn_times, plotX(:,2),plotX(:,3), 'VarNames', {'BAT','X','Y'});
    model = strcat('model_',filename);
    eval(strcat(model, ' = LinearModel.fit(ds, ''BAT~X+Y'');'));
    R2 = eval(strcat(model,'.Rsquared.Ordinary'));
    F = (eval(strcat(model,'.SSR / ',model,'.NumPredictors', '/ ',model,'.MSE')));
    pValue = 1 - fcdf(F, eval(strcat(model,'.NumPredictors')), eval(strcat(model,'.DFE')));

    %% Add Linear Model to results
    results{d,INDX_MODEL} = eval(model);
    results{d,INDX_R2} = R2;
    results{d,INDX_PValue} = pValue;    
    
    %% Simulation: Simulate to compute Vector Strength of empirical 
    'Simulating'
    NUM_SIMS = 1000;
    BIN = 50;
    TRIAL_RANGE = 1:NUM_TRIALS;
    %Simulation actually calculates BAM; estimation is necessary with
    %smaller BINS. 
    bootstrapped_BAM = simulate_BAT(all_electrodes,NUM_ELEC, TRIAL_RANGE, NUM_SIMS, BIN, START_TIME, END_TIME,[82]);
    
    %% Convert zeros to NaN
    bootstrapped_BAM(~bootstrapped_BAM) = nan;
    
    %% BAT Bootstrap Plot Prep
    plot_this_matrix_bootstrap = zeros(10,10,NUM_SIMS);
    for i=1:128
        for j=1:NUM_SIMS
            if i <= NUM_ELEC
                index = chan2rc(i,:);
                plot_this_matrix_bootstrap(index(2),-index(1)+11,j) = bootstrapped_BAM(i,j);
            end
        end
    end
    
    %% Fit Linear Model for Bootstrapped Data
    'Fitting Linear Model for Bootstrapped Data'
    plotB = ones(96,3);
    plotB(:,2) = chan2rc(1:96,1);
    plotB(:,3) = chan2rc(1:96,2);
    for i = 1:NUM_SIMS
        ds_b = dataset(bootstrapped_BAM(:,i), plotB(:,2),plotB(:,3), 'VarNames', {'BAT','X','Y'});
        simulated_model = LinearModel.fit(ds_b, 'BAT~X+Y');
        R2_b = simulated_model.Rsquared.Ordinary;
        coefs_b = simulated_model.Coefficients.Estimate;
        coefs_b_norm = norm([(coefs_b(3)/2),(coefs_b(2)/2)]);
        if ~exist('BAO_b_coef','var')
            BAO_b_coef = {coefs_b filename R2_b 'angle'};
        else
            BAO_b_coef = [BAO_b_coef; {coefs_b filename R2_b 'angle'}];
        end
    end
    
    %% Compute VS and WVS from Simulation.
    'Computing VS and WVS from Simulation'
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
        if i == NUM_SIMS
            vs = vectorStrength(cell2mat(BAO_b_coef(:,4)));
            wvs = weightedVectorStrength(cell2mat(BAO_b_coef(:,4)),cell2mat(BAO_b_coef(:,3)));
        end
    end
    
    %% Add Vector Strength and Weighted Vector Strength to results
    results{d,INDX_VS} = vs;
    results{d,INDX_WVS} = wvs;
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%% Repeat Analysis by breaking data into trial bins %%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BIN = 50;
    endpnts = 1:BIN:NUM_TRIALS;
    results{d,INDX_TRIAL_BIN} = cell(floor(NUM_TRIALS/BIN),8);
    %%
    for k = 1:size(endpnts,2)-1
        binned_trials = endpnts(k):endpnts(k+1);
        results{d,INDX_TRIAL_BIN}{k,INDX_FILENAME} = strcat('Trials',num2str(binned_trials));
        
        %%
        %% Compute empirical BAM
        'Computing BAM'
        beta_attn_med = compute_BAM(all_electrodes(:,binned_trials,:), NUM_ELEC, START_TIME, END_TIME, [82]);

        %% Compute empirical BAT
        'Computing BAT'
        beta_attn_times = compute_BAT(all_electrodes(:,binned_trials,:), NUM_ELEC, START_TIME, END_TIME, THRESHOLD, [82]);

        %% Add BAM and BAT to results

        results{d,INDX_TRIAL_BIN}{k,INDX_BAM} = beta_attn_med;
        results{d,INDX_TRIAL_BIN}{k,INDX_BAT} = beta_attn_times;

        %% Chan2RC preparation 1/2 (needs internet)
%         chan2rc = makechan2rc_mac('y','mio');

        %% Chan2RC preparation 2/2 (doesn't need internet)
        plot_this_matrix = zeros(10,10);

        %% BAT Plot Prep
        for i=1:128
            if i <= NUM_ELEC
                index = chan2rc(i,:);
                plot_this_matrix(index(2),-index(1)+11) = beta_attn_times(i);
            end
        end

        %% Convert 0 to NaN for color plot 
        beta_attn_times(~beta_attn_times) = nan; %Turn the zeros into nans so they don't screw up the color plot
        beta_attn_med(~beta_attn_med) = nan;

        %% Fit Linear Model and Compute Statistics
        'Fitting Linear Model'
        plotX = ones(96,3);
        plotX(:,2) = chan2rc(1:96,1);
        plotX(:,3) = chan2rc(1:96,2);
        ds = dataset(beta_attn_times, plotX(:,2),plotX(:,3), 'VarNames', {'BAT','X','Y'});
        model = strcat('model_',filename);
        eval(strcat(model, ' = LinearModel.fit(ds, ''BAT~X+Y'');'));
        R2 = eval(strcat(model,'.Rsquared.Ordinary'));
        F = (eval(strcat(model,'.SSR / ',model,'.NumPredictors', '/ ',model,'.MSE')));
        pValue = 1 - fcdf(F, eval(strcat(model,'.NumPredictors')), eval(strcat(model,'.DFE')));

        %% Add Linear Model to results
        results{d,INDX_TRIAL_BIN}{k,INDX_MODEL} = eval(model);
        results{d,INDX_TRIAL_BIN}{k,INDX_R2} = R2;
        results{d,INDX_TRIAL_BIN}{k,INDX_PValue} = pValue;    

        %% Simulation: Simulate
        'Simulating'
        NUM_SIMS = 1000;
        BIN = 30;
        TRIAL_RANGE = binned_trials;
        %Simulation actually calculates BAM; estimation is necessary with
        %smaller BINS. 
        bootstrapped_BAM = simulate_BAT(all_electrodes,NUM_ELEC, TRIAL_RANGE, NUM_SIMS, BIN, START_TIME, END_TIME,[82]);

        %% Convert zeros to NaN
        bootstrapped_BAM(~bootstrapped_BAM) = nan;

        %% BAT Bootstrap Plot Prep
        plot_this_matrix_bootstrap = zeros(10,10,NUM_SIMS);
        for i=1:128
            for j=1:NUM_SIMS
                if i <= NUM_ELEC
                    index = chan2rc(i,:);
                    plot_this_matrix_bootstrap(index(2),-index(1)+11,j) = bootstrapped_BAM(i,j);
                end
            end
        end

        %% Fit Linear Model for Bootstrapped Data
        'Fitting Linear Model for Bootstrapped Data'
        plotB = ones(96,3);
        plotB(:,2) = chan2rc(1:96,1);
        plotB(:,3) = chan2rc(1:96,2);
        for i = 1:NUM_SIMS
            ds_b = dataset(bootstrapped_BAM(:,i), plotB(:,2),plotB(:,3), 'VarNames', {'BAT','X','Y'});
            simulated_model = LinearModel.fit(ds_b, 'BAT~X+Y');
            R2_b = simulated_model.Rsquared.Ordinary;
            coefs_b = simulated_model.Coefficients.Estimate;
            coefs_b_norm = norm([(coefs_b(3)/2),(coefs_b(2)/2)]);
            if ~exist('BAO_bin_coef','var')
                BAO_bin_coef = {coefs_b filename R2_b 'angle'};
            else
                BAO_bin_coef = [BAO_bin_coef; {coefs_b filename R2_b 'angle'}];
            end
        end

        %% Compute VS and WVS from Simulation.
        'Computing VS and WVS from Simulation'
        for i = 1:NUM_SIMS
            coefs_b = cell2mat(BAO_bin_coef(i,1));
            t_filename = char(BAO_bin_coef(i,2));
            coefs_b_norm = norm([(coefs_b(3)/2),(coefs_b(2)/2)]);
            x = ((coefs_b(3)/2)/coefs_b_norm);
            y = ((coefs_b(2)/2)/coefs_b_norm);
            %index 3 is the R2 of the vector 
            r2_b = cell2mat(BAO_bin_coef(i,3));
            %index 4 is the angle of the vector in radians
            BAO_bin_coef(i,4) = {acos(x)};      
            if i == NUM_SIMS
                vs = vectorStrength(cell2mat(BAO_bin_coef(:,4)));
                wvs = weightedVectorStrength(cell2mat(BAO_bin_coef(:,4)),cell2mat(BAO_bin_coef(:,3)));
            end
        end
        clear BAO_bin_coef

        %% Add Vector Strength and Weighted Vector Strength to results
        results{d,INDX_TRIAL_BIN}{k,INDX_VS} = vs;
        results{d,INDX_TRIAL_BIN}{k,INDX_WVS} = wvs;

    end
    
    
    %%
    strcat('Finished Dataset ', filename, ' Press any key to continue.')
    pause
end
%%

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%                       PLOTS                      %%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CONSTANTS 
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
            c = plot_this_matrix_BAT(i,j);
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

    % set(get(colorbar),'Location','southoutside');\
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
    
    angles(d) = acos(x);  
    
    
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







