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
%     d = 2;
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
    chan2rc = makechan2rc_mac('y','mio');
    
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
    [bootstrapped_BAM,bootstrapped_trials] = simulate_BAT(all_electrodes,NUM_ELEC, TRIAL_RANGE, NUM_SIMS, BIN, START_TIME, END_TIME,[82]);
    
    %% Convert zeros to NaN
    bootstrapped_BAM(~bootstrapped_BAM) = nan;
    
    %% BAT Bootstrap Plot Prep (I don't think this does anything here)
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
    indx_trials_with_high_r2 = []
    num_qual = 0;
    for i = 1:NUM_SIMS
        ds_b = dataset(bootstrapped_BAM(:,i), plotB(:,2),plotB(:,3), 'VarNames', {'BAT','X','Y'});
        simulated_model = LinearModel.fit(ds_b, 'BAT~X+Y');
        R2_b = simulated_model.Rsquared.Ordinary;
%         if R2_b > .25
%             R2_b
%             indx_trials_with_high_r2 = [indx_trials_with_high_r2; i];
%             num_qual = num_qual + 1
%         end
        coefs_b = simulated_model.Coefficients.Estimate;
        coefs_b_norm = norm([(coefs_b(3)/2),(coefs_b(2)/2)]);
        if ~exist('BAO_b_coef','var')
            BAO_b_coef = {coefs_b filename R2_b 'angle'};
        else
            BAO_b_coef = [BAO_b_coef; {coefs_b filename R2_b 'angle'}];
        end
    end
    
    %% Testing VS
%     for q = 1:size(indx_trials_with_high_r2,1)
% %         z
%         z = indx_trials_with_high_r2(q);
%         coefs_b = cell2mat(BAO_b_coef(z,1));
%         t_filename = char(BAO_b_coef(z,2));
%         coefs_b_norm = norm([(coefs_b(3)/2),(coefs_b(2)/2)]);
%         x = ((coefs_b(3)/2)/coefs_b_norm);
%         y = ((coefs_b(2)/2)/coefs_b_norm);
%         %index 3 is the R2 of the vector 
%         r2_b = cell2mat(BAO_b_coef(z,3));
% %         r2_b
% %         'hi'
% %         x
%         %index 4 is the angle of the vector in radians
%         BAO_b_coef(z,4) = {atan2(y,x)};    
%     end
%     vs = vectorStrength(cell2mat(BAO_b_coef(indx_trials_with_high_r2,4)));
%     wvs = weightedVectorStrength(cell2mat(BAO_b_coef(indx_trials_with_high_r2,4)),cell2mat(BAO_b_coef(indx_trials_with_high_r2,3)));
    
    
    %% Testing VS
%     other_indices = [];
% 
%     for q = 1:1000
% %         z
%         if any(q == [indx_trials_with_high_r2])
%             continue
%             'pass'
%         end
%         other_indices = [other_indices; q];
%         z = q;
%         coefs_b = cell2mat(BAO_b_coef(z,1));
%         t_filename = char(BAO_b_coef(z,2));
%         coefs_b_norm = norm([(coefs_b(3)/2),(coefs_b(2)/2)]);
%         x = ((coefs_b(3)/2)/coefs_b_norm);
%         y = ((coefs_b(2)/2)/coefs_b_norm);
%         %index 3 is the R2 of the vector 
%         r2_b = cell2mat(BAO_b_coef(z,3));
% %         r2_b
% %         'hi'
% %         x
%         %index 4 is the angle of the vector in radians
%         BAO_b_coef(z,4) = {atan2(y,x)};    
%     end
%     vs = vectorStrength(cell2mat(BAO_b_coef(other_indices,4)));
%     wvs = weightedVectorStrength(cell2mat(BAO_b_coef(other_indices,4)),cell2mat(BAO_b_coef(other_indices,3)));
%     
%     
%     
%     
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
        BAO_b_coef(i,4) = {atan2(y,x)};    
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
            BAO_bin_coef(i,4) = {atan2(y,x)};      
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
end
%%

save('youke_results', 'results')






