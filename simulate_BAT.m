function bootstrapped_times = simulate_BAT(all_electrodes, NUM_ELEC,TRIAL_RANGE, NUM_SIMS,BIN,START_TIME, END_TIME,REMOVE)
    %Simulate NUM_SIMS experiments by randomly selecting BIN trials in
    %TRIAL_RANGE for each experiment. Yields a matrix of NUM_SIMS beta attenuation 
    %times for each electrode
    
    %TRIAL_RANGE: Array of trials to use in simulation
    %NUM_SIMS: Number of simulations to run (per bin)
    %BIN: Number of trials to include in each average

    bootstrapped_times = zeros(NUM_ELEC,NUM_SIMS);
    if BIN > max(TRIAL_RANGE) - min(TRIAL_RANGE)
        error('BIN must be smaller than the number of trials')
    end
    for j = 1:NUM_SIMS
        i = round(rand(BIN,1)*(max(TRIAL_RANGE)-min(TRIAL_RANGE)) + min(TRIAL_RANGE));
        BAM = compute_BAM(all_electrodes(:,i,:), NUM_ELEC, START_TIME, END_TIME, REMOVE);
        bootstrapped_times(:,j) = BAM;    
        j
    end        
