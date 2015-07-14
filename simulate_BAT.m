function bootstrapped_times = simulate_BAT(data, NUM_ELEC,TRIAL_RANGE, NUM_SIMS,BIN,REMOVE)
%TRIAL_RANGE is a column vector of trials [1 2 323 322 etc.] indicating the
%trials to sample from 
%BIN is the number of trials to use for each average.  
    if(length(TRIAL_RANGE) < BIN)
        error('Number of bins must be less than the number of trials to sample from');
    end
    
    bootstrapped_times = zeros(NUM_ELEC,NUM_SIMS);
    for k = 1:NUM_ELEC
        for j = 1:NUM_SIMS
            %Could also use randsample(1:num_trials, 50)
            %DO YOU WANT TO BE ABLE TO SPECIFY WHICH TRIALS YOU WANT, OR
            %ONLY SPECIFY THE RANGE FROM WHICH YOU DRAW THEM? i.e.
            %implement as 32 or 1:32?
            i = round(rand(BIN,1)*(NUM_TRIALS-1) + 1);
            trial_average = nanmean(data(:,i,k),2);
            clear min_ba max_ba cutoff tidx
            min_ba = min(trial_average(START_TIME:END_TIME));
            max_ba = max(trial_average(START_TIME:END_TIME));
            cutoff = min_ba + (max_ba - min_ba) * THRESHOLD; %cutoff is the voltage at which point we say beta attenuation occured at this time.
            tidx = find(trial_average(START_TIME:END_TIME) < cutoff, 1); %Find the first time the average amplitude is less than the cutoff voltage (because we expect beta to decrease)
            if ~isempty(tidx) && ~any(k == REMOVE) %tossing out electrode 82
                bootstrapped_times(k,j) = tidx + START_TIME - 1;
            else
                bootstrapped_times(k,j) = NaN;
            end        
        end
    end
end