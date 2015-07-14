function all_data = build_filtered_data(filename,TOTAL_NUM_ELEC,NUM_OBS, PEAK_FRQ, SAMPLE_FRQ, BANDWIDTH)
    all_data = build_data(filename, TOTAL_NUM_ELEC, NUM_OBS);
    %We need to find a way to index into this 3D matrix to only filter
    %non-NaN values.
    %I think this line below takes a long time. 
%     all_data(~arrayfun(@(x) all(isnan(x)), all_data)) = abs(hilbert(filterData(all_data(~arrayfun(@(x) all(isnan(x)),all_data)), PEAK_FRQ, BANDWIDTH, SAMPLE_FRQ)));
end

function all_data = build_data(filename,TOTAL_NUM_ELEC,NUM_OBS, PEAK_FRQ, BANDWIDTH, SAMPLE_FRQ)
    f = whos('-file',filename,'-regexp','^M1lfp_MB_elec');
    NUM_TRIALS = f(1).size(2);
    all_data = NaN(NUM_OBS, NUM_TRIALS, TOTAL_NUM_ELEC);
    
    j = 1;
    for i = 1:3
        electrode_name = strcat('M1lfp_MB_elec',num2str(i,'%03d'));
        %if electrode exists
        if strcmp(electrode_name,f(j).name)
            load(filename, electrode_name);
            electrode_struct = arrayfun(@(x) x.times(:), eval(electrode_name), 'un', 0);
            all_data(:,~cellfun(@(x) all(isnan(x)), electrode_struct),i) = abs(hilbert(filterData(cell2mat(electrode_struct(~cellfun(@(x) all(isnan(x)), electrode_struct))),PEAK_FRQ,BANDWIDTH,SAMPLE_FRQ)));
            clear(electrode_name)
            j = j + 1
        end
    end    
end   