function all_data = build_data(filename,TOTAL_NUM_ELEC,NUM_OBS, PEAK_FRQ, BANDWIDTH, SAMPLE_FRQ)
    f = whos('-file',filename,'-regexp','^M1lfp_MB_elec');
    NUM_TRIALS = f(1).size(2);
    all_data = NaN(NUM_OBS, NUM_TRIALS, TOTAL_NUM_ELEC);
    
    j = 1;
    for i = 1:TOTAL_NUM_ELEC
        electrode_name = strcat('M1lfp_MB_elec',num2str(i,'%03d'));
        %if electrode exists
        num_working_elec = size(f,1);
        if j <= num_working_elec && strcmp(electrode_name,f(j).name)
            load(filename, electrode_name);
            electrode_struct = arrayfun(@(x) x.times(:), eval(electrode_name), 'un', 0);
            all_data(:,~cellfun(@(x) all(isnan(x)), electrode_struct),i) = abs(hilbert(filterData(cell2mat(electrode_struct(~cellfun(@(x) all(isnan(x)), electrode_struct))),PEAK_FRQ,BANDWIDTH,SAMPLE_FRQ)));
            clear(electrode_name)
            j = j + 1, filename
        end
    end    
end   