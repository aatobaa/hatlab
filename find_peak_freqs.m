datasets = ['PowerSpecM1TM_20111010.mat'; 'PowerSpecM1TM_20111014.mat'; 'PowerSpecM1TM_20111017.mat'; 'PowerSpecM1TM_20111019.mat'; 'PowerSpecM1TM_20111021.mat'; 'PowerSpecM1TM_20111025.mat'];



for d = 1:6
    dataset = datasets(d,:);
    clearvars -except dataset datasets
    load(dataset)
    clearvars -except powerspectra freqs dataset datasets
    
    min_freq = 10;
    max_freq = 45;

    freq = freqs(1,:);
    peak_freq = [];

    for chan = 1:96
        YoukePowerSpectra = powerspectra(:,chan);

        freq_range = freq(1,freq > min_freq & freq < max_freq);

        ps_range = YoukePowerSpectra(freq > min_freq & freq < max_freq,:);

        local_max = max(ps_range);

            %Find the minimum before the local maximum
        %     indx_min = 1;
        %     while ps_range(indx_min,i) > ps_range(indx_min+1,i)
        %         indx_min = indx_min + 1;
        %     end

        %We know beta is always greater than 16 hz empirically.
        endpnt = find(freq_range > 20,1);
        min_val = min(ps_range(1:endpnt));
        indx_min = find(ps_range == min_val,1);
        max_val = max(ps_range(indx_min:end));
        indx_max = find(ps_range == max_val,1);

    %     [troughs, trough_indx] = findpeaks(-ps_range(:,i),'minpeakdistance', 4, )


    %     [peaks,peaks_indx] = findpeaks(ps_range(indx_min:end,i),'minpeakdistance',10);
    %     peaks
        pf = ps_range(indx_max);
        pf_indx = find(ps_range(indx_min:end) == pf,1) + indx_min;
        pf = freq_range(pf_indx);
    %     pf
        if sum(powerspectra(:,chan)) == 0
            pf = 0;
        end
        peak_freq = [peak_freq; pf];

        %Use function 'findpeaks', remove CHRONUX from path! rmpath()
        % SET XLIMITS on plot, rather than indexing
        % peak_freq
    end
    save(strcat('peak_freq_',dataset),'peak_freq');
end


%% PLOT: Note, the day with a higher beta is pre-train day
plot(freq_range, 10*log10(ps_range),'LineWidth',2)

%% TestPlot
for i = 1:96
    plot(freq, 10*log10(powerspectra(:,i)))
    peak_frequency = peak_freq(i);
    line([peak_frequency peak_frequency],get(gca,'ylim'));
    pause
end

