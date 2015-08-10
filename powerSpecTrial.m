function [powerspectra, freqs] = powerSpecTrial(matfilename)

% matfilename
load(matfilename);

powerspectra = [];
freqs = [];

%%
if exist('M1lfp_MB_elec001','var')
    num_trials = size(M1lfp_MB_elec001,2);
elseif exist('M1lfp_MB_elec002','var')
    num_trials = size(M1lfp_MB_elec002,2);
elseif exist('M1lfp_MB_elec003','var')
    num_trials = size(M1lfp_MB_elec003,2);
elseif exist('M1lfp_MB_elec004','var')
    num_trials = size(M1lfp_MB_elec004,2);
end
num_electrodes = 96;
all_electrodes = zeros(4001,num_trials,num_electrodes);




for j = 1:num_electrodes
    trial_data = zeros(4001,num_trials);
    if j < 10
        if exist(strcat('M1lfp_MB_elec00',num2str(j)),'var') %Can also say Z = whos; and do Z.name; 
            electrode = eval(['M1lfp_MB_elec00',num2str(j)]);
        else
            electrode = 'null';
        end
    else
        if exist(strcat('M1lfp_MB_elec0',num2str(j)),'var')
            electrode = eval(['M1lfp_MB_elec0',num2str(j)]);
        else
            electrode = 'null';
        end
    end
    %deprecated%
    %num_trials = size(electrode,2);
    %if ~isstr(electrode)
    %    num_trials = size(electrode, 2);
    %else
    %    num_trials = 0;
    %end;
    if isstruct(electrode)
        for i = 1:num_trials
            trial = getfield(electrode,{i},'times')';
            if isnan(trial)
                trial_data(:,i) = 0;
            else
%                 trial_data(:,i) = abs(hilbert(filterData(trial, 18, 3, 1000)));   
                trial_data(:,i) = trial;
            end
        end
        all_electrodes(:,:,j) = trial_data;
    end 
end


params.tapers = [5 9];
params.pad = 0;
params.fpass = [.1 100];
params.Fs = 1000;
params.trialave = 1;

for i = 1:96
    data = all_electrodes(:,:,i);
    [S, f] = mtspectrumc(data, params);
    powerspectra = [powerspectra S];
    freqs = [freqs; f];
end

save(strcat('PowerSpec',matfilename));


