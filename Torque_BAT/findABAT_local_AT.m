function [aveBetaAttenuationTimes, attMag, attAmp ampAtBAT] = findABAT(H, floorThresh, makeFigure, fs, method, lambda)
% [aveBetaAttenuationTimes, attMag, attAmp, ampAtBAT] = findABAT(H, floorThresh, makeFigure, fs, method, lambda)
if nargin < 6
    lambda = 500;
end;
if nargin < 5
    method = 1;
end;
if nargin < 4
    fs = 2000;
end;
if nargin < 3
    makeFigure = 0;
end;
if nargin < 2
    floorThresh = 1.2;
end;

H = abs(H);
numChannels = size(H, 3);
numTrials = size(H,2);
aveBetaAttenuationTimes = zeros(numChannels, 1);
attMag = zeros(numChannels, 1);
attAmp = zeros(numChannels, 1);
ampAtBAT   = nan(numChannels, 1);

if method == 1
    for j = 1:numChannels
        meanEnvelope = mean(H(:,:,j), 2);
        [betaFloor, idx] = min(meanEnvelope(round(.5*fs):round(1.5*fs)));
        aveBetaAttenuationTimes(j) = find(meanEnvelope(round(.05*fs):idx+round(fs/2)+1) < floorThresh * betaFloor, 1)+round(0.05*fs)+1;
        
        if makeFigure
            clf
            plot((-fs:fs)./fs, meanEnvelope, 'color', 'b', 'lineWidth', 2);
            set(gca, 'fontSize', 14);
            line(([aveBetaAttenuationTimes(j) aveBetaAttenuationTimes(j)] - fs) ./ fs, get(gca, 'ylim'), 'color', 'k');
            line(([idx + round(fs/2)+1, idx + round(fs/2)+1] - fs) ./ fs, get(gca, 'ylim'), 'color', 'r')
            line([-1, 1], [floorThresh * betaFloor, floorThresh * betaFloor], 'color', 'r');
            %xlabel('time relative to movement onset (s)');
            %ylabel('LFP amplitude (au)');
            title(num2str(j))
            %xlabel(num2str(goodTrials(j)))
            pause;
        end;
    end;
elseif method == 2
    for j = 1:numChannels
        meanEnvelope = nanmean(H(:,:,j), 2);
        [betaFloor, fidx] = min(meanEnvelope(round(.5*fs):round(1.5*fs)));
        [betaCeil,  cidx] = max(meanEnvelope(round(.5*fs):round(1.5*fs)));
        tabat = find(meanEnvelope(cidx+fs/2:fidx+fs/2) < floorThresh * (betaCeil - betaFloor) + betaFloor, 1) + cidx+fs/2 + 1;
        if ~isempty(tabat)
            aveBetaAttenuationTimes(j) = tabat;
            ampAtBAT(j) = meanEnvelope(tabat);
        else
            aveBetaAttenuationTimes(j) = nan;
            disp(['no attenuation time found for channel ', num2str(j)]);
        end;
        attMag(j) = meanEnvelope(cidx + round(fs/2)+1) - meanEnvelope(fidx + round(fs/2)+1);
        attAmp(j) = meanEnvelope(fidx + round(fs/2)+1);
        if makeFigure
            clf
            plot((-fs:fs)./fs, meanEnvelope, 'color', 'b', 'lineWidth', 2);
            set(gca, 'fontSize', 14);
            line(([aveBetaAttenuationTimes(j) aveBetaAttenuationTimes(j)] - fs) ./ fs, get(gca, 'ylim'), 'color', 'k');
            line(([fidx + round(fs/2)+1, fidx + round(fs/2)+1] - fs) ./ fs, get(gca, 'ylim'), 'color', 'r')
            line(([cidx + round(fs/2)+1, cidx + round(fs/2)+1] - fs) ./ fs, get(gca, 'ylim'), 'color', 'r')
            %line([-1, 1], [(1+floorThresh) * betaFloor, (1+floorThresh) * betaFloor], 'color', 'r');
            line([-1, 1], [floorThresh * (betaCeil - betaFloor) + betaFloor, floorThresh * (betaCeil - betaFloor) + betaFloor], 'color', 'r');
            %xlabel('time relative to movement onset (s)');
            %ylabel('LFP amplitude (au)');
            title(num2str(j))
            %xlabel(num2str(goodTrials(j)))
            pause;
        end;
    end;
elseif method == 3
    for j = 1:numChannels
        meanEnvelope = tvd(mean(H(:,:,j), 2), lambda);
        [betaFloor, fidx] = min(meanEnvelope(round(.5*fs):round(1.5*fs)));
        [betaCeil,  cidx] = max(meanEnvelope(round(.5*fs):round(1.5*fs)));
        tabat = find(meanEnvelope(cidx+fs/2:fidx+fs/2) < floorThresh * (betaCeil - betaFloor) + betaFloor, 1) + cidx+fs/2 + 1;
        if ~isempty(tabat)
            aveBetaAttenuationTimes(j) = tabat;
            ampAtBAT(j) = meanEnvelope(tabat);
        else
            aveBetaAttenuationTimes(j) = nan;
            disp(['no attenuation time found for channel ', num2str(j)]);
        end;
        attMag(j) = meanEnvelope(cidx + round(fs/2)+1) - meanEnvelope(fidx + round(fs/2)+1);
        attAmp(j) = meanEnvelope(fidx + round(fs/2)+1);
        if makeFigure
            clf
            plot((-fs:fs)./fs, meanEnvelope, 'color', 'b', 'lineWidth', 2);
            set(gca, 'fontSize', 14);
            line(([aveBetaAttenuationTimes(j) aveBetaAttenuationTimes(j)] - fs) ./ fs, get(gca, 'ylim'), 'color', 'k');
            line(([fidx + round(fs/2)+1, fidx + round(fs/2)+1] - fs) ./ fs, get(gca, 'ylim'), 'color', 'r')
            line(([cidx + round(fs/2)+1, cidx + round(fs/2)+1] - fs) ./ fs, get(gca, 'ylim'), 'color', 'r')
            %line([-1, 1], [(1+floorThresh) * betaFloor, (1+floorThresh) * betaFloor], 'color', 'r');
            line([-1, 1], [floorThresh * (betaCeil - betaFloor) + betaFloor, floorThresh * (betaCeil - betaFloor) + betaFloor], 'color', 'r');
            %xlabel('time relative to movement onset (s)');
            %ylabel('LFP amplitude (au)');
            title(num2str(j))
            %xlabel(num2str(goodTrials(j)))
            pause;
        end;
    end;
elseif method == 4
    for j = 1:numChannels
        meanEnvelope = mean(H(:,:,j), 2);
        [betaFloor, fidx] = min(meanEnvelope(round(.5*fs):round(1.05*fs)));
        [betaCeil,  cidx] = max(meanEnvelope(round(.5*fs):round(1.05*fs)));
        tabat = find(meanEnvelope(cidx+fs/2:fidx+fs/2) < floorThresh * (betaCeil - betaFloor) + betaFloor, 1) + cidx+fs/2 + 1;
        if ~isempty(tabat)
            aveBetaAttenuationTimes(j) = tabat;
            ampAtBAT(j) = meanEnvelope(tabat);
        else
            aveBetaAttenuationTimes(j) = nan;
            disp(['no attenuation time found for channel ', num2str(j)]);
        end;
        attMag(j) = meanEnvelope(cidx + round(fs/2)+1) - meanEnvelope(fidx + round(fs/2)+1);
        attAmp(j) = meanEnvelope(fidx + round(fs/2)+1);
        if makeFigure
            clf
            plot((-fs:fs)./fs, meanEnvelope, 'color', 'b', 'lineWidth', 2);
            set(gca, 'fontSize', 14);
            line(([aveBetaAttenuationTimes(j) aveBetaAttenuationTimes(j)] - fs) ./ fs, get(gca, 'ylim'), 'color', 'k');
            line(([fidx + round(fs/2)+1, fidx + round(fs/2)+1] - fs) ./ fs, get(gca, 'ylim'), 'color', 'r')
            line(([cidx + round(fs/2)+1, cidx + round(fs/2)+1] - fs) ./ fs, get(gca, 'ylim'), 'color', 'r')
            %line([-1, 1], [(1+floorThresh) * betaFloor, (1+floorThresh) * betaFloor], 'color', 'r');
            line([-1, 1], [floorThresh * (betaCeil - betaFloor) + betaFloor, floorThresh * (betaCeil - betaFloor) + betaFloor], 'color', 'r');
            %xlabel('time relative to movement onset (s)');
            %ylabel('LFP amplitude (au)');
            title(num2str(j))
            %xlabel(num2str(goodTrials(j)))
            pause;
        end;
    end;
elseif method == 5
    for j = 1:numChannels
        numTrials = size(H, 2);
        meanEnvelope = mean([zeros(500,numTrials); zscore(abs(H(500:1500,:,j))); zeros(500,numTrials)], 2);
        
        [betaFloor, fidx] = min(meanEnvelope(round(.5*fs):round(1.5*fs)));
        [betaCeil,  cidx] = max(meanEnvelope(round(.5*fs):round(1.5*fs)));
        tabat = find(meanEnvelope(cidx+fs/2:fidx+fs/2) < floorThresh * (betaCeil - betaFloor) + betaFloor, 1) + cidx+fs/2 + 1;
        if ~isempty(tabat)
            aveBetaAttenuationTimes(j) = tabat;
            ampAtBAT(j) = meanEnvelope(tabat);
        else
            aveBetaAttenuationTimes(j) = nan;
            disp(['no attenuation time found for channel ', num2str(j)]);
        end;
        attMag(j) = meanEnvelope(cidx + round(fs/2)+1) - meanEnvelope(fidx + round(fs/2)+1);
        attAmp(j) = meanEnvelope(fidx + round(fs/2)+1);
        if makeFigure
            clf
            plot((-fs:fs)./fs, meanEnvelope, 'color', 'b', 'lineWidth', 2);
            set(gca, 'fontSize', 14);
            line(([aveBetaAttenuationTimes(j) aveBetaAttenuationTimes(j)] - fs) ./ fs, get(gca, 'ylim'), 'color', 'k');
            line(([fidx + round(fs/2)+1, fidx + round(fs/2)+1] - fs) ./ fs, get(gca, 'ylim'), 'color', 'r')
            line(([cidx + round(fs/2)+1, cidx + round(fs/2)+1] - fs) ./ fs, get(gca, 'ylim'), 'color', 'r')
            %line([-1, 1], [(1+floorThresh) * betaFloor, (1+floorThresh) * betaFloor], 'color', 'r');
            line([-1, 1], [floorThresh * (betaCeil - betaFloor) + betaFloor, floorThresh * (betaCeil - betaFloor) + betaFloor], 'color', 'r');
            %xlabel('time relative to movement onset (s)');
            %ylabel('LFP amplitude (au)');
            title(num2str(j))
            %xlabel(num2str(goodTrials(j)))
            pause;
        end;
    end;
elseif method == 6
    'in method 6'
%METHOD 6 averages across electrodes to attempt to get a BAT for EACH
%TRIAl.
aveBetaAttenuationTimes = zeros(numTrials, 1);
attMag = zeros(numTrials, 1);
attAmp = zeros(numTrials, 1);
ampAtBAT   = nan(numTrials, 1);

    for j = 1:numTrials
        meanEnvelope = nanmean(squeeze(H(:,j,:)), 2);
        [betaFloor, fidx] = min(meanEnvelope(round(.5*fs):round(1.5*fs)));
        [betaCeil,  cidx] = max(meanEnvelope(round(.5*fs):round(1.5*fs)));
        tabat = find(meanEnvelope(cidx+fs/2:fidx+fs/2) < floorThresh * (betaCeil - betaFloor) + betaFloor, 1) + cidx+fs/2 + 1;
        if ~isempty(tabat)
            aveBetaAttenuationTimes(j) = tabat;
            ampAtBAT(j) = meanEnvelope(tabat);
        else
            aveBetaAttenuationTimes(j) = nan;
            disp(['no attenuation time found for trial ', num2str(j)]);
        end;
        attMag(j) = meanEnvelope(cidx + round(fs/2)+1) - meanEnvelope(fidx + round(fs/2)+1);
        attAmp(j) = meanEnvelope(fidx + round(fs/2)+1);
        if makeFigure
            clf
            plot((-fs:fs)./fs, meanEnvelope, 'color', 'b', 'lineWidth', 2);
            set(gca, 'fontSize', 14);
            line(([aveBetaAttenuationTimes(j) aveBetaAttenuationTimes(j)] - fs) ./ fs, get(gca, 'ylim'), 'color', 'k');
            line(([fidx + round(fs/2)+1, fidx + round(fs/2)+1] - fs) ./ fs, get(gca, 'ylim'), 'color', 'r')
            line(([cidx + round(fs/2)+1, cidx + round(fs/2)+1] - fs) ./ fs, get(gca, 'ylim'), 'color', 'r')
            %line([-1, 1], [(1+floorThresh) * betaFloor, (1+floorThresh) * betaFloor], 'color', 'r');
            line([-1, 1], [floorThresh * (betaCeil - betaFloor) + betaFloor, floorThresh * (betaCeil - betaFloor) + betaFloor], 'color', 'r');
            %xlabel('time relative to movement onset (s)');
            %ylabel('LFP amplitude (au)');
            title(num2str(j))
            %xlabel(num2str(goodTrials(j)))
            pause;
        end;
    end;
end;