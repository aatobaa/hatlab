function [powerSpectra, f] = estimatePowerSpectrum(fileStr, varargin)
%ESTIMATEPOWERSPECTRUM using multi-taper estimators and chunking
%   [powerSpectra, f] = estimatePowerSpectrum(fileStr) loads an NSx file
%   specified by the string fileStr and computes a power spectrum for each
%   channel in the NSx file.  This function assumes that you have both Blackrock
%   tools and Chronux currently in your matlab path.  To verify that Blackrock
%   tools are in your path, you can type "which openNSx" on the command line.
%   If Blackrock tools are not in the path, it will return "'openNSx' not
%   found."  Similarly, you can use the command "which mtspectrumc" to verify
%   that Chronux is currently in the path.  
%
%   [powerSpectra, f] = estimatePowerSpectrum(fileStr, 'NAME', val, ...)
%   specifies one or more of the following name/value pairs:
%       Parameter       Value
%       'numFolds'      For memory efficiency, each time series is divided into
%                       numFolds chunks and then the power spectra for all
%                       chunks are averaged.  The default value is 1000.
%       'chans'         A vector specifying for which channels to estimate power
%                       spectra.  The default is 1:96.
%       'timeSamples'   A vector specifying which time samples to include in
%                       power spectrum estimates.  The default is
%                       1:size(LFP.Data, 2)
%       'params'        The params structure as defined by the function
%                       mtspectrumc in Chronux.
%
%   See also MTSPECTRUMC

% Created by Matt Best (mattbest@uchicago.edu) on 7/6/15
warning('Use openNSx 5.6.1.0 or newer');
try
    LFP = openNSx('report','read',fileStr);
catch
    error(['error loading ', fileStr]);
end;

p = inputParser;

numFolds = 1000;

chans = 1:96;
params.tapers = [5 9];
params.pad = 0;
try
    params.Fs = LFP.MetaTags.SamplingFreq;
catch
    params.Fs = 2000;
    warning('assuming a sample rate of 2 kHz')
end;
params.fpass = [.1 100];
params.trialave = 1;


if iscell(LFP.Data)
    addParameter(p, 'timeSamples', 1:size(LFP.Data{end}, 2));
else
    addParameter(p, 'timeSamples', 1:size(LFP.Data{end}, 2));
end;
addParameter(p, 'numFolds', numFolds, @(x) assert(isnumeric(x) && isscalar(x) && (x > 0), 'numFolds must be a positive, numeric scalar'));
addParameter(p, 'params', params);
addParameter(p, 'chans', chans, @(x) assert(isnumeric(x) && all(x > 0), 'chans must be a vector of positive numbers'));


parse(p, varargin{:});

numFolds = p.Results.numFolds;
params = p.Results.params;
chans = p.Results.chans;
timeSamples = p.Results.timeSamples;
e = chans(1);


if iscell(LFP.Data)
    foldSize = floor(numel(timeSamples) / numFolds);
    numTimePts = numel(timeSamples);
    timeSamples = timeSamples(1:(foldSize * numFolds));
    numTruncPts = numTimePts - numel(timeSamples);
    if numTruncPts
        warning(['The last ' num2str(numTruncPts), ' time points were truncated (', ...
            num2str(100 * (numTruncPts / numTimePts), '%2.1f'),'% of the data cut)'])
    end;
    data = double(LFP.Data{end}(e,timeSamples))';
else
    foldSize = floor(size(LFP.Data, 2) / numFolds);
    data = double(LFP.Data(e,1:(numFolds*foldSize)))';
end
data = reshape(data, [], numFolds);


tic
[S,f] = mtspectrumc(data,params);
toc

powerSpectra = nan(numel(S), numel(chans));
powerSpectra(:,1) = S;
for e = 2:numel(chans)
    data = double(LFP.Data{end}(chans(e),1:(numFolds*foldSize)))';
    data = reshape(data, [], numFolds);
    [S,~] = mtspectrumc(data,params);
    powerSpectra(:, e) = S;
    toc
end;
disp('finished computing powerSpectra')






