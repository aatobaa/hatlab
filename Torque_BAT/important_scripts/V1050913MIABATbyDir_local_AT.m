
load('P-V1050913_M1Ab_PMvAc001.mat', 'conditions', 'events')
%load('Y:\Aaron\Premotor Corex Tuning\S-Files\Velma\S-V1050913_M1Ab_PMvAc001.mat')
V1050913LFP = openNSx('report', 'read', 'V1050913_M1Ab_PMvAc001.ns2');
disp('done');
%%
% good trials are defined as trials in which a reward was delivered
numTrials = size(conditions(1).epochs, 1);
mo = zeros(numTrials, 1);
me = zeros(numTrials, 1);

% This will spit out an error because k goes to 1104
k = 1;
for i = 1:numTrials
    while conditions(1).epochs(i,2) > conditions(5).epochs(k,1)
        mo(i) = conditions(5).epochs(k,1);
        me(i) = conditions(5).epochs(k,2);
        k = k + 1;
    end;
end;


%%
X = V1050913LFP.Data(1:64, :)';
chan2rc = makechan2rc('v', 'mi');
H = zeros(2001, numTrials, 64);
for j = 1:64
    Y = zeros(2001, numTrials);
    for i = 1:numTrials
        Y(:,i) = X(round(mo(i)*1000)-1000:round(mo(i)*1000)+1000,j);
    end;
    H(:, :, j) = hilbert(filterData(Y, 21, 3, 1000));
end;
goodChans = setdiff(1:64, [40 44]);
%%
movementDirections = nan(numTrials, 1);


[~, idx, ~] = intersect(conditions(1).epochs(:,1), ...
                        conditions(10).epochs(:,1));
movementDirections(idx) = 1;
                    
[~, idx, ~] = intersect(conditions(1).epochs(:,1), ...
                        conditions(11).epochs(:,1));
movementDirections(idx) = 2;

[~, idx, ~] = intersect(conditions(1).epochs(:,1), ...
                        conditions(12).epochs(:,1));
movementDirections(idx) = 3;                    

[~, idx, ~] = intersect(conditions(1).epochs(:,1), ...
                        conditions(13).epochs(:,1));
movementDirections(idx) = 4;

[~, idx, ~] = intersect(conditions(1).epochs(:,1), ...
                        conditions(14).epochs(:,1));
movementDirections(idx) = 5;

[~, idx, ~] = intersect(conditions(1).epochs(:,1), ...
                        conditions(15).epochs(:,1));
movementDirections(idx) = 6;

[~, idx, ~] = intersect(conditions(1).epochs(:,1), ...
                        conditions(16).epochs(:,1));
movementDirections(idx) = 7;

[~, idx, ~] = intersect(conditions(1).epochs(:,1), ...
                        conditions(17).epochs(:,1));
movementDirections(idx) = 8;
                    


fastTrials = (me - mo) < median(me - mo);
%% 8 directions
[abatLeft, amLeft, aaLeft, adaLeft] = findABAT(...
    H(:, fastTrials & (movementDirections == 5), :), ...
    .15, 0, 1000, 2);

[abatRight, amRight, aaRight, adaRight] = findABAT(...
    H(:, fastTrials & (movementDirections == 1), :), ...
    .15, 0, 1000, 2);

[abatTop, amTop, aaTop, adaTop]    = findABAT(...
    H(:, fastTrials & (movementDirections == 3), :), ...
    .15, 0, 1000, 2);

[abatBottom, amBottom, aaBottom, adaBottom] = findABAT(...
    H(:, fastTrials & (movementDirections == 7), :), ...
    .15, 0, 1000, 2);

[abatNE, amNE, aaNE, adaNE] = findABAT(...
    H(:, fastTrials & (movementDirections == 2), :), ...
    .15, 0, 1000, 2);

[abatNW, amNW, aaNW, adaNW] = findABAT(...
    H(:, fastTrials & (movementDirections == 4), :), ...
    .15, 0, 1000, 2);

[abatSW, amSW, aaSW, adaSW] = findABAT(...
    H(:, fastTrials & (movementDirections == 6), :), ...
    .15, 0, 1000, 2);

[abatSE, amSE, aaSE, adaSE] = findABAT(...
    H(:, fastTrials & (movementDirections == 8), :), ...
    .15, 0, 1000, 2);



%% blocked directions 
[abatLeft, amLeft, aaLeft, adaLeft] = findABAT(...
    H(:, fastTrials & (movementDirections == 4 | movementDirections == 5 | movementDirections == 6), :), ...
    .15, 0, 1000, 2);

[abatRight, amRight, aaRight, adaRight] = findABAT(...
    H(:, fastTrials & (movementDirections == 2 | movementDirections == 1 | movementDirections == 8), :), ...
    .15, 0, 1000, 2);

[abatTop, amTop, aaTop, adaTop]    = findABAT(...
    H(:, fastTrials & (movementDirections == 4 | movementDirections == 3 | movementDirections == 2), :), ...
    .15, 0, 1000, 2);

[abatBottom, amBottom, aaBottom, adaBottom] = findABAT(...
    H(:, fastTrials & (movementDirections == 8 | movementDirections == 7 | movementDirections == 6), :), ...
    .15, 0, 1000, 2);

[abatNE, amNE, aaNE, adaNE] = findABAT(...
    H(:, fastTrials & (movementDirections == 3 | movementDirections == 2 | movementDirections == 1), :), ...
    .15, 0, 1000, 2);

[abatNW, amNW, aaNW, adaNW] = findABAT(...
    H(:, fastTrials & (movementDirections == 3 | movementDirections == 4 | movementDirections == 5), :), ...
    .15, 0, 1000, 2);

[abatSW, amSW, aaSW, adaSW] = findABAT(...
    H(:, fastTrials & (movementDirections == 5 | movementDirections == 7 | movementDirections == 6), :), ...
    .15, 0, 1000, 2);

[abatSE, amSE, aaSE, adaSE] = findABAT(...
    H(:, fastTrials & (movementDirections == 8 | movementDirections == 7 | movementDirections == 1), :), ...
    .15, 0, 1000, 2);
