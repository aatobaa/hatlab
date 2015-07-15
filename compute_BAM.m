%REMOVE is an array of all electrodes that should be removed from analysis.
%
function BAM = compute_BAM(data, NUM_ELEC, START_TIME, END_TIME, REMOVE)
    BAM = zeros(96,1);
%     BAMT = zeros(96,1);
    block = round((END_TIME - START_TIME) / 3);
%     block
    for i = 1:NUM_ELEC
        trial_average = nanmean(data(START_TIME:END_TIME,:,i),2);
%         tatest = nanmean(data(START_TIME:END_TIME,1:30,i),2);
%         size(trial_average)

        LB = [min(trial_average(1:block)) min(trial_average(block*2:block*3)) -1, 1];
        UB = [max(trial_average(1:block)) max(trial_average(block*2:block*3)) -.01, 3000];
        
%         LBT = [min(tatest(1:block)) min(tatest(block*2:block*3)) -1, 1];
%         UBT = [max(tatest(1:block)) max(tatest(block*2:block*3)) -.01, 2000];
        
        
%         [LB LBT]
%         x = START_TIME:END_TIME;
        x = START_TIME:END_TIME;
        clear cHat
        if(~isnan(sum(trial_average)))
            cHat = fitSigmoid2(x, trial_average, LB, UB);
            %The fourth argument of cHat is BAM
        end
%         if(~isnan(sum(tatest)))
%             cHatTest = fitSigmoid2(x,tatest, LBT, UBT);
%         end
        
        if exist('cHat','var') && ~any(i == REMOVE) %tossing out electrodes
            BAM(i) = cHat(4); %I don't think I need this here+ START_TIME - 1;
%             cHat(4)
        end
%         if exist('cHatTest','var') && ~any(i == REMOVE)
%             BAMT(i) = cHatTest(4);
%         end
    end
end
