%close all;
figure;
set(gcf, 'position', [100 100 300 round(5/3 * 300)]);
subplot(2, 1, 1);

hold on; 
plot(reshape(ones(numel(MIchans), 1) * (1:8), [], 1), ...
    [adaRight(MIchans); adaNE(MIchans); adaTop(MIchans); adaNW(MIchans); ...
     adaLeft(MIchans); adaSW(MIchans); adaBottom(MIchans); adaSE(MIchans)], ...
     '.', 'color', [.6 .6 .7])
 
plot(1:8, mean([adaRight(MIchans), adaNE(MIchans), adaTop(MIchans), adaNW(MIchans), ...
     adaLeft(MIchans), adaSW(MIchans), adaBottom(MIchans), adaSE(MIchans)]), ...
     'lineWidth', 1, 'color', [0 0 0])
xlim([.5 8.5])
set(gca, 'xtick', 1:8);
ylabel('amplitude (au)')

subplot(2, 1, 2);
hold on;
plot(reshape(ones(numel(MIchans), 1) * (1:8), [], 1), ...
    ([abatRight(MIchans); abatNE(MIchans); abatTop(MIchans); abatNW(MIchans); ...
     abatLeft(MIchans); abatSW(MIchans); abatBottom(MIchans); abatSE(MIchans)]-1000)./1000, ...
     '.', 'color', [.6 .6 .7])
 
plot(1:8, (mean([abatRight(MIchans), abatNE(MIchans), abatTop(MIchans), abatNW(MIchans), ...
     abatLeft(MIchans), abatSW(MIchans), abatBottom(MIchans), abatSE(MIchans)])-1000)./1000, 'lineWidth', 1, 'color', [0 0 0])
xlim([.5 8.5])

x = reshape(ones(numel(MIchans), 1) * (1:8), [], 1);
y = ([abatRight(MIchans); abatNE(MIchans); abatTop(MIchans); abatNW(MIchans); ...
     abatLeft(MIchans); abatSW(MIchans); abatBottom(MIchans); abatSE(MIchans)]-1000)./1000;
 
myFun = @(x, z) z(1)*cos((x-1)*2*pi/8 + z(2)) + z(3);
[z, r] = lsqcurvefit(@(z, x) myFun(x, z), ones(3,1), x, y);

plot(1:.01:8, myFun(1:.01:8, z), 'lineWidth', 1, 'color', rgbStr2Mat('6b7ce0'))

disp(['R2 = ', num2str(1 - (r ./ sum((y- mean(y)).^2)))]);
set(gca, 'xtick', 1:8);
xlabel('movement direction');
ylabel('beta attenuation time');

%% 
close all;

hold on;
plot(reshape(ones(numel(MIchans), 1) * (1:8), [], 1), ...
    [amRight(MIchans); amNE(MIchans); amTop(MIchans); amNW(MIchans); ...
     amLeft(MIchans); amSW(MIchans); amBottom(MIchans); amSE(MIchans)], ...
     '.', 'color', [.6 .6 .7])
 
plot(1:8, mean([amRight(MIchans), amNE(MIchans), amTop(MIchans), amNW(MIchans), ...
     amLeft(MIchans), amSW(MIchans), amBottom(MIchans), amSE(MIchans)]), 'lineWidth', 1.5, 'color', [0 0 0])
%%

clc
disp('modulation depth test')
lm = LinearModel.fit( ordinal(reshape(ones(numel(MIchans), 1) * (1:8), [], 1)), ...
    [amRight(MIchans); amNE(MIchans); amTop(MIchans); amNW(MIchans); ...
     amLeft(MIchans); amSW(MIchans); amBottom(MIchans); amSE(MIchans)]);
 
disp(anova(lm))

disp(' ');
disp(' ');
disp('trough amplitude test');
lm = LinearModel.fit( ordinal(reshape(ones(numel(MIchans), 1) * (1:8), [], 1)), ...
    [aaRight(MIchans); aaNE(MIchans); aaTop(MIchans); aaNW(MIchans); ...
     aaLeft(MIchans); aaSW(MIchans); aaBottom(MIchans); aaSE(MIchans)]);
 
disp(anova(lm))

disp(' ');
disp(' ');
disp('amplitude at BAT test');
lm = LinearModel.fit( ordinal(reshape(ones(numel(MIchans), 1) * (1:8), [], 1)), ...
    [adaRight(MIchans); adaNE(MIchans); adaTop(MIchans); adaNW(MIchans); ...
     adaLeft(MIchans); adaSW(MIchans); adaBottom(MIchans); adaSE(MIchans)]);
 
disp(anova(lm))

disp(' ');
disp(' ');
disp('BAT test')

lm = LinearModel.fit( ordinal(reshape(ones(numel(MIchans), 1) * (1:8), [], 1)), ...
    [abatRight(MIchans); abatNE(MIchans); abatTop(MIchans); abatNW(MIchans); ...
     abatLeft(MIchans); abatSW(MIchans); abatBottom(MIchans); abatSE(MIchans)]);
 
disp(anova(lm));
%%
clc
amplitudeAtMO = zeros(96, 8);
for i = 1:8
    amplitudeAtMO(:,i) = squeeze(mean(abs(H(1001, i==beh(:,8), :)), 2));
end;

lm = LinearModel.fit(ordinal(reshape(ones(numel(MIchans), 1) * (1:8), [], 1)), ...
    reshape(amplitudeAtMO(MIchans,:), [], 1));
disp(anova(lm))


plot(ordinal(reshape(ones(numel(MIchans), 1) * (1:8), [], 1)), ...
    reshape(amplitudeAtMO(MIchans,:), [], 1), '.');


%%
clc
disp('amplitude at BAT k-w test')
[p, atab] = kruskalwallis(...
    [adaRight(MIchans); adaNE(MIchans); adaTop(MIchans); adaNW(MIchans); ...
     adaLeft(MIchans); adaSW(MIchans); adaBottom(MIchans); adaSE(MIchans)], ...
     ordinal(reshape(ones(numel(MIchans), 1) * (1:8), [], 1)))
