clearvars
close all
clc

cd('C:\Users\arnau\Desktop\Groove_final\Data')

% Load decoding results
load('Ps_az_decoding_goodchan_goodtrial_GrooveS_nfold10_alpha2.mat')

% remove bad subject
s(10,:,:) = []; % remove 011 if not removed during computation
size(s) % must be 29 participants

for z0 = 1:length(fth)
    
    pval = [];
    % regressors
    pval = 0.005;
    [h, p] = ttest(squeeze(s(:,1,z0)),0, pval, 'both');
    p_tot1(z0) = p;

    [h, p] = ttest(squeeze(s(:,2,z0)),0, pval, 'both');
    p_tot2(z0) = p; 
    
    [h, p] = ttest(squeeze(s(:,3,z0)),0, pval, 'both');
    p_tot3(z0) = p; 
    
    % contrast syncop vs groove
    pval = 0.05;
    [h, p] = ttest(squeeze(s(:,3,z0)),squeeze(s(:,2,z0)), pval, 'both');
    p_tot4(z0) = p; 
    
    clear h p
end

% regressors
pval = 0.005;
[h] = mafdr(p_tot1); %, 'BHFDR', 'true');
z1 = NaN(1, length(fth)); z1(h<pval) = 0; % fth(z1 == 0)

[h] = mafdr(p_tot2); %, 'BHFDR', 'true');
z2 = NaN(1, length(fth)); z2(h<pval) = 0;

[h] = mafdr(p_tot3); %, 'BHFDR', 'true');
z3 = NaN(1, length(fth)); z3(h<pval) = 0;

% contrast
pval = 0.05;
[h] = mafdr(p_tot4); %, 'BHFDR', 'true');
z4 = NaN(1, length(fth)); z4(h<pval & z3==0 ) = 0;

fig = figure; set(gcf,'color','w'); colormap([.5 .5 .5; 1 1 1]);
fig.Position = fig.Position./sqrt(2);

hold on
Ip = 100;
xn = linspace(1, length(fth), length(fth) * Ip);

y = squeeze(s(:, 2,:)); % syncope
yn = spline(fth, y, xn);
    plot(xn, mean(yn), 'Color',[.3 .3 .8], 'Linewidth',2)
    boundedline(xn, mean(yn), std(yn)/sqrt(size(yn,1)), ...
        'cmap',[.3 .3 .8], 'transparency', 'alpha', .25);
        
sync = plot(fth, z2-0.09, 'Color', [.3 .3 .8],'Linewidth', 3, 'DisplayName','Degree of syncopation');

y = squeeze(s(:, 3,:)); % Groove
yn = spline(fth, y, xn);
    plot(xn, mean(yn),'Color',[.8 .3 .3], 'Linewidth',2)
    boundedline(xn, mean(yn), std(yn)/sqrt(size(yn,1)), ...
        'cmap', [.8 .3 .3], 'transparency', 'alpha', .25); 
        
g = plot(fth, z3-0.08, 'Color', [.8 .3 .3], 'Linewidth', 3, 'DisplayName','Groove ratings');

d = plot(fth, z4-0.07, '.', 'Color', 'k', 'Linewidth', 3,  'DisplayName','Difference');
plot(fth, z4-0.07,       'Color', 'k', 'Linewidth', 3)

yline(0,':');
nbTicks = 6 ;
xlabel('Neural frequency (Hz)')
ylabel('Coding precision (z)')
xticks([2 5 10 20 50 100])
xticklabels([2 5 10 20 50 100])
yticks([-.1 0 .1 .2 .3  .4])

ylim([-.1 .35])
set(gca,'XScale', 'log')
set(gca, 'FontSize', 10, 'FontName', 'Arial')
set(gca,'Layer','top','Box','off','TickLength',[.001 .001])
legend([sync g d], 'EdgeColor', 'w')

helper_ILL(gcf, 1.29, 6.43);
print('C:\Users\arnau\Desktop\Groove_final\Figures\Fig4d','-painters', '-dpdf');

