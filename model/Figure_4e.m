% clean workspace
clearvars
close all
clc

% path
cd('C:\Users\arnau\Dropbox\Groove\5b_Ed_model')

% model data
load('Mvpa_chanPs_Edall.mat')
fq2 = 29;

% Plot 
fig = figure; set(gcf,'color','w'); colormap([.5 .5 .5; 1 1 1]);
fig.Position = fig.Position./sqrt(sqrt(4));

cmap  = [.3 .3 .8 ; .5 .5 .5; .8 .3 .3];
hold on

Ip = 100;
xn = linspace(1, length(fth), length(fth) * Ip);

h = [];
high_stat = [0.09 0.1 0.08];
pval = 0.005;

leg = ["Layer 1" "Layer 2" "Layer 3"];
l = {};
for i0 = 1:3

    %stats
    p = [];
    for z0 = 1:length(fth)
        [~, p] = ttest(squeeze(OUT(:,i0,z0,fq2)),0, pval, 'both');
        p_tot(z0) = p;
    end

    [h] = mafdr(p_tot);%, 'BHFDR', 'true');
    z = NaN(1, length(fth)); z(h<pval) = 0;

    x = squeeze(OUT(:,i0,:,fq2));
    yn = spline(fth, x, xn);
    
    l{i0} = plot(xn, mean(yn), 'Linewidth', 2, 'Color', cmap(i0,:), 'Displayname',leg(i0));
    boundedline(xn, mean(yn), std(yn)/sqrt(size(yn,1)),'cmap',cmap(i0,:), 'transparency', 'alpha', .25);
        
    plot(fth, z-high_stat(i0), 'Color', cmap(i0,:), 'Linewidth', 3);

end


% contrast layer 1 vs layer 3
% stats
p = [];
pval = 0.05;

for z0 = 1:length(fth)
    [~, p] = ttest(squeeze(OUT(:,1,z0,fq2)),squeeze(OUT(:,3,z0,fq2)), pval, 'both');
    p_tot(z0) = p;
end

[h] = mafdr(p_tot);%, 'BHFDR', 'true');
z_con = NaN(1, length(fth)); z_con(h<pval & z == 0 ) = 0;

d = plot(fth, z_con-0.07, '.', 'Color', 'k', 'Linewidth', 3,  'DisplayName','Difference');
plot(fth, z_con-0.07,       'Color', 'k', 'Linewidth', 3)

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
legend([l{1} l{2} l{3} d], 'EdgeColor', 'w')

% legend([l{1} l{2} l{3}],'FontSize', 10, 'FontName', 'Arial');
% legend('boxoff');

helper_ILL(gcf, 1.29, 6.43);
print('C:\Users\arnau\Desktop\Groove_final\Figures\Fig4e','-painters', '-dpdf');


% saveas(gcf, 'C:\Users\arnau\Desktop\Groove_final\Figures\Fig4e.png'); 
