clearvars
close all
clc

cd('C:/Users/arnau/Desktop/Groove_final')

% via PAC spatial profile
% only delta and beta left parietal
load('./Data/Ps_brut_spatial_profile_scout_nbin8_p3_zscored_0.mat')
load('./Data/scouts_Pac_bn.mat', 'scout', 'scout_Name')

LabelY = {'\delta' '\beta'};
Frqs = [1.3 1.5; 20 30];
limy = [9.4e-3 10.4e-3; 1.75e-3 1.9e-3];
ticky = [9.4e-3 9.9e-3 10.4e-3 ; 1.75e-3 1.82e-3 1.9e-3];
ireg = 3; % groove
iscout = 2; % left parietal
regs = {'Syncope' 'Groove'};

fig = figure;
fig.Position = fig.Position./sqrt(2);

left_color = [.5 .5 .5];
right_color = [0 0 0];

set(fig,'defaultAxesColorOrder',[left_color; right_color]);
set(gcf,'color','w');
hold on

for f0 = 1:size(Frqs,1)

    if f0 == 1     
        yyaxis left
    elseif f0 == 2
        yyaxis right
    end
    
    % define index of frequency of interest
    [~,ifreq(1)] = min(abs(fth - Frqs(f0,1)));
    [~,ifreq(2)] = min(abs(fth - Frqs(f0,2)));

    x = squeeze(mean(mean(xbin(:,iscout,ireg,ifreq(1):ifreq(2),:),1),4));
    y = squeeze(mean(ybin(:,iscout,ireg,ifreq(1):ifreq(2),:),4));
    y = bsxfun(@plus,bsxfun(@minus, y, mean(y,2)), mean(y(:)));

    plot(x, mean(y),  'Marker', '.', 'MarkerSize', 8);
    errorbar(x, mean(y), std(y)/sqrt(size(y,1)), ...
        '.', 'MarkerSize', 15, 'Linewidth', 1.5,'Capsize', 1)
    

    xlim([ 1, 5 ])
    xticks([1 3 5])
    ylim(limy(f0,:))
    yticks([])
    ylabel([ LabelY{f0} ' amplitude (a.u.)'])
end

xlabel('Groove')
set(gca, 'FontSize', 11, 'FontName', 'Calibri')
set(gca,'Layer','top','Box','off','TickLength',[.001 .001])
set(gca,'box','off')


helper_ILL(gcf, 1.1, 5, false);
print('./Figures/Fig5d','-painters', '-dpdf');
