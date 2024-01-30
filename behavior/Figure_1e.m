%% Figure 1e : loading files
clearvars
close all
clc

cd('C:/Users/arnau/Desktop/Groove_final')
addpath(genpath('./tools'))

% analyse all subjects
num_std = 0; % 1 or 2 ?

SUJETS = {'S01','S02','S03','S04','S05', 'S06', 'S07', 'S08','S09', 'S10', 'S11', 'S12', 'S13','S14'};
stim_order = {'chris_' , 'chris2_', 'cosmic_', 'Danno_', 'Dano_', 'Danu_', 'Dino_', 'LouDano_', 'meters_', 'monk_', 'please_', 'rocky_'};

% load data: 1 session per participant
IN = [];
Dif = NaN(length(SUJETS),1);

for isuj = 1:length(SUJETS)

    ifile = dir(['./Data/behavior/tapping_exp/TA18_' SUJETS{isuj} '_exp_control*.mat']);
    load([ifile.name])
    
    IN{isuj,1}  = X;    
    
    nstim = numel(Stim)/3;
    cond = 2*ones(numel(Stim),1);
    cond(contains(X.stim, '_reg_')) = 1;
    cond(contains(X.stim, '_hsync_')) = 3;
    
    D = 1./diff(X.tapp,[],2); % compute instantaneous frequency (in Hz)
    
    % remove outlier durations
    if num_std ~= 0
        D(D > nanmean(D(:))+num_std*nanstd(D(:)) ) = NaN;
        D(D < nanmean(D(:))-num_std*nanstd(D(:)) ) = NaN;
    else
    end
 
    for i0 = 1:3
        IN{isuj}.freq_instant(:,i0,:) = D(cond==i0,:);        
    end 
    
    D = nanmean(D,2); % avg along trial
    IN{isuj}.out = [D(cond==1), D(cond==2), D(cond==3)];
    
    out(isuj,:,:) = IN{isuj}.out;
    out_instant(isuj,:,:,:) = IN{isuj}.freq_instant;
    
    clear X subject ifile D cond 
end  

% reshape
nstim = numel(Stim)/3;
out = squeeze(mean(out, 1));
out_instant = reshape(permute(out_instant, [1 2 4 3]), nstim*599*length(SUJETS), 3);

%% Figure 1e : display figure

color = [0 0 0 ; .4 .4 .4; .8 .8 .8];
cond = {'Low', 'Medium', 'High'};

bin = linspace(0, 13, 131);
Ip = 100;
xn = linspace(bin(2)-0.05, bin(end)-0.05, length(bin) * Ip);

set(gcf,'color','w'); colormap([.5 .5 .5; 1 1 1])

hold on
for i0 = 1:3    
    h = histogram(out_instant(:,i0),bin, 'FaceColor', color(i0,:), 'EdgeColor', color(i0,:));
    h.FaceAlpha = 0;
    h.EdgeAlpha = 0;
    
    yn = spline(h.BinEdges(2:end)-0.05, h.Values, xn);
    p(i0) = plot(xn, yn, 'Color',color(i0,:), 'Linewidth',1.2, 'DisplayName',cond{i0});
    
end

clear h yn

ylabel( 'Number of events' )
xlabel( 'Instantaneous frequency (Hz)' )
xticks([1 2 4 6 8])
yticks([0 600 1200])
yticklabels([0 600 1200])
xlim([0 8.5])
ylim([0 1200])

set(gca, 'FontSize', 11, 'FontName', 'Calibri')
set(gca,'Layer','top','Box','off','TickLength',[.001 .001])


helper_ILL(gcf, 1.6, 3,false);
print('./Figures/Fig1e','-painters', '-dpdf');