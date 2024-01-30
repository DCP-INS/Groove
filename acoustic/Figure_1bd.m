%% Figure 1bd : loading files
clearvars
close all
clc

cd('C:/Users/arnau/Desktop/Groove_final')
addpath(genpath('./tools'))

melodies   = {'chris', 'chris2', 'cosmic', 'Danno', 'Dano', 'Danu', ...
              'Dino', 'LouDano', 'meters', 'monk', 'please', 'rocky'};
conditions = {'_reg_NM_2Hz', '', '_hsync_M'};
Fs  = 32;     % [default: 48]
dur = [0 16]; % in sec

% load acoustic (Ed)
load('./Data/stimuli/accenv_Ed.mat'); % accenv_norm1.mat

%% Figure 1b

x =[]; x_mean =[]; x_std =[]; 
Ip = 100;
xn = linspace(out.fq(1), length(out.fq), length(out.fq) * Ip);

figure; set(gcf,'color','w'); colormap([.5 .5 .5; 1 1 1])
colors = [0 0 0 ; .4 .4 .4; .8 .8 .8];
set(gca, 'ColorOrder', colors);% legend(Cond)    

% acoustic env mod spectrum
x = (out.D * 1e2); % rescale
x_mean = spline(out.fq, squeeze(mean(x,2)) , xn);
x_std = spline(out.fq, squeeze(std(x,[],2) /sqrt(size(x,2))) , xn);

h = plot(xn, x_mean,'Linewidth',1.2);

ylim( [0, 12] )
yticks([0 6 12])
xticks([1 2 4 6 8])
xlim([0 8.5])  
xlabel('Frequency (Hz)')
ylabel('Amplitude (a.u.)')

set(gca,'Layer','top','Box','off','TickLength',[.001 .001])
set(gca, 'FontSize', 11, 'FontName', 'Calibri')
set(h, {'color'}, num2cell(colors,2));

helper_ILL(gcf, 1.6, 3,false);
print('./Figures/Fig1b','-painters', '-dpdf');

%% Figure 1d

% load Syncopation index data
X = xlsread('./Data/stimuli/Index de syncopation.xlsx');
X = X(:,2); % use poly-phonic index (2)
X = reshape(X, 3, length(X)/3);
sync = X([3 1 2],:); % R-G-S (Low|Med|High)

% load acoustic (new: Ed) - extract 2 Hz
fq = 2;
load('./Data/stimuli/accenv_Ed.mat'); % accenv_norm1.mat

[~,b] = min(abs(out.fq-fq));
acc = 1e2*out.D(:,:,b);
acc = mag2db(acc);% 2hz acoustic

clear X out b

cmap = [0 0 0 ; .4 .4 .4; .8 .8 .8];

figure; set(gcf,'color','w'); 
hold on

for i0 = 1:3    
    scatter(sync(i0,:), acc(i0,:), 20, ...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',cmap(i0,:),...
    'MarkerFaceAlpha',.9,...
    'MarkerEdgeAlpha',0,...
    'LineWidth',.05)
end

% compute adjusted r-value of fits
x1 = reshape(sync', [], 1 ); 
x2 = reshape(acc', [], 1);
[~,s] = sort(x1);

rsq_adj = [];

for j0 = 1:4 %!
    n = numel(x1);
    p = polyfit(x1(s), x2(s), j0); % set order of fit
    f = polyval(p, x1(s));
    SSE1 = sum((x2(s) -f).^2);
    SSE2 = (n-1) *var(x2(s));

    rsq_adj(j0) = 1 -SSE1/SSE2 *(n-1)/(n-length(p));
    clear n p f SSE1 SSE2 rsq
end    

fprintf('1st order fit: r2_adj= %5.2f \n2nd order fit: r2_adj= %5.2f\n3rt order fit: r2_adj= %5.2f \n', ...
         rsq_adj(1), rsq_adj(2),rsq_adj(3));   
         
% Pearson correlation
[r,~,p] = Pearson(x1, x2, 0);
r = r.^2;
fprintf('Pearson: r2 = %1.2f p = %1.3f df = %2.0f \n', r, p, numel(x1)-2);   

% plot fit(s)
n = 1; % n order fit
p = polyfit(x1(s), x2(s), n);
    tmp_vec = min(x1(s)):.1:max(x1(s));
    f = polyval(p, tmp_vec);        
    
hold on;
plot(tmp_vec, f, 'k', 'Linewidth', 1.2)     

ylabel('2 Hz acoustic (dB)')
    ylim( [5 25] )
    yticks( [5 15 25]  )
    
xlabel('Degree of syncopation')
    xlim( [-2 17] )
    xticks( 0:5:15 )
    
text(min(xlim)+max(xlim)*.1,min(ylim)+max(ylim)*.07,['r^2 = ' num2str(round(r, 2))], ...
            'FontSize', 10, 'FontName', 'Calibri');
        
set(gca,'Layer','top','Box','off','TickLength',[.001 .001])
set(gca, 'FontSize', 11, 'FontName', 'Calibri')             


helper_ILL(gcf, 1.1, 3, false);
print('./Figures/Fig1d','-painters', '-dpdf');

clear p f tmp_vec
clear x1 x2 s rsq_adj
