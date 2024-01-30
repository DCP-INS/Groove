%% Figure 2e&f & S1a&b : loading files

clearvars
close all
clc

cd('C:/Users/arnau/Desktop/Groove_final')

% frequency of interest (model)
fq = 2; % Hz (for both acoustic and Ed fft data) - (Hz: 2|6.01|14.03|8.1623)

% labels
SUJETS  = {'002','003','004','005','006','007','008','009','010','012'  ...
           '014','015','016','017','018','019','020','021','022','023','024','025'...
           '026','027','028','029','030','031','032'}; % (001&013out) 11 removed
                      
Cond    = {'Low', 'Medium', 'High'};
Network = {'2 Hz layer 1 (a.u.)', '2 Hz layer 2 (a.u.)', '2 Hz layer 3 (a.u.)', 'Layer 2-1'}; % Layer 3
Regress = {'2 Hz Acoustic (dB)', 'Degree of syncopation', 'Groove ratings (MEG)', 'Groove ratings (online)'};
nstim = 12;

% load acoustic
load('./Data/stimuli/accenv_Ed.mat');
[~,b] = min(abs(out.fq-fq));
out.D = out.D*1e2; %!
IN{1} = out.D(:,:,b);
clear out b

% load syncopation
X = xlsread('./Data/stimuli/Index de syncopation.xlsx');
X = X(:,2); % use poly-phonic index (2)
X = reshape(X, 3, length(X)/3);
sync = [X(3,:); X(1:2,:)]; % R-G-S
IN{2} = sync;
clear sync X

% load behavior
Rpf = nan(length(SUJETS),3, 4, 12); Rid = Rpf;
for isuj = 1:length(SUJETS)
   
    % load data
    for irun = 1:4
        raw = importdata(['./Data/behavior/MEG/logfiles/'...
                            'groove_' SUJETS{isuj} '/'...
                            SUJETS{isuj} '-groove' num2str(irun) '.log']);
                        
        stim = raw.textdata(:,2); %!
            stim = str2num(cell2mat(stim(3:end)));
            stim(stim == 100) = [];
        resp = raw.data(:,1); %!
            resp(isnan(resp)) = [];       
            
        % !organise fn condition + sort fn stim id!
        for idata = 1:length(Cond)
            x = stim(round(stim/1e2)==idata);
            [x,var] = sort(x);
            Rid(isuj,idata,irun,1:length(x)) = x;            
            x = resp(round(stim/1e2)==idata);
            Rpf(isuj,idata,irun,1:length(x)) = x(var);          
        end
        clear stim resp idata x raw var
    end
end

% avg repetitions & subjects
Rpf = squeeze(nanmean(nanmean(Rpf,3),1));
IN{3} = Rpf;
clear Rpf Rid x


% load online experiment data
stimuli = {'chrisG','chrisR','chrisS','chris2G','chris2R','chris2S', ...
            'cosmicG','cosmicS','cosmicR','dannoG','dannoS','dannoR','danoG', ...
            'danoS','danoR','danuG','danuS','danuR','dinoG','dinoS','dinoR', ...
            'funkadelicG','funkadelicS','funkadelicR','loudanoG','loudanoS', ...
            'loudanoR','metersG','metersS','metersR','monkG','monkS','monkR', ...
            'motoG','motoS','motoR','pleaseG','pleaseS','pleaseR','rockyG', ...
            'rockyS','rockyR','telegirlG','telegirlS','telegirlR'};
        
stimrank = nan(1,length(stimuli));
for i0 = 1:length(stimuli)
    if strcmp(stimuli{i0}(end), 'R'), stimrank(i0) = 1;
    elseif strcmp(stimuli{i0}(end), 'G'), stimrank(i0) = 2;
    elseif strcmp(stimuli{i0}(end), 'S'), stimrank(i0) = 3;
    end
end

% cd([origDir '2 onlinexp/'])
files = dir('./Data/behavior/onlinexp/*/Musical_Groove*.txt');
for i0 = 1:length(files)
    % load & sort data
    x = readtable([files(i0).folder '/' files(i0).name ]);
    x = table2array(x(:,1:3)); % data: pf|rt|id
    [~,b] = sort(x(:,3)); % sort fn stim id
    x = x(b,:);    
    % initialize ouput
    if (~exist('out', 'var'))
        out = nan(length(files), 3, size(x,1)/3, size(x,2));
    end    
    % divide fn stim type (R,G,S)
    for j0 = 1:3, out(i0,j0,:,:) = x(stimrank==j0,:); end
end
IN{4} = squeeze(nanmean(out(:,:,[1:7 9:11 13:14],1), 1)); % ! remove 3 melodies of non-interest
clear x files b stimrank stimuli


% load Ed model data
load('./Data/model/model_3l_v2_n29_fs60_dup0_tmin2_230106_final.mat')

Ed = squeeze(mean(out.D,4)); % avg nrun
[~,b] = min(abs(out.fq-fq));
Ed = Ed(:,:,:,b);

%% Figure S1b   : display figure

% show each melody (avg subjects)
figure; set(gcf,'color','w');
cmap = [0 0 0 ; .4 .4 .4; .8 .8 .8];

for i0 = 1:3 % layer
    x = squeeze(Ed(i0,:,:));
    % x = bsxfun(@plus, bsxfun(@minus, x, mean(x,2)), mean(x(:))); % repeated-measure sem
    
    n = length(x);    
    subplot(3,3,3*i0-1)    
    hold on
    
    for z0 = 1:3 % condition
        boxplot(x,'boxstyle','outline', 'colors','k','OutlierSize',4,'Symbol','','Widths',0.2)
        scatter(ones(n,1).*(z0+(rand(n,1)-0.5)/10), x(:,z0), ...
            'MarkerFaceColor',cmap(z0,:),'MarkerEdgeColor','k',...
            'MarkerFaceAlpha',.9,'MarkerEdgeAlpha',.8)
    end

    if i0 ==1;      ylabel( '2 Hz layer 1 (a.u.)' ); ylim( [0, 2.5] ); yticks([0 1 2])
    elseif i0 == 2; ylabel( '2 Hz layer 2 (a.u.)' ); ylim( [0, 2.5] ); yticks([0 1 2])
    elseif i0 == 3; ylabel( '2 Hz layer 3 (a.u.)' ); ylim( [0, .2] ); yticks([0 .1 .2]) %!
    end

    set(gca, 'FontSize', 8, 'FontName', 'Calibri')
    xlim( [.5, 3.5] )
    set(gca,'XTickLabel', [] ) % Cond
%         xlabel( 'Condition' )
    set(gca,'Layer','top','Box','off','TickLength',[.01 .01])
end

helper_ILL(gcf, 1.1, 9, false);
print('./Figures/FigS1b','-painters', '-dpdf');

%% Figure S1a   : display figure

% load Ed model data
load('./Data/model/model_3l_v2_n29_fs60_dup0_tmin2_230106_final.mat')
Ed = squeeze(mean(out.D,4)); % avg nrun
[~,b] = min(abs(out.fq-fq));
Ed = Ed(:,:,:,b);
Ed(isnan(Ed)) = 0; % remove NaN if any (but see above)

x_reg = [0 12; -2 17; 1 5; 1 7];
y_reg = [0 2.5; 0 2.5; 0 .2; -.1 .15]; % 3+1 layers


hfig = figure; set(gcf,'color','w');
col = [repmat([0 0 0],nstim,1); repmat([.4 .4 .4],nstim,1); repmat([.8 .8 .8],nstim,1)]; % b|g|r

% Ed Layers
f0 = 0;
for y0 = 1:size(Ed,1)
    x1 = reshape(Ed(y0,:,:), [],1);

    % Regressors
    for x0 = 1:numel(IN) 
        x2 = reshape(permute(IN{x0}, [2 1]), [],1);

        % Pearson correlation (df = n-2)
        [r,~,p] = Pearson(x1, x2, 0);
        r = r.^2;
        fprintf('   %-15s r= %5.2f, p= %4.4f\n', [Network{y0} ':'], r,p);
        
        % figure
        f0 = f0+1; set(0, 'CurrentFigure', hfig); 
        subplot(size(Ed,1), numel(IN) , f0); hold on    
        
        scatter(x2, x1, 25, col, 'filled', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
        h1 = lsline; set(h1, 'Color','k', 'Linewidth', 1)
        set(gca, 'FontSize', 8, 'FontName', 'Calibri')
        set(gca,'Layer','top','Box','off','TickLength',[.01 .01])
        axis('square')
        
        xlim(x_reg(x0,:))
        ylim(y_reg(y0,:))               
        
        if x0 == 1; ylabel(Network{y0}); end
        if y0 == size(Ed,1); xlabel(Regress{x0}); end
 
        text( min(xlim)+max(xlim)*.5, min(ylim)+max(ylim)*.1, sprintf('r^2 = %5.2f', r), ...
            'FontSize', 8, 'FontName', 'Calibri') % ['r^2 = ' num2str(round(r, 2))]        
    end
end
helper_ILL(gcf, 1.1, 9, false);
print('./Figures/FigS1a','-painters', '-dpdf');

%% Figure 2e
% figure article 1/2
col = [repmat([0 0 0],nstim,1); repmat([.4 .4 .4],nstim,1); repmat([.8 .8 .8],nstim,1)]; % b|g|r

figure; set(gcf,'color','w'); 
% subplot(2,2,1)

x2 = reshape(permute(IN{2}, [2 1]), [],1); % syncope
x1 = reshape(Ed(1,:,:), [],1);  % Layer 1

    % Pearson correlation
    [r,~,p] = Pearson(x1, x2, 0);
    r = r*r;

    % figure     
    scatter(x2, x1, 20, col,'filled', 'MarkerFaceAlpha',.7,'MarkerEdgeAlpha',.8)    
    h1 = lsline; set(h1, 'Color','k', 'Linewidth', 1.2)

    set(gca, 'FontSize', 11, 'FontName', 'Calibri')
    set(gca,'Layer','top','Box','off','TickLength',[.01 .01])

    xlabel('Degree of syncopation')
        xlim( [-2 17] )
        xticks( 0:5:15 )    

    ylabel('2 Hz layer 1 (a.u.)')
        a = ylim;
        ylim([0 3])
        yticks([1 2])
    
    text(min(xlim)+max(xlim)*.1,min(ylim)+max(ylim)*.1, ['r^2 = ' num2str(round(r, 2))], ...
            'FontSize', 10, 'FontName', 'Calibri')
        
helper_ILL(gcf, 1.2, 4,false);
print('./Figures/Fig2e','-painters', '-dpdf');

% saveas(gcf, [origDir '5_bhv&acc/figures/fig2e_Edmodel_corr_sync.tif'])

%% Figure 2f

col = [repmat([0 0 0],nstim,1); repmat([.4 .4 .4],nstim,1); repmat([.8 .8 .8],nstim,1)]; % b|g|r

% figure article 2/2
figure; set(gcf,'color','w'); 
% subplot(2,2,1)

x2 = reshape(permute(IN{4}, [2 1]), [],1); % groove
x1 = reshape(Ed(3,:,:), [],1);  % Layer 3

    % Pearson correlation
    [r,~,p] = Pearson(x1, x2, 0);
    r = r*r;

    % figure     
    scatter(x2, x1, 20, col,'filled', 'MarkerFaceAlpha',.7,'MarkerEdgeAlpha',.8)    
    h1 = lsline; set(h1, 'Color','k', 'Linewidth', 1.2)

    set(gca, 'FontSize', 11, 'FontName', 'Calibri')
    set(gca,'Layer','top','Box','off','TickLength',[.01 .01])

    xlabel('Groove ratings')
        xlim( [1 7] )
        xticks( 1:2:7 )    

    ylabel('2 Hz layer 3 (a.u.)')
        a = ylim;
        ylim([0 .17])
        yticks([0 .1 .2])
    
    text(min(xlim)+max(xlim)*.1,min(ylim)+max(ylim)*.1, ['r^2 = ' num2str(round(r, 2))], ...
            'FontSize', 10, 'FontName', 'Calibri')

helper_ILL(gcf, 1.2, 4,false);
print('./Figures/Fig2f','-painters', '-dpdf');
% saveas(gcf, [origDir '5_bhv&acc/figures/fig2f_Edmodel_corr_groove.tif'])
