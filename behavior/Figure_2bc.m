%% Figure 2b
clearvars
close all
clc

cd('C:/Users/arnau/Desktop/Groove_final')
addpath(genpath('./tools'))

% parameters
% warning off
Cond    = {'Low', 'Medium', 'High'};
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
nstim = 12;

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
clear x files b

% ! remove 3 melodies of non-interest
stimuli = {'chrisG','chrisR','chrisS','chris2G','chris2R','chris2S', ...
    'cosmicG','cosmicS','cosmicR','dannoG','dannoS','dannoR','danoG', ...
    'danoS','danoR','danuG','danuS','danuR','dinoG','dinoS','dinoR', ...
    'loudanoG','loudanoS', 'loudanoR','metersG','metersS','metersR', ...
    'monkG','monkS','monkR', 'pleaseG','pleaseS','pleaseR','rockyG', ...
    'rockyS','rockyR'};

clear stimrank
out = out(:,:,[1:7 9:11 13:14],:);
pf = out(:,:,:,1);
IN{1} = pf; % performance accuracy
clear pf rt

% load Syncopation index data
X = xlsread('./Data/stimuli/Index de syncopation.xlsx');
X = X(:,2); % use poly-phonic index (2)
X = reshape(X, 3, length(X)/3);
sync = X([3 1 2],:); % R-G-S (Low|Med|High)
IN{2} = sync; % syncopation
clear sync X

% load acoustic (Ed): 2 hz - *100
fq_acc = 2;
t = load('./Data/stimuli/accenv_Ed.mat'); % accenv_norm1.mat
[~,b] = min(abs(t.out.fq-fq_acc));
IN{3} = 1e2*t.out.D(:,:,b); 
fq_acc = t.out.fq;
clear t

% figure    
colors = [0 0 0 ; .4 .4 .4; .8 .8 .8];
col = [repmat([0 0 .5],nstim,1); repmat([0 .6 0],nstim,1); repmat([.8 0 0],nstim,1)]; % b|g|r
set(gcf,'color','w'); colormap([.5 .5 .5; 1 1 1])

for i0 = 1:3
    hold on
    h = errorbar(IN{2}(i0,:), squeeze(mean(IN{1}(:,i0,:))), squeeze(std(IN{1}(:,i0,:)))/sqrt(size(IN{1},1)), ...
    '.', 'color', colors(i0,:), 'MarkerSize',9.6, 'Linewidth', 1.2,'Capsize', 1);

    % Set transparency (undocumented)
    alpha = 0.65;        
    set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha]); 
end

% compute adjusted r-value of fits
x1 = IN{2}; 
x2 = squeeze(mean(IN{1},1)); % avg subjects
x1 = reshape(x1', [], 1); x2 = reshape(x2', [], 1);
[~,s] = sort(x1);

rsq_adj = [];
for j0 = 1:2 %!
    n = numel(x1);
    p = polyfit(x1(s), x2(s), j0); % set order of fit
    f = polyval(p, x1(s));
    SSE1 = sum((x2(s) -f).^2);
    SSE2 = (n-1) *var(x2(s));
    rsq_adj(j0) = 1 -SSE1/SSE2 *(n-1)/(n-length(p));
    clear n p f SSE1 SSE2 rsq
end
fprintf('1st order fit: r2_adj= %5.2f \n2nd order fit: r2_adj= %5.2f \n', rsq_adj(1),rsq_adj(2));
        
% plot fit(s)
p = polyfit(x1(s), x2(s), 2); % 2nd order fit
    tmp_vec = min(x1(s)):.1:max(x1(s));
    f = polyval(p, tmp_vec);                
hold on;

plot(tmp_vec, f, 'k', 'Linewidth', 1.2)  

set(gca, 'FontSize', 11, 'FontName', 'Calibri')

ylabel('Groove ratings')
ylim( [1 7] ) %!
yticks( [1 3 5 7] )

xlabel('Degree of syncopation')
xlim( [-2 17] )
xticks( 0:5:15 )

set(gca,'Layer','top','Box','off','TickLength',[.001 .001])
clear p f tmp_vec

text(min(xlim)+max(xlim)*.15, min(ylim)+max(ylim)*.07,['r^2 = ' num2str(round(rsq_adj(2), 2))], ...
    'FontSize', 10, 'FontName', 'Calibri');

helper_ILL(gcf, 1.2, 4,false);
print('./Figures/Fig2b','-painters', '-dpdf');

% clear x1 x2 s rsq_adj

%% Figure 2c

clearvars
close all
clc

% origDir = '/Users/Bn/Dropbox/Groove/';
% origDir = 'C:/Users/zalta/Dropbox/Groove/';
cd('C:/Users/arnau/Desktop/Groove_final')
addpath(genpath('./tools'))

% addpath(genpath([origDir '5_bhv&acc/functions/']))
SUJETS  = {'002','003','004','005','006','007','008','009','010','012'  ...
           '014','015','016','017','018','019','020','021','022','023','024','025'...
           '026','027','028','029','030','031','032'}; % (001&013out) without 11

Cond    = {'Low', 'Medium', 'High'};
nsuj    = length(SUJETS);
nstim   = 12;

% load: stimulus index (Rid), Grooviness (Rpf), Syncopation (polyphonic), acoustic
IN = [];
Rpf = nan(nsuj,3, 4, 12); Rid = Rpf;
for isuj = 1:nsuj
%     subject = ['groove_' SUJETS{isuj}];
%     cd([origDir '3b logfiles/' subject])
%     
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

% avg repetitions
Rpf = squeeze(nanmean(Rpf,3));
IN{1} = Rpf;
clear Rpf Rid

% load Syncopation index data (Vuust v. is not working)
% cd([origDir '1 stimuli/'])
X = xlsread('./Data/stimuli/Index de syncopation.xlsx');
X = X(:,2); % use poly-phonic index (2)
X = reshape(X, 3, length(X)/3);
sync = X([3 1 2],:); % R-G-S (Low|Med|High)
IN{2} = sync;
clear sync X

% load acoustic (Ed): 2 hz - *100
fq_acc = 2;
t = load('./Data/stimuli/accenv_Ed.mat');
% t = load([origDir '5_bhv&acc/data_tmp/accenv_Ed.mat']); % accenv_norm1.mat
[~,b] = min(abs(t.out.fq-fq_acc));
IN{3} = 1e2*t.out.D(:,:,b); 
fq_acc = t.out.fq;
clear t

% or load acoustic (morlet, bn)
%cd([origDir '5_bhv&acc/data_tmp'])
%t = load([origDir '5_bhv&acc/data_tmp/accenv_norm1.mat']);
%fq_acc = t.fq; % stim = t.stim;
%IN{3} = t.OUT;
%clear t

x1 = IN{2}; % 2|3: syncope|acoustic
x2 = squeeze(mean(IN{1},1)); % avg subjects
x1 = reshape(x1', [], 1); x2 = reshape(x2', [], 1);

[~,s] = sort(x1);
% colors = [.3 .3 .8 ; .3 .8 .3; .8 .3 .3];
cmap = [0 0 0 ; .4 .4 .4; .8 .8 .8];

% figure for article    
%     col = [repmat([0 0 .5],nstim,1); repmat([0 .6 0],nstim,1); repmat([.8 0 0],nstim,1)]; % b|g|r
%     % col = [repmat([0 0 0],nstim,1); repmat([0 0 0],nstim,1); repmat([0 0 0],nstim,1)]; % k
%     fig = figure; set(gcf,'color','w'); colormap([.5 .5 .5; 1 1 1])
%     fig.Position =  fig.Position./sqrt(2); 

% for poster 
figure; set(gcf,'color','w');
% subplot(2,2, 1)
    
%%%%% Moyenne au travers des sujets%%%%%
%     boxplot(x2(s),'boxstyle','outline', 'colors',[0 0 0]) %!
%     scatter(x1(s), x2(s), 60, col(s,:), 'filled',...
%         'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5) %!
        

%%%% Avec Error bars %%%%%%%
%  x = IN{1} - mean(IN{1},2) + mean(IN{1}(:)); % Normalisation Bizarre...
for i0 = 1:3
    hold on
    h = errorbar(IN{2}(i0,:), squeeze(mean(IN{1}(:,i0,:))), squeeze(std(IN{1}(:,i0,:)))/sqrt(size(IN{1},1)), ...
    '.', 'color', cmap(i0,:), 'MarkerSize',9.6, 'Linewidth', 1.2,'Capsize', 1);
    % Set transparency (undocumented)
    alpha = 0.65;   
    set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
end

                  
    % compute adjusted r-value of fits
    rsq_adj = [];
    for j0 = 1:2 %!
        n = numel(x1);
        p = polyfit(x1(s), x2(s), j0); % set order of fit
        f = polyval(p, x1(s));
        SSE1 = sum((x2(s) -f).^2);
        SSE2 = (n-1) *var(x2(s));
        % rsq = 1 - SSE1/SSE2;
        rsq_adj(j0) = 1 -SSE1/SSE2 *(n-1)/(n-length(p));
        clear n p f SSE1 SSE2 rsq
    end    
    fprintf('1st order fit: r2_adj= %5.2f \n2nd order fit: r2_adj= %5.2f \n', rsq_adj(1),rsq_adj(2));   
        
    % plot fit(s)
    % h1 = lsline; set(h1,'Color',[.5 .5 .5], 'Linewidth', 2) % 1st order fit
    p = polyfit(x1(s), x2(s), 2); % 2nd order fit
        tmp_vec = min(x1(s)):.1:max(x1(s));
        f = polyval(p, tmp_vec);                
    hold on; 
    plot(tmp_vec, f, 'k', 'Linewidth', 1.2) 
    
text(3,1.5,['r^2 = ' num2str(round(rsq_adj(2), 2))], 'FontSize', 10, 'FontName', 'Calibri');

    set(gca, 'FontSize', 11, 'FontName', 'Calibri')
    ylabel('Groove ratings')
        ylim( [1 5] ) %!
        yticks( [1 3 5] )
    xlabel('Degree of syncopation')
        xlim( [-2 17] )
        xticks( 0:5:15 )
    set(gca,'Layer','top','Box','off','TickLength',[.001 .001])
    clear p f tmp_vec

helper_ILL(gcf, 1.2, 4,false);
print('./Figures/Fig2c','-painters', '-dpdf');

% saveas(gcf, [origDir '5_bhv&acc/figures/fig2c_bhvb_MEG.tif'])
clear x1 x2 s rsq_adj
