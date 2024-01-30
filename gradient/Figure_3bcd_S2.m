%% Figure 3b&c create gradient & figures
clearvars
close all
clc

cd('C:\Users\arnau\Desktop\Groove_final')

% load PS for gradient 
load('.\Data\Ps_task\Ps_task_gradient_subj_3cond.mat', 'out', 'fth');
load('.\Data\Ps_rest\Ps_rest_gradient_subj_3cond.mat', 'out_rest');

nvox = size(out, 2);

% avg 3 conditions 
out = mean(out,3);
out_rest = mean(out_rest,3);

% % avg 30 subjects
% out = mean(out,1);
% out_rest = mean(out_rest,1);

% zscore to see the spectral gradient
zout = squeeze(zscore(out, [], 2));
zout_rest = squeeze(zscore(out_rest, [], 2));

% avg subjects (after zscore, so that each subject is normalised)
zout = squeeze(mean(zout,1));
zout_rest = squeeze(mean(zout_rest,1));

% estimate per subject and voxel the preferred frequency
TF = zeros(nvox, 1);
TF_rest = zeros(nvox, 1);
for ivox = 1:nvox
    fq = find(zout(ivox,:) == max(zout(ivox,:)));
    fq_rest = find(zout_rest(ivox,:) == max(zout_rest(ivox,:)));
    TF(ivox) = fth(fq);
    TF_rest(ivox) = fth(fq_rest);
end

% compare rest and data
dim = 5; % degree of fit
fmax = 45;

% load data and rest from bst
X = load('./Data/Ps_task/timefreq_fft_211015_1605.mat'); % source_Ps_az Frequency gradient complete data (histogram) in bst
X.TF = TF; % data or rest, loaded above

% X = load('./Data/Ps_rest/timefreq_fft_211015_1606.mat'); % source_Ps_rest_az Frequency gradient complete data (histogram)  in bst
% X.TF = TF_rest; 

fig = figure; set(gcf,'color','w');
fig.Position(3) =  250; fig.Position(4) =  250;


rsq = [];
X.GridLoc(X.TF>fmax,:) == 0;
X.TF(X.TF>fmax) = 0;

cubicFit = []; cubicCoef = [];
for i0 = 1:3
    IN = X.GridLoc(:,i0);
    [cubicCoef,stats,ctr] = polyfit(X.TF, IN, dim);
    cubicFit(i0, :) = polyval(cubicCoef, X.TF, [],ctr);
    
    n = length(X.GridLoc); % 1673
    p = polyfit(X.TF, IN, dim);
    f = polyval(p, X.TF);
    df = length(p);
    
    SSE1 = sum((IN-f).^2); 
    SSE2 = (n-1) * var(IN);  
    rsq(i0)  = 1 - SSE1/SSE2 * (n-1) /(n-df); % R-square adjusted
    clear n p f df SSE1 SSE2 IN
end

ax1 = axes;
scatter3(ax1, cubicFit(1, :), cubicFit(2, :), cubicFit(3, :), 35, X.TF, 'filled','MarkerEdgeColor', 'k');
axis equal

ax1.YLim = [-0.03 0.03];
ax1.XLim = [-0.07 0.07];
ax1.ZLim = [0.02 0.1];

ax1.XLabel.String = 'Y (mm)';
ax1.YLabel.String = 'X (mm)'; ax1.YLabel.Rotation = -35; ax1.YLabel.Position = [-0.1 -0.04 0.01];
ax1.ZLabel.String = 'Z (mm)';

view(-27,20)

ax1.GridLineStyle = ':';
ax1.YTick = [-0.03 0 0.03];
ax1.XTick = [-0.06 -.03 0 0.03 0.06];
ax1.ZTick = [0.03 0.06 0.09];
ax1.YTickLabel = {'-30' '0' '30'};
ax1.XTickLabel = {'-60' '-30' '0' '30' '60'};
ax1.ZTickLabel = {'30' '60' '90'};
ax1.TickLength = [0.002 0.002];
colormap(ax1,'parula');

set(gca, 'FontSize', 11, 'FontName', 'Calibri') % 10 11 or 12 ?

% helper_ILL(gcf, 1.2, 3.5, false);
% print('./Figures/Fig3c','-painters', '-dpdf');

%% Figure 3d1
clearvars
close all
clc

cd('C:\Users\arnau\Desktop\Groove_final')

% load PS for gradient 
load('.\Data\Ps_task\Ps_task_gradient_subj_3cond.mat', 'out', 'fth');
load('.\Data\Ps_rest\Ps_rest_gradient_subj_3cond.mat', 'out_rest');

nvox = 1673;

% avg 3 conditions - zscore to see the spectral gradient
zout = zscore(mean(out,3), [], 2);
zout_rest = zscore(mean(out_rest,3), [], 2);

% estimate per subject and voxel the preferred frequency
TF = zeros(size(out,1), nvox, 1);
TF_rest = zeros(size(out,1), nvox, 1);
for isuj = 1:size(out,1)
    for ivox = 1:nvox
        fq = find(zout(isuj, ivox,:,:) == max(zout(isuj, ivox,:,:)));
        fq_rest = find(zout_rest(isuj, ivox,:,:) == max(zout_rest(isuj, ivox,:,:)));

        TF(isuj, ivox) = fth(fq);
        TF_rest(isuj, ivox) = fth(fq_rest);
    end
end



dim = 5; % degree of fit
fmax = 45;

% load data and rest and compute R 
rsq = [];
r2z   = @(r)0.5*log((1+r)./(1-r));

for z0 = 1:2

    if z0 == 1 % task
        X = load('./Data/Ps_task/timefreq_fft_211015_1605.mat'); % source_Ps_az Frequency gradient complete data (histogram) in bst
    elseif z0 == 2 % resting
        X = load('./Data/Ps_rest/timefreq_fft_211015_1606.mat'); % source_Ps_rest_az Frequency gradient complete data (histogram)  in bst
    end
    
    for isuj = 1:size(out,1)

        if z0 == 1 % task
            D = TF(isuj,:)'; % TF or TF_rest (loaded above in section Create Gradient per subject)
        elseif z0 == 2 % resting
            D = TF_rest(isuj,:)'; % TF or TF_rest (loaded above in section Create Gradient per subject)
        end

        X.GridLoc(D>fmax,:) == 0;
        D(D>fmax) = 0;

        cubicFit = []; cubicCoef = [];
        for i0 = 1:3
            IN = X.GridLoc(:,i0);
            [cubicCoef,stats,ctr] = polyfit(D, IN, dim);
            cubicFit(i0, :) = polyval(cubicCoef, D, [], ctr);

            n = length(X.GridLoc); % 1673
            p = polyfit(D, IN, dim);
            f = polyval(p, D);
            df = length(p);
            SSE1 = sum((IN-f).^2); 
            SSE2 = (n-1) * var(IN);  
            rsq(z0,isuj,i0)  = r2z(1 - SSE1/SSE2 * (n-1) /(n-df)); % R-square adjusted to z
            clear n p f df SSE1 SSE2 IN
        end
    end
end

% Stats
fprintf('\n1.Parametric repeated-measure Anova on performance corr - incorr \n')
tmp = permute(rsq,[2 3 1]);
[efs,F,cdfs,p,eps,dfs,b,y2,sig] = repanova(reshape(tmp, size(tmp,1), []), [2, 3], {'task/rest', 'axes'});

% Figure
figure; set(gcf,'color','w'); colormap([.5 .5 .5; 1 1 1])
subplot(2,4, 1); hold on

x = squeeze(rsq(1,:,:)); % music
    errorbar(1:3, mean(x), std(x)/sqrt(size(x,1)), ...
    'b.', 'color', [0 0 0], 'MarkerSize',15, 'Linewidth', 1.5,'Capsize', 1)
    
x = squeeze(rsq(2,:,:)); % rest
    errorbar(1:3, mean(x), std(x)/sqrt(size(x,1)), ...
    'r.', 'color', [.5 .5 .5], 'MarkerSize',15, 'Linewidth', 1.5,'Capsize', 1)

% esthetics      
hold off
set(gca, 'FontSize', 11, 'FontName', 'Calibri') % 10 or 11
xlabel( [] )
ylabel( 'r²' )
set(gca,'Layer','top','Box','off','TickLength',[.01 .01])
ylim( [0, .6] ) %!
yticks( 0:.2:1 )  
    xticks( 1:3 )
    xticklabels( {'Y', 'X', 'Z'} )
    xlim( [.5, 3.5] )
    
% legend({'Music', 'Rest'}, 'box', 'off')
    
helper_ILL(gcf, 0.8, 12, false);
print('./Figures/Fig3d1','-painters', '-dpdf');

%% Figure 3d2

clearvars
close all
clc

cd('C:\Users\arnau\Desktop\Groove_final')

% load PS for gradient 
load('.\Data\Ps_task\Ps_task_gradient_subj_3cond.mat', 'out', 'fth');
G = load('./Data/Ps_task/timefreq_fft_211015_1605.mat'); % source_Ps_az Frequency gradient complete data (histogram) in bst

% load('.\Data\Ps_rest\Ps_rest_gradient_subj_3cond.mat', 'out_rest');

nvox = 1673;
nsuj = size(out,1);
MEG = nan(size(out,1), 3 , nvox);

for isuj = 1:size(out,1)

    TF_suj = squeeze(out(isuj,:,:,:));
    TF_suj = permute(TF_suj, [1 3 2]);
    
    for icond = 1:3
        tmp = squeeze(zscore_az(TF_suj(:,:,icond), [], 1));
        for ivox = 1:nvox
            fq = [];
            fq = find(tmp(ivox,:) == max(tmp(ivox,:)));
            MEG(isuj,icond,ivox) = fth(fq);
        end
    end
end

% parameters
dim = 5; % degree of fit
fmax = 45;

FIT = nan(nsuj, 3, 3);
for isuj = 1:nsuj
    disp(['participant ' num2str(isuj) '/' num2str(nsuj)]);
    for icond = 1:3
        D = squeeze(MEG(isuj,icond,:)); % MEG (loaded above)

        G.GridLoc(D>fmax,:) == 0;
        D(D>fmax) = 0;

        cubicFit = []; cubicCoef = [];
        for i0 = 1:3
            IN = G.GridLoc(:,i0);
            [cubicCoef,stats,ctr] = polyfit(D, IN, dim);
            cubicFit(i0, :) = polyval(cubicCoef, D, [], ctr);

            n = length(G.GridLoc); % 1673
            p = polyfit(D, IN, dim);
            f = polyval(p, D);
            df = length(p);
            SSE1 = sum((IN-f).^2); 
            SSE2 = (n-1) * var(IN);  
            FIT(isuj,icond,i0)  = 1 - SSE1/SSE2 * (n-1) /(n-df); % R-square adjusted
            clear n p f df SSE1 SSE2 IN
        end
        clear cubicCoef stats ctr cubicFit D
    end
end
clear G dim fmax

% stat
[~,F,~,~] = repanova(reshape(FIT, nsuj, []), [3, 3], {'Spatial dim', 'Condition'}); % suj,cond,YXZ

% Figure
figure; set(gcf,'color','w'); colormap([.5 .5 .5; 1 1 1])
subplot(2,4, 1); 
hold on
 
x = squeeze(FIT(:,1,:)); % cond1
    errorbar( -.2+(1:3), mean(x), std(x)/sqrt(size(x,1)), ...
    'b.', 'color', [0 0 0], 'MarkerSize',15, 'Linewidth', 1.5,'Capsize', 1)
    
x = squeeze(FIT(:,3,:)); % cond2
    errorbar( 0+(1:3), mean(x), std(x)/sqrt(size(x,1)), ...
    'b.', 'color', [.4 .4 .4], 'MarkerSize',15, 'Linewidth', 1.5,'Capsize', 1)
    
x = squeeze(FIT(:,2,:)); % cond3
    errorbar( +.2+(1:3), mean(x), std(x)/sqrt(size(x,1)), ...
    'b.', 'color', [.8 .8 .8], 'MarkerSize',15, 'Linewidth', 1.5,'Capsize', 1)
    
% esthetics      
hold off
set(gca, 'FontSize', 11, 'FontName', 'Calibri') % 10 or 11
xlabel( [] )
ylabel( 'r²' )
set(gca,'Layer','top','Box','off','TickLength',[.01 .01])
ylim( [0, .6] ) %!
yticks( 0:.2:1 )  
xticks( 1:3 )
xticklabels( {'Y', 'X', 'Z'} )
xlim( [.5, 3.5] )
    

helper_ILL(gcf, 0.8, 12, false);
print('./Figures/Fig3d2','-painters', '-dpdf');

%% Figure S2a loading data

% initialize
clearvars
close all
clc
SUJETS  = {'002','003','004','005','006','007','008','009','010','012'...
           '014','015','016','017','018','019','020','021','022','023','024','025' ...
           '026','027','028','029','030','031','032'}; % (001&013out) without 011
       
Cond    = {'Low', 'Medium', 'High'};

[~, hname] = system('hostname');
if strcmp(hname(1:2),'bn')  == 1
    Dir = '/Users/Bn/';
elseif strcmp(hname(1),'R') == 1
    Dir = 'C:\Users\arnau\';
elseif strcmp(hname(1),'C') == 1
    Dir = 'C:\Users\zalta\';
end
addpath([ Dir 'Dropbox\tools\'])

% parameters
nvox = 1673; ntrial = 144; nsuj = numel(SUJETS); nmel = 12;

% load MEG data
out = NaN(numel(SUJETS), nvox, 100, 3);

z = '\source_Ps_az.mat';

for isuj = 1:nsuj
    disp(['participant ' num2str(isuj) '/' num2str(nsuj)]);
    
    % 1. Load TF data (Ps_az)
    load([ DriveName('Maxtor') ':\Groove_script\Results\groove_' SUJETS{isuj} z]); %!
    
    % reshape - avg 3 components
    y_cur = permute(TF, [2 3 1]);
    y_cur = reshape(y_cur, 3, nvox, length(fth), 3, ntrial/3);
    y_cur = squeeze(mean(y_cur,1)); %  voxel|freq|cond|trials
    
    % avg melodies (n=48)
    out(isuj,:,:,:) = squeeze(mean(y_cur,4));
    
    clear TF y_cur
end

% extact gradient (group avg)
MEG = nan(nvox, 3);

%  zscore - avg subjects - extract freq
for icond = 1:3
    tmp = out(:,:,:,icond);
    tmp = squeeze(zscore_az(tmp, [], 2));
    tmp = squeeze(mean(tmp,1));
    for ivox = 1:nvox
        fq = find(tmp(ivox,:) == max(tmp(ivox,:)));
        MEG(ivox,icond) = fth(fq);
    end
    clear tmp fq
end

dim = 5; % degree of fit
fmax = 45;
cond = ["Low" "Medium" "High"];

% load data from bst : subject groove 011 present
X = load('C:/Users/arnau/Desktop/Groove_final/Data/Ps_task/timefreq_fft_211015_1605.mat');

%% Figure S2a display figure

for icond = 1:length(cond) % select 1 condition
    
    X.TF = MEG(:,icond);

    fig = figure; set(gcf,'color','w');
    fig.Position(3) =  250; fig.Position(4) =  250;

    rsq = [];
    X.GridLoc(X.TF>fmax,:) == 0;
    X.TF(X.TF>fmax) = 0;

    cubicFit = []; cubicCoef = [];
    for i0 = 1:3
        IN = X.GridLoc(:,i0);
        [cubicCoef,stats,ctr] = polyfit(X.TF, IN, dim);
        cubicFit(i0, :) = polyval(cubicCoef, X.TF, [],ctr);

        n = length(X.GridLoc); % 1673
        p = polyfit(X.TF, IN, dim);
        f = polyval(p, X.TF);
        df = length(p);

        SSE1 = sum((IN-f).^2); 
        SSE2 = (n-1) * var(IN);  
        rsq(i0)  = 1 - SSE1/SSE2 * (n-1) /(n-df); % R-square adjusted
        clear n p f df SSE1 SSE2 IN
    end

    ax1 = axes;
    scatter3(ax1, cubicFit(1, :), cubicFit(2, :), cubicFit(3, :), 35, X.TF, 'filled','MarkerEdgeColor', 'k');
    axis equal

    ax1.YLim = [-0.03 0.03];
    ax1.XLim = [-0.07 0.07];
    ax1.ZLim = [0.02 0.1];

    ax1.XLabel.String = 'Y (mm)';
    ax1.YLabel.String = 'X (mm)'; ax1.YLabel.Rotation = -35; ax1.YLabel.Position = [-0.1 -0.04 0.01];
    ax1.ZLabel.String = 'Z (mm)';

    view(-27,20)
    
    title(cond(icond))

    % for poster
    ax1.GridLineStyle = ':';
    ax1.YTick = [-0.03 0 0.03];
    ax1.XTick = [-0.06 -.03 0 0.03 0.06];
    ax1.ZTick = [0.03 0.06 0.09];
    ax1.YTickLabel = {'-30' '0' '30'};
    ax1.XTickLabel = {'-60' '-30' '0' '30' '60'};
    ax1.ZTickLabel = {'30' '60' '90'};
    ax1.TickLength = [0.002 0.002];
    colormap(ax1,'parula');

    set(gca, 'FontSize', 11, 'FontName', 'Calibri') % 10 11 or 12 ?

helper_ILL(gcf, 1.2, 4, false);
print(['C:/Users/arnau/Desktop/Groove_final/Figures/FigS2_' Cond{icond}] ,'-painters', '-dpdf');
end

%% Figure S2b loading data

% initialize
clearvars
close all
clc
SUJETS  = {'002','003','004','005','006','007','008','009','010','012'...
           '014','015','016','017','018','019','020','021','022','023','024','025' ...
           '026','027','028','029','030','031','032'}; % (001&013out) without 011
       
Cond    = {'Low', 'Medium', 'High'};

[~, hname] = system('hostname');
if strcmp(hname(1:2),'bn')  == 1
    Dir = '/Users/Bn/';
elseif strcmp(hname(1),'R') == 1
    Dir = 'C:\Users\arnau\';
elseif strcmp(hname(1),'C') == 1
    Dir = 'C:\Users\zalta\';
end
addpath([ Dir 'Dropbox\tools\'])

% parameters
nvox = 1673; ntrial = 144; nsuj = numel(SUJETS); nmel = 12;


% Load regressors (outvar = IN)
% IN1 = groove ; IN2 = sync ; IN3 = indexes

% 1. load Bhv: stimulus index, Grooviness, Syncopation (polyphonic)
Rpf = nan(nsuj, 3, 4, 12); Rid = Rpf;
for isuj = 1:nsuj
    subject = ['groove_' SUJETS{isuj}];
    cd([Dir 'Dropbox/Groove/3b logfiles/' subject])
    
    % load data
    for irun = 1:4
        raw = importdata([SUJETS{isuj} '-groove' num2str(irun) '.log']);
        stim = raw.textdata(:,2); %!
        stim = str2double((stim(3:end)));
        stim = stim(~isnan(stim));
        stim(stim == 100) = [];
        resp = raw.data(:,1); %!
        resp(isnan(resp)) = [];
        
        % organise fn condition
        for idata = 1:length(Cond)
            x = resp(round(stim/1e2)==idata);
            Rpf(isuj,idata,irun,1:length(x)) = x;
            x = stim(round(stim/1e2)==idata);
            Rid(isuj,idata,irun,1:length(x)) = x;
        end
        clear stim resp idata x raw
    end
end

% reshape to fit with Bs3 orga: bloc then trials
Rpf = reshape(Rpf, nsuj, length(Cond), []);
Rid = reshape(Rid, nsuj, length(Cond), []);

% transform id to get trial index (remove condition index = hundreds)
Rid = mod(Rid,100);

% fill NaNs of participant 029
Rpf(isnan(Rpf)) = 1; Rid(isnan(Rid)) = 3; %!

% load Syncopation index data
cd([Dir 'Dropbox/Groove/1 stimuli/'])
X = xlsread('Index de syncopation.xlsx');
X = X(:,2); % use poly-phonic index (2)
X = reshape(X, 3, length(X)/3);
sync = [X(3,:); X(1:2,:)]; % R-G-S
% sync(2,id2(1,2,1))

% organize data fn suject-specific presentation
IN{1} = Rpf;
for isuj = 1:nsuj
    for idata = 1:length(Cond)
        for itrial = 1:size(Rpf,3)
            x = Rid(isuj,idata,itrial);
            IN{2}(isuj,idata,itrial) = sync(idata,x);
        end
    end
end
IN{3} = Rid;
clear Rpf Rid Rid2 sync Ed  X temp fq out x  irun idata isuj itrial subject ac


% load MEG data & define gradient (outvar = MEG)

% loop across subjects
z = '\source_Ps_az.mat';
MEG = nan(nsuj, 3, nmel , nvox);
for isuj = 1:nsuj
    disp(['participant ' num2str(isuj) '/' num2str(nsuj)]);
    
    % 1. Load TF data (Ps_az)
    load([ DriveName('Maxtor') ':\Groove_script\Results\groove_' SUJETS{isuj} z]); %!
    
    % reshape - avg 3 components
    y_cur = permute(TF, [2 3 1]);
    y_cur = reshape(y_cur, 3, nvox, length(fth), 3, ntrial/3);
    y_cur = squeeze(mean(y_cur,1)); %  voxel|freq|cond|trials
    
    % avg 4 repetitions - zscore - extract freq
    for icond = 1:3
        for imel = 1:nmel
            out = mean(y_cur (:, :, icond, squeeze(IN{3}(isuj,icond,:))==imel ), 4);
            out = squeeze(zscore_az(out, [], 1));
            for ivox = 1:nvox
                fq = find(out(ivox,:) == max(out(ivox,:)));
                MEG(isuj,icond,imel,ivox) = fth(fq);
            end
            clear out fq
        end
    end
    clear TF y_cur
end

GRO = nan(nsuj, 3, nmel);
SYN = nan(nsuj, 3, nmel);
for isuj = 1:nsuj
for icond = 1:3
for imel = 1:nmel        
    GRO(isuj,icond,imel) = mean(IN{1} (isuj, icond, squeeze(IN{3}(isuj,icond,:))==imel ), 3);
    SYN(isuj,icond,imel) = mean(IN{2} (isuj, icond, squeeze(IN{3}(isuj,icond,:))==imel ), 3);
end
end
end
clear IN

% fit MEG data (outvar = rsq)

% parameters
dim = 5; % degree of fit
fmax = 45;

% load data from bst (GridLoc)
G = load('C:/Users/arnau/Desktop/Groove_final/Data/Ps_task/timefreq_fft_211015_1605.mat'); % source_Ps_az Frequency gradient complete data (histogram) in bst

FIT = nan(nsuj, 3, nmel, 3);
for isuj = 1:nsuj
    disp(['participant ' num2str(isuj) '/' num2str(nsuj)]);
    for icond = 1:3
    for imel  = 1:nmel
        D = squeeze(MEG(isuj,icond,imel,:)); % MEG|TF (loaded above)

        G.GridLoc(D>fmax,:) == 0;
        D(D>fmax) = 0;

        cubicFit = []; cubicCoef = [];
        for i0 = 1:3
            IN = G.GridLoc(:,i0);
            [cubicCoef,stats,ctr] = polyfit(D, IN, dim);
            cubicFit(i0, :) = polyval(cubicCoef, D, [], ctr);

            n = length(G.GridLoc); % 1673
            p = polyfit(D, IN, dim);
            f = polyval(p, D);
            df = length(p);
            SSE1 = sum((IN-f).^2); 
            SSE2 = (n-1) * var(IN);  
            FIT(isuj,icond,imel,i0)  = 1 - SSE1/SSE2 * (n-1) /(n-df); % R-square adjusted
            clear n p f df SSE1 SSE2 IN
        end
        clear cubicCoef stats ctr cubicFit D
    end
    end
end
clear G dim fmax

[~,F,~,~] = repanova(reshape(FIT, nsuj, []), [3, 36], {'Spatial dim', '36 melodies'}); % suj,cond,YXZ

% repanova( squeeze(mean(FIT,2)), [3], {'X'}); % suj,cond,melo,YXZ

%% Figure S2b display figure

figure; set(gcf,'color','w');
cmap = [0 0 0 ; .4 .4 .4; .8 .8 .8];

cond = ["Y" "X" "Z"];

for i0 = 1:length(cond)

    subplot(2,3,i0)
    hold on
   
x = permute(mean(FIT(:,:,:,i0) ,1), [3 2 1]); % suj avg   
        for z0 = 1:3
            % boxplot(x,'boxstyle','outline', 'colors','k','OutlierSize',4,'Symbol','','Widths',0.2)
            scatter(ones(nmel,1).*(z0+(rand(nmel,1)-0.5)/10), x(:,z0), ...
                'MarkerFaceColor',cmap(z0,:),'MarkerEdgeColor','k',...
                'MarkerFaceAlpha',.9,'MarkerEdgeAlpha',.8)
        end

% % per subject        
%     for z0 = 1:3
%         x = squeeze(FIT(:,z0,:,i0)); % per suj
    
%         h = errorbar(ones(nmel,1).*(z0+(rand(nmel,1)-0.5)/2), ...
%                      squeeze(mean(x)), squeeze(std(x))/sqrt(nsuj), ...
%         '.', 'color', cmap(z0,:), 'MarkerSize',15, 'Linewidth', 1.5,'Capsize', 1);
%         alpha = 0.65;   
%         set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
%     end

title(cond(i0))

ylim( [0, .5])
yticks(.1:.2:.5)
ylabel( 'r²' )

set(gca, 'FontSize', 11, 'FontName', 'Calibri')
xlim( [.5, 3.5] )
    xticks( 1:3 )
    xticklabels( {'Low', 'Medium', 'High'} )
    xtickangle(45)
    set(gca,'Layer','top','Box','off','TickLength',[.01 .01])

end

helper_ILL(gcf, 0.8, 12, false);
print('C:/Users/arnau/Desktop/Groove_final/Figures/FigS2b','-painters', '-dpdf');
% saveas(gcf,'C:\Users\zalta\Dropbox\Groove\6_Arnaud\Groove_au_propre\Figures\Fig_supp_gradient_statistics_per_cond_per_melodies.png')


