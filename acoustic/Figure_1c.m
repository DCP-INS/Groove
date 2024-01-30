clearvars
close all
clc

cd('C:/Users/arnau/Desktop/Groove_final')
addpath(genpath('./tools'))

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

% % load Syncopation index data (Vuust v. is not working)
X = xlsread('./Data/stimuli/Index de syncopation.xlsx');

X = X(:,2); % use poly-phonic index (2)
X = reshape(X, 3, length(X)/3);
sync = X([3 1 2],:); % R-G-S (Low|Med|High)
IN{2} = sync;
clear sync X

% load acoustic (Ed): 2 hz - *100
fq_acc = 2;
t = load('./Data/stimuli/accenv_Ed.mat'); % accenv_norm1.mat
[~,b] = min(abs(t.out.fq-fq_acc));
IN{3} = 1e2*t.out.D(:,:,b); 
fq_acc = t.out.fq;
clear t

% show each melody (ie avg subjects)
figure; set(gcf,'color','w');
cmap = [0 0 0 ; .4 .4 .4; .8 .8 .8];

for i0 = 2
    if i0 == 1; x = permute(IN{3}, [2 1]); % 2hz acoustic        
    elseif i0 == 2; x = permute(IN{2}, [2 1]); % syncope
    elseif i0 == 3; x = permute(nanmean(IN{1},1), [3 2 1]); % groove: avg subjects!
    end
    n = length(x);    
    
    hold on
        for z0 = 1:3
            boxplot(x,'boxstyle','outline', 'colors','k','OutlierSize',4,'Symbol','','Widths',0.2)
            scatter(ones(n,1).*(z0+(rand(n,1)-0.5)/10), x(:,z0), 20, ...
                'MarkerFaceColor',cmap(z0,:),'MarkerEdgeColor','k',...
                'MarkerFaceAlpha',.9,'MarkerEdgeAlpha',0)
        end

        if i0 ==1;      ylabel( 'Amplitude (2 Hz)' ); ylim( [0, 12] ); yticks([0 4 8 12])
        elseif i0 == 2; ylabel( 'Degree of syncopation' )
        elseif i0 == 3; ylabel( 'Groove' ); ylim( [1, 5] ); yticks(1:2:5) %!
        end
        
        set(gca, 'FontSize', 11, 'FontName', 'Calibri')
        xlim( [.5, 3.5] )
        xlabel( 'Condition' )
        set(gca,'Layer','top','Box','off','TickLength',[.001 .001])
end

helper_ILL(gcf, 1.1, 3, false);
print('./Figures/Fig1c','-painters', '-dpdf');
