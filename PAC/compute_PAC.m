%% Compute the 1.4 or 2 Hz Phase / whole spectrum amplitude coupling
% clear workspace
clearvars
close all
clc

% cd data path
cd G:\GrooveS_downsample\

SUJETS  = {'002','003','004','005','006','007','008','009','010','011','012', ...
           '014','015','016','017','018','019','020','021','022','023','024','025', ...
           '026','027','028','029','030','031','032'};
       
Rname   = {'Low', 'Medium', 'High'};

% for control with PAC 2 Hz
savename = 'source_Pac_1.4Hz';

fth = [1.4, 2, 2.9, 3.9, 4.9, 6.1, 7.3, 8.6, 10.1, 11.7, 13.4, 15.2, 17.2, 19.3, 21.7, 24.2, 26.9, 29.9,... 
33.1, 36.6, 40.3, 44.4, 45]; % Pac, 30fq

% brainstorm
addpath('C:\Users\zalta\Downloads\brainstorm3')
brainstorm nogui; warning off

% computation parameters
ntrial = 48; nsourc = 5019; % ? 1673
O  = []; O.Method = 'morlet'; O.Output = 'all'; O.Comment = 'bidule';
    O.ListFiles = []; O.iTargetStudy = []; O.SensorTypes = [];
    O.TimeVector = []; O.RowNames = []; O.Freqs = []; O.Measure = [];
    O.TimeBands = []; O.MorletFc = 1; O.MorletFwhmTc = 3; O.ClusterFuncTime = 'none';

for isuj = 1:length(SUJETS)
    
    % load data
    subject = ['signal_brut_downsample_groove_' SUJETS{isuj}];
    disp(subject);  

    load([subject '.mat']); % s, tim
    
    % load source inversion matrix
    cd(['G:/GrooveS/data/groove_' SUJETS{isuj} '\TRIGGER_16']) 
    files = dir('results_dSPM*.mat');
    
    for i0 = 1:length(files)
        source = load(files(i0).name);
        if strcmp(source.Comment, 'dSPM: MEG(Unconstr) 2016'); break; end
    end
    s = s(:,:,source.GoodChannel,:);
    
    % initialise TF parameters
    O.Freqs = fth; O.RowNames = {1:length(fth)}; O.Measure = 'none';
    O.TimeVector = tim;

    for idata = 1:length(Rname)
        for i0 = 1:ntrial
            Data = squeeze(s(idata,i0,:,:));

            % initialise outputs
            if (~exist('Pac', 'var'))
                Pac.x = zeros(length(Rname), ntrial, nsourc, length(fth));
                Pac.n = zeros(length(Rname), ntrial, nsourc, length(fth));
                Pac.m = zeros(length(Rname), ntrial, 1);
                Pac.f = fth;
            end

            % loop over vertex
            for v0 = 1:nsourc
                
                % TF decomposition
                tf = bst_timefreq( source.ImagingKernel(v0,:)*Data, O); bst_progress('stop');
                tf = squeeze(tf{1}.TF);

                % cut [.5 16] sec - avg time
                [~,deb] = min(abs(tim - 0.5));
                [~,fin] = min(abs(tim - 16));

                % load complex data matrix
                Phi = angle(tf(deb:fin, 1)); % [-pi, +pi] SELECT 1.4 Hz FREQ

                % kernel - loop on vertices - compute TF power - Pac (sum otime)
                Pac.x(idata,i0,v0,:) = sum( bst_bsxfun(@times, abs(tf(deb:fin, :)), exp(1i*Phi) ), 1); % Phase of 1.4 or 2 Hz
                Pac.n(idata,i0,v0,:) = sum( abs(tf(deb:fin, :).^2), 1);                                % Power of all the frequencies
                Pac.m(idata,i0) = length(tf(deb:fin, 1));
            end
        end
    end
    
    % save data
    save(['G:\GrooveS\results\groove_' SUJETS{isuj} '\' savename '.mat'], 'Pac', 'fth')
    clear Pac  s tim source files  idata i0 v0
end

% load files in bst
% zscore over vertex then ttest vs whole Brain mean PAC by frequency

%% spatial profil of the 1.4 Hz Phase / whole spectrum amplitude coupling 

X = load('G:/GrooveS/data/Group_analysis/Pac_bn_1.4Hz/timefreq_fft_pthresh_211001_1426.mat'); % 'Pac_bn_1.4Hz: zovertex ttest trial by trialavgthencenter | alpha=0.005 (FDR:1,3)' in bst

X.TF = X.TF.* 0;
for i0 = 1:length(scout)
    X.TF(scout{i0},:,i0) = 1;
end

X.TF(:,:,(length(scout)+1):end) = [];
X.Comment = 'PAC 1.4 ROI in the good order';
X.Freqs = 1:4;

% Load in bst then import surfaces then plot the blob then add fig to spectral profil
scout_Name = {'Right Frontal' 'Left Parietal'};

% Procedure to orient figure correctly in bst
% % right click in the font of the figure -> "Figure" -> "Change Backgroung color" -> "RVB"
% % Image Font color -> [0 0 0] = Black
% % 
% % cortex_15002V
% %     Sulci on
% %     color        -> [240 240 240]
% %     Transparency -> 20%
% %     Smooth       -> 10%
% %     
% % The scouts
% %     color        -> [100 100 100]
% %     Transparency -> 20%
% %     Smooth       -> 75%
% %     
% % Orientation for each scout
% % 
% %     'lPAC'      -> view(172,5)
% %     'rPAC'      -> view(8,5)
% %     'SMA'       -> view(270,90)
% %     'rMotor'    -> view(8,5)
% %     'lParietal' -> view(172,5)
% % 
% % saveas(gcf,'C:\Users\zalta\Desktop\Blob_Groove\blob scouts\Contrast_Frites_Atlas\Image scouts\lPAC.png')

%% plot spectral profil of the ROI

clearvars
close all
clc

cd('C:/Users/arnau/Desktop/Groove_final')

% The good voxels (voir plus haut)
load('./Data/scouts_Pac_bn.mat', 'scout', 'scout_Name')

% load Pac bn & Pac bn 1.4 Hz
Pac14 = load('./Data/Pac_bn_1.4Hz_3_computes_brut.mat', 'fth', 'OUT_onebyone');
Pac2  = load('./Data/Pac_bn_3_computes_brut.mat','OUT_onebyone');

fth = Pac14.fth(3:end);

% Choose Pac computing you want to plot
out14 = Pac14.OUT_onebyone(:,:,:,3:end);
out14 = nanmean(out14, 2);
out14 = permute(out14, [1 3 2 4]);

out2 = Pac2.OUT_onebyone(:,:,:,2:end);
out2 = nanmean(out2, 2);
out2 = permute(out2, [1 3 2 4]);

% Smooth the lines
Ip = 100;
xn = linspace(1, fth(end), fth(end) * Ip);

% Blobs without zovertex
for i0 = 2 % 1:length(scout) %[7 8 10 12]

    % Statistics
    pval = 0.005; %!

    for z0 = 1:length(fth)
        [h, p] = ttest(squeeze(mean( out2(:,scout{i0}, :,z0), 2)) ,squeeze(mean( out2(:,: , :,z0), 2)) , pval, 'both');
        p_tot1(z0) = p;

        [h, p] = ttest(squeeze(mean( out14(:,scout{i0}, :,z0), 2)) ,squeeze(mean( out14(:,: , :,z0), 2)) , pval, 'both');
        p_tot2(z0) = p; 
        clear h p
    end

    [h] = mafdr(p_tot1); % , 'BHFDR', 'true'
    z1 = NaN(1, length(fth)); z1(h<pval) = 0; % fth(z1 == 0)

    [h] = mafdr(p_tot2); % , 'BHFDR', 'true'
    z2 = NaN(1, length(fth)); z2(h<pval) = 0;

    % Figure
    fig = figure; set(gcf,'color','w'); colormap([.5 .5 .5; 1 1 1]);
    fig.Position = fig.Position./sqrt(2);

    % image of the scout blob
%     img = imread(['./Data/Blob_Groove/blob PAC/images blob\ROIs/' scout_Name{i0} '.png']);
%     image(img, 'XData', [22 130], 'YData', [.05 .035])

    y = squeeze(mean( out2(:,scout{i0}, :,:), 2));
    y = spline(fth, y, xn);
    boundedline(xn, mean(y,1), nanstd(y,[],1)/sqrt(size(y,1)), ...
            'k.', 'transparency', 'alpha', .25);

    y = [];
    y = squeeze(mean( out2(:,:, :,:), 2));
    y = spline(fth, y, xn);
    boundedline(xn, mean(y,1), nanstd(y,[],1)/sqrt(size(y,1)), ...
            'k:', 'transparency', 'alpha', .25);   


    stat2 = plot(fth, z1+0.011, 'Color', 'k', 'Linewidth', 3, 'DisplayName','Phase 2 Hz');


    x = squeeze(mean( out14(:,scout{i0}, :,:), 2));
    x = spline(fth, x, xn);
    boundedline(xn, mean(x,1), nanstd(x,[],1)/sqrt(size(x,1)),'.', ...
            'cmap', [.94 .5 .04], 'transparency', 'alpha', .25);
            
    x = [];      
    x = squeeze(mean( out14(:,:, :,:), 2));
    x = spline(fth, x, xn);
    [hl,hp] = boundedline(xn, mean(x,1), nanstd(x,[],1)/sqrt(size(x,1)), ':',...
            'cmap', [.94 .5 .04], 'transparency', 'alpha', .25);
            
    stat14 = plot(fth, z2+0.012, 'Color', [.94 .5 .04], 'Linewidth', 3, 'DisplayName','Phase 1.4 Hz');

    set(gca, 'FontSize', 11, 'FontName', 'Calibri')
    xlabel( 'Freq. for amplitude (Hz)' )
    ylabel('PAC')         
          
    yticks([.01 .03 .05])
    xticks([3 5 10 20 45])

    ylim([0.01 .05])
    xlim([fth(1) 45])
    
    set(gca,'Layer','top','Box','off','TickLength',[.002 .002])
    set(gca,'XScale', 'log')
    legend([stat14 stat2], 'EdgeColor', 'w','Location','SouthEast')

    helper_ILL(gcf, 1.29, 6.43);
    print('./Figures/Fig5e','-painters', '-dpdf');

%     saveas(gcf, ['C:\Users\zalta\Dropbox\Groove\6_Arnaud\Groove_au_propre\Figures\Fig5_PAC_' scout_Name{i0} '_under_45.png'])
end