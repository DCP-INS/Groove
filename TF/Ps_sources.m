clearvars
close all
clc


origDir = 'F:\GrooveS_downsample\';
addpath('C:\Users\zalta\Dropbox\tools')
SUJETS  = {'002','003','004','005','006','007','008','009','010','011','012','014', '015',...
        '016','017''018','019','020','021','022','023','024','025' ...
        '026','027','028','029','030','031','032'}; % (001&013out)

istask = true;

if logical(istask) 
    savename= 'source_Ps_az'; % PsdB|HGdB
else
    savename= 'source_Ps_rest_az'; % PsdB|HGdB
end

fth = logspace(0,2,100); % 100 values, 1-100 Hz log scale

cd C:\Users\zalta\Downloads\brainstorm3
brainstorm nogui; warning off
% ------------------------------------------------------------------------
% computation parameters
nchan  = 299; ntrial = 144; nsourc = 5019;
O  = []; O.Method = 'morlet'; O.Output = 'all'; O.Comment = 'bidule';
    O.ListFiles = []; O.iTargetStudy = []; O.SensorTypes = [];
    O.TimeVector = []; O.RowNames = []; O.Freqs = []; O.Measure = []; %!
    O.TimeBands = []; O.MorletFc = 1; O.MorletFwhmTc = 3; O.ClusterFuncTime = 'none';
    
    % initialise TF parameters
    O.Freqs = fth; O.RowNames = {1:length(fth)}; O.Measure = 'power';
    
    if logical(istask) 
        load([origDir 'signal_brut_downsample_groove_002.mat'], 'tim' );
    else
        load([origDir 'signal_brut_downsample_rest_groove_002.mat'], 'tim' );
    end

    O.TimeVector = tim;
    clear tim
    
% loops
parfor isuj = 1:length(SUJETS)
    
tic    
    TF = zeros(ntrial, nsourc, numel(fth));

    if logical(istask) 
        subject = ['signal_brut_downsample_groove_' SUJETS{isuj}];
    else
        subject = ['signal_brut_downsample_rest_groove_' SUJETS{isuj}];
    end
    
    disp(subject);  
    
    % Load data
    X = load([ origDir subject '.mat'], 's');
    s = reshape(X.s, [], size(X.s, 3), size(X.s, 4));
    
    % load source inversion matrix
    source = [];
    cd(['G:/GrooveS/data/groove_' SUJETS{isuj} '\TRIGGER_16']) 
    files = dir('results_dSPM*.mat');
    for i0 = 1:length(files)
        source = load(files(i0).name);
        if strcmp(source.Comment, 'dSPM: MEG(Unconstr) 2016'); break; end
    end
    
    s = s(:,source.GoodChannel',:);
   
    disp(isuj);
    for i0 = 1:ntrial % 20sec|trial
        Data = squeeze(s(i0,:,:));
       
        % loop over vertex
        for v0 = 1:nsourc
            % TF decomposition
            tf = bst_timefreq( source.ImagingKernel(v0,:)*Data, O); bst_progress('stop');
            tf = squeeze(tf{1}.TF);

            % cut [.5 16] sec - avg time
            [~,deb] = min(abs(O.TimeVector - 0.5));
            [~,fin] = min(abs(O.TimeVector - 16));
            TF(i0,v0,:) = mean( tf(deb:fin,:), 1);
        end
    end
    
% save data
parsave(['G:\GrooveS\results\groove_' SUJETS{isuj} '\' savename '.mat'], TF, fth);
toc
end