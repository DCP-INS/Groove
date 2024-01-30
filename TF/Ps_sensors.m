%% Compute Spectral decompostion at the channel level for task

clearvars
close all
clc

origDir = 'C:\Users\zalta\Documents\GrooveS_downsample\';

SUJETS  = {'002','003','004','005','006','007','008','009','010','011','012'...
        '014','015','016','017','018','019','020','021','022','023','024','025' ...
        '026','027','028','029','030','031','032'}; % (001&013out)
    
Rname   = {'Low', 'Medium', 'High'};

% Select task (istask = true) or resting data (istask = false)
istask = true; 

if logical(istask)
    savename = 'Ps_az_down';
else
%     savename = 'Ps_rest_az_down';
end
    
% Frequencies to compute
fth = logspace(0,2,100); % 100 values, 1-100 Hz log scale

brainstorm nogui; warning off
% ------------------------------------------------------------------------

% computation parameters
nchan  = 299;
ntrial = 48;
O  = []; O.Method = 'morlet'; O.Output = 'all'; O.Comment = 'bidule';
    O.ListFiles = []; O.iTargetStudy = []; O.SensorTypes = [];
    O.TimeVector = []; O.RowNames = []; O.Freqs = []; O.Measure = [];
    O.TimeBands = []; O.MorletFc = 1; O.MorletFwhmTc = 3; O.ClusterFuncTime = 'none';

% Compute Ps for each subject
for isuj = 1:length(SUJETS)
    
    if logical(istask)
        subject = sprintf('groove_%s',SUJETS{isuj});
    else
        subject = sprintf('rest_groove_%s',SUJETS{isuj});
    end
    
    s_nam = ['signal_brut_downsample_' subject];
    disp(s_nam);  
    
    load([ origDir s_nam '.mat']);
    tic
    for idata = 1:length(Rname)         
        for i0 = 1:ntrial % 20sec|trial

            % A: time-frequency - cut - avg otime
            % initialise outputs
            O.TimeVector = tim;
            O.Freqs = fth; O.RowNames = {1:length(fth)}; O.Measure = 'power';   
            
            % TF power
            tf = bst_timefreq(squeeze(s(idata,i0,:,:)) ,O); bst_progress('stop');
            tf = squeeze(tf{1}.TF);            
            
            % cut [.5 16] sec - avg time
            [~,deb] = find(round(tim*5e2)== .5 *5e2);
            [~,fin] = find(round(tim*5e2)== 16 *5e2);
            TF(idata,i0,:,:) = squeeze(mean( tf(:,deb:fin,:), 2));  
            clear tf deb fin
        end

        save(['G:\GrooveS\results\groove_' SUJETS{isuj} '\' savename '.mat'], 'TF')
    end
    clear TF s tim idata i0
    toc
end

clear z0

close brainstorm