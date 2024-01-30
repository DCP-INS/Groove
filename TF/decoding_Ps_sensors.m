%% Compute decoding on Ps sensor

clearvars
close all
clc

% Path
Dir = 'C:\Users\arnau\';
addpath([Dir 'Dropbox\tools\'])

% load MEG data and regressors
load('C:\Users\arnau\Desktop\Groove_final\Data\Chan_data_7regr.mat', 'X'); % X = behav; Y = Meg; fth = freqs

% load bad channels & bad trials
load('C:\Users\arnau\Desktop\Groove_final\Data\bad.mat') % chan
load('C:\Users\arnau\Desktop\Groove_final\Data\bad_bn_process.mat'); % trials

SUJETS  = {'002','003','004','005','006','007','008','009','010','011','012'...
           '014','015','016','017','018','019','020','021','022','023','024','025' ...
           '026','027','028','029','030','031','032'}; % (001&013out)
 
% parameters
ntrial = 144;
nvox = 1673; 
s = [];

for isuj = 1:size(X,1)
    disp(['participant ' num2str(isuj) '/' num2str(size(X,1))]);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x_cur0 -> low(1) vs. med|high (0) -> Periodic
%x_cur1 -> sync
%x_cur2 -> bhv 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x_cur = squeeze(X(isuj,:,:));
    y_cur = squeeze(Y(isuj,:,:,:));
    
    % load source inversion matrix
    cd([DriveName('Groove') ':\GrooveS\data\groove_' SUJETS{isuj} '\TRIGGER_16']) 
    files = dir('results_dSPM*.mat');
    
    for i0 = 1:length(files)
        source = load(files(i0).name);
        if strcmp(source.Comment, 'dSPM: MEG(Unconstr) 2016')
            break
        end
    end
    
    % remove bad trials
    goodtrial = badtrial_bnprocess(isuj,:)';
    goodtrial = find(goodtrial == 1);
    y_cur = y_cur(goodtrial,:,:); 
    x_cur = x_cur(goodtrial,:); % include/exclude bads
 
    % multivariate analysis subject by subject
    nfold = 10 ; %size(x_cur,1); % 10
    alpha = 2; % 1 ridge parameter (to adjust?)
    r2z   = @(r)0.5*log((1+r)./(1-r));
    
    for ifq = 1:length(fth)
        disp(['freq ' num2str(ifq) '/' num2str(length(fth))]);

        for i0 = 1:size(x_cur,2)
            out_dec = decode_ridg(x_cur(:,i0), y_cur(:,:,ifq), nfold, alpha);
            s(isuj,i0,ifq) = r2z(out_dec.rdec); % regressors|fq
            clear out_dec
        end
    end
end
clear iStudy R y_cur x_cur r2z

save(['C:\Users\arnau\Desktop\Groove_final\Data\MEG\Ps\Ps_az_decoding_goodchan_goodtrial_GrooveS_nfold' num2str(nfold) '_alpha' num2str(alpha) '.mat'], 's', 'fth');

clear nfold alpha
disp('TERMINARES');