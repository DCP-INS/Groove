%% load PS and PS rest from bst (except suj 23 which has no PS rest)
clearvars
close all
clc

% Choose your target
analysis = 'Ps';

switch analysis
    case 'Ps'
        name = 'source_Ps_az';
        tag  = '/source_Ps_az';
        tag_rest = '/source_Ps_rest_az';
    case 'Ps_zscored'
        name = 'source_Ps_az';
        tag  = '/source_Ps_az';
end        

SUJETS  = {'002','003','004','005','006','007','008','009','010','012'...
           '014','015','016','017','018','019','020','021','022','023','024','025' ...
          '026','027','028','029','030','031','032'}; % (001&013out) without 011
warning off

% les 1001 paths 
[~, hname] = system('hostname');
if strcmp(hname(1),'R') == 1
    Dir = 'C:\Users\arnau\';
    addpath(genpath([Dir 'Dropbox\tools\']))
    origDir = [DriveName('Maxtor') ':\GrooveS_downsample\'];
    load([DriveName('Maxtor') ':/GrooveS_downsample/bad_bn_process.mat']); % trials

elseif strcmp(hname(1:4),'COG8') == 1
    Dir = 'C:\Users\zalta\';
    addpath(genpath([Dir 'Dropbox\tools\']))
    origDir = [DriveName('Groove') ':\\GrooveS\results\'];
    load([Dir 'Desktop\Groove\bad_bn_process.mat']); % trials
    load([Dir 'Desktop\Groove\GoodChannel_ImagingKernel.mat'])
end 

% % Bst 
cd ([Dir '\Downloads\brainstorm3'])
brainstorm nogui

% Export files to matlab
sFiles = bst_process('CallProcess', 'process_select_files_timefreq', [], [], ... % results/timefreq
    'subjectname',   'Group analysis', ...
    'condition',     '', ...
    'tag',           tag, ...
    'includebad',    0, ...
    'includeintra',  1, ...
    'includecommon', 0);
sFiles_rest = bst_process('CallProcess', 'process_select_files_timefreq', [], [], ... % results/timefreq
    'subjectname',   'Group analysis', ...
    'condition',     '', ...
    'tag',           tag_rest, ...
    'includebad',    0, ...
    'includeintra',  1, ...
    'includecommon', 0);  
export_matlab({sFiles.FileName}','Y')
export_matlab({sFiles_rest.FileName}','X')

% parameters
fth = X01.Freqs;
nvox = size(X01.TF, 1);

% Load in Matlab
out = NaN(numel(SUJETS), nvox, 3, numel(fth));
out_rest = out;
for isuj = 1:length(SUJETS)
    if isuj < 10
        y = eval(['Y0' num2str(isuj)]);
        x = eval(['X0' num2str(isuj)]);
    elseif isuj >9 && isuj < 21
        y = eval(['Y' num2str(isuj)]);
        x = eval(['X' num2str(isuj)]);
    elseif isuj > 20
        y = eval(['Y' num2str(isuj)]);
        x = eval(['X' num2str(isuj-1)]); % hence suj20_rest is duplicated
    end
    out(isuj,:,:,:)      = y.TF;
    out_rest(isuj,:,:,:) = x.TF;
    % y.Comment
    clear y x
end
clear Y* X*

save('C:\Users\arnau\Desktop\Groove_final\Data\Ps_task\Ps_task_gradient_subj_3cond.mat', 'out', 'fth');
save('C:\Users\arnau\Desktop\Groove_final\Data\Ps_rest\Ps_rest_gradient_subj_3cond.mat', 'out_rest', 'fth');
