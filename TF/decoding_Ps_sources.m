clearvars
close all
clc

SUJETS  = {'002','003','004','005','006','007','008','009','010','011','012'...
           '014','015','016','017','018','019','020','021','022','023','024','025' ...
           '026','027','028','029','030','031','032'}; % (001&013out)


[~, hname] = system('hostname');
if strcmp(hname(1:2),'bn')  == 1
    Dir = '/Users/Bn/';
elseif strcmp(hname(1),'R') == 1
    Dir = 'C:\Users\arnau\';
elseif strcmp(hname(1),'C') == 1
    Dir = 'C:\Users\zalta\';
end

addpath(genpath('/Applications/Work/Softwares/Toolbox/'))
addpath([ Dir 'Dropbox\Groove\5b scripts MotDatdc copy\decoders']);
addpath([ Dir 'Dropbox\tools\'])
addpath([ Dir 'Downloads\brainstorm3\toolbox\io'])

cd([Dir 'Downloads\brainstorm3']);
brainstorm

% parameters
nreg   = 4;
nvox   = 1673;
ntrial = 144;
kclust = 50;

% savename = 'Fft16_k';
savename = 'Source_PS_az_tot_k';
% savename = 'Source_Pac_Loto_k';
% savename = 'Source_Pac_az_light_Loto_k';


%% Compute decoding all freq in one (figure 4a)

z = '\source_Ps_az.mat';

    % do not zscore MEG data, but zscore regressors (xcur)
    nfold = 10; % 10
    alpha = 2; % 1 ridge parameter (to adjust?)
    r2z   = @(r)0.5*log((1+r)./(1-r));

% load regressors
load([DriveName('Maxtor') ':\GrooveS_downsample\MNE_MVPA\Chan_data_7regr.mat'], 'X'); % only regressors 4 sync; 5 groove

%%%%%%%%%%%%%%%%%%%%%POUR RAPPEL%%%%%%%%%%%%%%%%%%%%%%%
%x_cur0 -> low(1) vs. med|high (0)
%x_cur1 -> sync
%x_cur2 -> bhv 
%x_cur3 -> bhv orthogonalized with sync
%x_cur4 -> sync orthogonalized with periodic (xcur0)
%x_cur5 -> bhv othogonalized with synco orthogonalized with periodic (xcur4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% load bad channels & bad trials
load([DriveName('Maxtor') ':\GrooveS_downsample\bad.mat']) % chan
load([DriveName('Maxtor') ':\GrooveS_downsample\bad_bn_process.mat']); % trials

for isuj = 1:length(SUJETS)
    Q = [];
    R = [];
    source = []; y_cur = []; x_cur = [];
    iStudy = [];

    disp(['participant ' num2str(isuj) '/' num2str(length(SUJETS))]);
    
    % 1. Load TF data
    Q =  load([DriveName('Maxtor') ':\Groove_script\Results\groove_' SUJETS{isuj} z]); %!

    fth = Q.fth;
    y_cur = Q.TF;  
    y_cur = permute(y_cur, [2 3 1]);
    y_cur = reshape(y_cur, 3, nvox, length(fth), []);
    
    % output
    OUT = zeros(nvox, nreg);
    
    % regressors
    x_cur = squeeze(X(isuj,:,:));
        x_cur(:,1) = x_cur(:,1)*-1; %!

    % load source inversion matrix (to estimate neighbour voxels)
    cd([DriveName('Groove') ':\GrooveS\data\groove_' SUJETS{isuj} '\TRIGGER_16']) 
    files = dir('results_dSPM*.mat');
    for i0 = 1:length(files)
        source = load(files(i0).name);
        if strcmp(source.Comment, 'dSPM: MEG(Unconstr) 2016'); break; end
    end
    
    % loop over vertices (2sec/vertex)
    for ivox = 1:nvox
        if mod(ivox,100) == 0
            disp(['voxel ' num2str(ivox) '/' num2str(nvox)]);
        end
        iloc = [];  idist = [];  id_cur = []; 
        iloc = source.GridLoc(ivox,:);
        idist = sqrt( sum( bsxfun(@minus,source.GridLoc, iloc) .^ 2, 2) );
        [~,id_cur] = mink(idist, kclust); % n closest voxels
        
            
            out_dec_periodic = []; out_dec_sync = []; out_dec_groove = []; y = []; 
            
            y = reshape(y_cur(:,id_cur,:,:), [], size(y_cur,4));

            out_dec_periodic = decode_ridg(x_cur(:,1), y', nfold, alpha);
            out_dec_sync     = decode_ridg(x_cur(:,2), y', nfold, alpha);
            out_dec_groove   = decode_ridg(x_cur(:,3), y', nfold, alpha);
           
            OUT(ivox,1) = r2z(out_dec_periodic.rdec);
            OUT(ivox,2) = r2z(out_dec_sync.rdec);
            OUT(ivox,3) = r2z(out_dec_groove.rdec);

    end
    
    % save to Bs3 (voxel-time-frequency structure)
    [~, iStudy] = bst_get('Study',file_fullpath(['groove_' SUJETS{isuj} '/@intra/brainstormstudy.mat']));
    cd([DriveName('Groove') ':\GrooveS\data\groove_' SUJETS{isuj} '/@intra/'])
    R = dir('timefreq*.mat'); R = load(R(1).name);
    R = rmfield(R, {'Atlas', 'GridAtlas', 'History', 'Options', 'Std', ...
        'TimeBands', 'RefRowNames', 'TFmask'});
    R.nAvg = 1; R.DataFile = []; 
    R.Time = 1:size(OUT,2); 
    R.Freqs = Q.fth;
    R.TF = OUT;
    R.Comment = ['Decod_3reg_' savename num2str(kclust)]; %!
    db_add(iStudy, R);  
end
fprintf('done decoding.\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Build figure of decoding for each reg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% *param ttest Decod_3reg_Source_Ps_az_all_freqs_in_one_tot_k6 (30s) | alpha=0.005 (FDR:1,3) in bst
sFiles = {'Group_analysis/@intra/timefreq_fft_pthresh_210726_1420.mat'}; 
RawFiles = {'C:\Users\zalta\Desktop\Blob_Groove\blob reg all in one'};
name_reg = {'periodic' 'syncop' 'groove'};
% Start a new report
bst_report('Start', sFiles);

for i0 = 2:3
    % Process: Export to SPM8/SPM12 (volume)
    sFiles = bst_process('CallProcess', 'process_export_spmvol', sFiles, [], ...
        'outputdir',      {RawFiles{1}, 'Nifti1'}, ...
        'filetag',        name_reg{i0}, ...
        'isconcat',       0, ...
        'timewindow',     [i0, i0], ...
        'timedownsample', 0, ...
        'timemethod',     2, ...  % Keep time dimension (4D volume)
        'freq_export',    1, ...
        'voldownsample',  0, ...
        'iscut',          0);
end

% then import your blob with "import surfaces" (right click) in bst
% then display the blob and "cortex_15002V"
% then select color (in RGB [255 255 255], transparency and smoothing and orientation (using view function) for a nice picture

% % % % right click in the font of the figure -> "Figure" -> "Change Backgroung color" -> "RVB"
% % % % Image Font color -> [0 0 0] = Black
% % % % 
% % % % cortex_15002V
% % % %     Sulci on
% % % %     color        -> [240 240 240]
% % % %     Transparency -> 20%
% % % %     Smooth       -> 10%
% % % %     
% % % % Syncop
% % % %     color        -> [75 75 204]
% % % %     Transparency -> 20%
% % % %     Smooth       -> 40%
% % % %     
% % % % Groove
% % % %     color        -> [204 75 75]
% % % %     Transparency -> 20%
% % % %     Smooth       -> 40%
% % % %     
% % % % Orientation
% % % %     Left         -> view(172,5)
% % % %     Right        -> view(8,5)
% % % %     Add internal views -> remove the controlateral hemisphere and corresponding blobs
% % % %     
% % % % 
% % % % saveas(gcf,'C:\Users\zalta\Desktop\Blob_Groove\blob reg all in one\images blob\Left_groove.png')
% % % % 
% % % % 
% % % % saveas(gcf,'C:\Users\zalta\Desktop\Blob_Groove\blob contrast reg all in one\images blob\Contrast_Left_Extern.png')
% % % % 
% % % % Then remove the black background with Gimp routine



%% compute contrast btw syncope and groove all in one

clear; clc

SUJETS  = {'002','003','004','005','006','007','008','009','010','012'...
           '014','015','016','017','018','019','020','021','022','023','024','025' ...
           '026','027','028','029','030','031','032'}; % (001, 011 & 013 out)


[~, hname] = system('hostname');
if strcmp(hname(1:2),'bn')  == 1
    Dir = '/Users/Bn/';
elseif strcmp(hname(1),'R') == 1
    Dir = 'C:\Users\arnau\';
elseif strcmp(hname(1),'C') == 1
    Dir = 'C:\Users\zalta\';
end

addpath(genpath('/Applications/Work/Softwares/Toolbox/'))
addpath([ Dir 'Dropbox\Groove\5b scripts MotDatdc copy\decoders']);
addpath([ Dir 'Dropbox\tools\'])
addpath([ Dir 'Downloads\brainstorm3\toolbox\io'])
% cd [ Dir 'Downloads\brainstorm3'] ; brainstorm

for isuj = 1:length(SUJETS)

% Export files to matlab
    sFiles = bst_process('CallProcess', 'process_select_files_timefreq', [], [], ... % results/timefreq
        'subjectname',   ['groove_' SUJETS{i0}], ...
        'condition',     '', ...
        'tag',           'Decod_3reg_Source_Ps_az_all_freqs_in_one_tot_k6', ...
        'includebad',    0, ...
        'includeintra',  1, ...
        'includecommon', 0);
    
    export_matlab({sFiles.FileName} ,'X')
    out(isuj,:,:) = X.TF(:,3) - X.TF(:,2);

end

% Statistics

[h,p,ci,stats] = ttest(out);

% Export files to matlab
% sFiles = bst_process('CallProcess', 'process_select_files_timefreq', [], [], ... % results/timefreq
%     'subjectname',   'Group analysis', ...
%     'condition',     '', ...
%     'tag',           'paired ttest Decod_Source_Ps_az_allfreqs_Dif_Sync_Groove', ...
%     'includebad',    0, ...
%     'includeintra',  1, ...
%     'includecommon', 0);

sFiles.FileName = 'Group_analysis/@intra/ptimefreq_fft_no_210301_1306.mat';
export_matlab({sFiles.FileName} ,'Y')

Y.pmap = p';
Y.tmap = stats.tstat';
Y.df = stats.df';
Y.Comment = 'paired ttest Decod_Source_Ps_az_allfreqs_Dif_Sync_Groove_without_groove_011';

% then, Import Y in Bst


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Build figure of decoding for each reg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute 2 absolute files to separate syncop to groove specific blobs
% then plot the blobs ('contrast_groove' & 'contrast_syncope') in bst
% for syncop blobs -> import surfaces 'contrast_syncope_right' &'contrast_syncope_left' then merge then in bst

% *paired ttest Decod_Source_Ps_az_allfreqs_Dif_Sync_Groove | alpha=0.005 in bst
sFiles = {'Group_analysis/@intra/timefreq_fft_pthresh_210803_1543.mat'}; 
RawFiles = {'C:\Users\zalta\Desktop\Blob_Groove\blob contrast reg all in one'};
name_reg = {'contrast'};
% Start a new report
bst_report('Start', sFiles);


    % Process: Export to SPM8/SPM12 (volume)
    sFiles = bst_process('CallProcess', 'process_export_spmvol', sFiles, [], ...
        'outputdir',      {RawFiles{1}, 'Nifti1'}, ...
        'filetag',        name_reg, ...
        'isconcat',       0, ...
        'timewindow',     [3, 3], ...
        'timedownsample', 0, ...
        'timemethod',     2, ...  % Keep time dimension (4D volume)
        'freq_export',    1, ...
        'voldownsample',  0, ...
        'iscut',          0);



% then import your blob with "import surfaces" (right click) in bst
% then display the blob and "cortex_15002V"
% then select color (in RGB [255 255 255], transparency and smoothing and orientation (using view function) for a nice picture

% saveas(gcf,'C:\Users\zalta\Desktop\Blob_Groove\blob contrast reg all in one\images blob\Left_groove.png')


%% Decoding sources by frequency bands (figures 

clear; clc

% start pmode (parallel processing toolbox)
pmode start local 15
isuj = labindex+15

% 1. Searchlight decoding
SUJETS = {'002','003','004','005','006','007','008','009','010','011','012'...
    '014','015','016','017','018','019','020','021','022','023','024','025' ...
    '026','027','028','029','030','031','032'}; % (001&013out)

% parameters
% savename = 'Source_PS_az_6fqbands_3reg_k';
dataname = '\source_Ps_az.mat';

% do not zscore MEG data, but zscore regressors (xcur)
nvox = 1673; ntrial = 144; nreg = 3; fqbands = 6;
kclust = 50; % base = 6
nfold = 10; % 10
alpha = 2; % 1 ridge parameter (to adjust?)
r2z   = @(r)0.5*log((1+r)./(1-r));

% computer
[~, hname] = system('hostname');
if strcmp(hname(1),'R') == 1;     Dir = 'C:\Users\arnau\';
elseif strcmp(hname(1),'C') == 1; Dir = 'C:\Users\M\Desktop\Groove';
    
end    
% DropDir = [Dir 'Dropbox\Groove\'];
%     addpath(genpath('/Applications/Work/Softwares/Toolbox/'))
%     addpath([Dir 'Dropbox\Groove\5b scripts MotDatdc copy\decoders']);
%     addpath([Dir 'Dropbox\tools\'])
%     addpath([Dir 'Dropbox\Groove\6_Arnaud\glmfst'])

addpath(genpath('C:\Users\M\Desktop\Groove\Tools\'))
addpath(genpath('C:\Users\M\Desktop\Groove\Groove_script\tools'))
% load data
if strcmp(hname(1),'R') == 1
    load([DriveName('Maxtor') ':\GrooveS_downsample\MNE_MVPA\Chan_data_7regr.mat'], 'X'); % regressors
    load([DriveName('Maxtor') ':/GrooveS_downsample/bad_bn_process.mat'], 'badtrial_bnprocess'); % trials
    load([DriveName('Maxtor') ':/GrooveS_downsample/Freqs_Ps_perband.mat']); % Freqs -> fqbands = [5 9; 14 18; 30 32; 36 49; 50 63; 66 74]; % [1.2 1.5; 1.8 2.2; 3.8 4.2; 5 10; 10 17; 20 30] in Hz; % indexes for fth

else
    load([Dir '\GrooveS_downsample\MNE_MVPA\Chan_data_7regr.mat'], 'X'); % regressors
    load([Dir '/GrooveS_downsample/bad_bn_process.mat'], 'badtrial_bnprocess'); % trials
    load([Dir '/GrooveS_downsample/Freqs_Ps_perband.mat']); % Freqs -> fqbands = [5 9; 14 18; 30 32; 36 49; 50 63; 66 74]; % [1.2 1.5; 1.8 2.2; 3.8 4.2; 5 10; 10 17; 20 30] in Hz; % indexes for fth
end
% Load GridLoc for searchlight
load([ Dir '\Groove_script\GridLoc.mat']);

% Bst 
cd C:\Users\M\Documents\toolboxes\brainstorm_200203\brainstorm3
brainstorm nogui
%%%%%%%%%%%%%%%%%%%%%POUR RAPPEL%%%%%%%%%%%%%%%%%%%%%%%
%x_cur0 -> low(1) vs. med|high (0)
%x_cur1 -> sync
%x_cur2 -> bhv 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% output

for isuj = 1:length(SUJETS)
    disp(['participant ' num2str(isuj) '/' num2str(length(SUJETS))]);

    % 1. Load TF data (Ps_az)
    Y = load([Dir '\Groove_script\Results\groove_' SUJETS{isuj} dataname]);
    TF = Y.TF;
    fth = Y.fth;
    y_cur = permute(TF, [2 3 1]);
    y_cur = reshape(y_cur, 3, nvox, length(fth), []);
    OUT = zeros(nvox, nreg, length(fqbands));

    % regressors
    x_cur = squeeze(X(isuj,:,:));
    x_cur(:,1) = x_cur(:,1)*-1; %!

    % remove bad trials
    goodtrial = badtrial_bnprocess(isuj,:)';
    goodtrial = find(goodtrial == 1);
    y_cur = y_cur(:,:,:,goodtrial); 
    x_cur = x_cur(goodtrial,:); % include/exclude bads

%     % load source inversion matrix (to estimate neighbour voxels)
%     cd(['G:\GrooveS\data\groove_' SUJETS{isuj} '\TRIGGER_16']) 
%     files = dir('results_dSPM*.mat');
%     for i0 = 1:length(files)
%         source = load(files(i0).name);
%         if strcmp(source.Comment, 'dSPM: MEG(Unconstr) 2016'); break; end
%     end

    % loop over vertices (2sec/vertex)
    for ivox = 1:nvox
        iloc = [];  idist = [];  id_cur = []; 
        iloc = GridLoc(ivox,:,isuj);
        idist = sqrt( sum( bsxfun(@minus,GridLoc(:,:,isuj), iloc) .^ 2, 2) );
        [~,id_cur] = mink(idist, kclust); % n closest voxels
        
        for ifq = 1:length(fqbands)
            y = y_cur(:,id_cur,fqbands(ifq,1):fqbands(ifq,2),:);
            y = reshape(y, [], length(goodtrial));
            out_dec_periodic = decode_ridg(x_cur(:,1), y', nfold, alpha);
            out_dec_sync     = decode_ridg(x_cur(:,2), y', nfold, alpha);
            out_dec_groove   = decode_ridg(x_cur(:,3), y', nfold, alpha);

            OUT(ivox,1,ifq) = r2z(out_dec_periodic.rdec);
            OUT(ivox,2,ifq) = r2z(out_dec_sync.rdec);
            OUT(ivox,3,ifq) = r2z(out_dec_groove.rdec);
        end
        
        if mod(ivox,100) == 0
            disp(['voxel ' num2str(ivox) '/' num2str(nvox)]);
            save(['C:/Users/M/Desktop/Groove/Temporaire_Decoding_Ps_per_band/groove_' SUJETS{isuj} '_perband.mat'], 'OUT' , 'Freqs');
        end         
    end
    save(['E:/GROOVE_A_CONSERVER/Temporaire_Decoding_Ps_per_band/groove_' SUJETS{isuj} '_perband.mat'], 'OUT' , 'Freqs');
    
    
end
fprintf('done decoding.\n');

%% Save in bst on AZ computer

savename = 'Source_PS_az_6fqbands_3reg_k';
kclust = 50; % base = 6
OUT = [];
for isuj = 1:length(SUJETS)
    load(['E:/GROOVE_A_CONSERVER/Temporaire_Decoding_Ps_per_band/groove_' SUJETS{isuj} '_perband.mat']);
    
    % save to Bs3 (voxel-time-frequency structure)
    [~, iStudy] = bst_get('Study',file_fullpath(['groove_' SUJETS{isuj} '/@intra/brainstormstudy.mat']));
    cd(['G:\GrooveS\data\groove_' SUJETS{isuj} '/@intra/'])
    R = dir('timefreq*.mat'); R = load(R(1).name);
    R = rmfield(R, {'Atlas', 'GridAtlas', 'History', 'Options', 'Std', ...
        'TimeBands', 'RefRowNames', 'TFmask'});
    R.nAvg = 1; R.DataFile = []; R.Measure = 'other';
    R.Time = 1:nreg; 
    R.Freqs = Freqs;
    R.TF = OUT;
    R.Comment = ['Decod_' savename num2str(kclust)];
    db_add(iStudy, R);

end

%% Then...: Project - statistics
% savename = 'Source_PS_az_tot_k';
% savename = 'Source_PS_az_tot_9reg_k'; % no low 4 reg (2 main & 2 contrasts with control decoding)


savename2 = 'Decod_6reg_Source_PS_az_tot_9reg_k50';
% savename2 = ['Decod_6reg_' savename '_k' num2str(kclust)];
% savename2 = ['Decod_6reg_' savename num2str(kclust)];
% savename2 = 'Decod_6reg_Source_PS_az_light_nolow_k6';

% parameters
Headtype = 'volume'; % volume|surface
isgroup  = 1;   % 1 = create avg data file
statkeep = 1;   % 1 = keep individual data for later statistics
avg_fq   = 1;   % avg by fq bands


% Loop on conditions
file = [];
savename3 = [savename2 ' (' num2str(length(SUJETS)) 's)' ];

% B/ Find filename: timefreq data
out = [];
for isuj = 1:length(SUJETS)
    subject = ['groove_' SUJETS{isuj}];
    [sStudy, ~] = bst_get('Study',file_fullpath( ...
        [subject '/@intra/brainstormstudy.mat']));
    varok = 0;
    for x = 1:length(sStudy.Timefreq)
        if strcmp(sStudy.Timefreq(x).Comment, savename2) == 1
            % display([subject ' data ok']);
            varok = 1;
            break
        end
    end
    if varok == 0; display('filename problem'); pause; end
    out = [out, {sStudy.Timefreq(x).FileName}];
end
if numel(out) ~= numel(SUJETS); display('filename problem'); pause; end


% Process: Project on default anatomy
out3 = bst_process('CallProcess', 'process_project_sources', out, [], ...
    'headmodeltype', Headtype);

% Process: Group in time or frequency bands
if avg_fq == 1
%     out3 = bst_process('CallProcess', 'process_tf_bands', out3, [], ...
%         'isfreqbands', 1, ...
%         'freqbands',  {'delta', '1, 3', 'mean'; 'theta', '4, 9','mean';...
%                        'alpha', '9, 11', 'mean'; 'beta1', '15, 22', 'mean'; ...
%                        'beta2', '24, 27', 'mean'; 'gamma', '35, 45', 'mean'; 'high-gamma', '60, 80', 'mean'}, ...
%         'istimebands', 0, 'timebands', '', 'overwrite', 1);
%    savename3 = [savename2 '_7bands_(3) (' num2str(length(SUJETS)) 's)' ];

    out3 = bst_process('CallProcess', 'process_tf_bands', out3, [], ...
        'isfreqbands', 1, ...
        'freqbands',  {'low-delta', '1.2, 1.4', 'mean'; '2Hz', '1.8, 2.3', 'mean'; '4Hz', '3.8, 4.2','mean';...
                       'theta', '5, 10', 'mean'; 'alpha_low-beta', '10, 17', 'mean'; 'beta', '20, 30', 'mean'}, ...
        'istimebands', 0, 'timebands', '', 'overwrite', 1);
    
    savename3 = [savename2 '_6bands_chanFDR005_ (' num2str(length(SUJETS)) 's)' ];
    
end
% if avg_fq == 1
%     out3 = bst_process('CallProcess', 'process_tf_bands', out3, [], ...
%         'isfreqbands', 1, ...
%         'freqbands',  {'delta', '1, 4', 'mean'; 'theta', '4, 8', 'mean'; 'alpha', '8, 13', 'mean'; 'beta', '13, 30', 'mean'; ...
%                        'gamma', '30, 50', 'mean'; 'high-gamma', '50, 100', 'mean'}, ...
%         'istimebands', 0, 'timebands', '', 'overwrite', 1);
%     savename3 = [savename2 '_6bands (' num2str(length(SUJETS)) 's)' ];
% end

if isgroup == 1
    % Process: Average: Everything
    out4 = bst_process('CallProcess', 'process_average', out3, [], ...
        'avgtype', 1, 'avg_func', 1, 'keepevents', 0); % Arithmetic avg
    
    % Process: Set comment (rename)
    % x = str_split(savename, ' ');
    bst_process('CallProcess', 'process_set_comment', ...
        out4, [], 'tag', savename3, 'isindex', 1); % [Condname{icond} x{2}]
end

% keep - delete
if statkeep == 1; file = out3;
else; bst_process('CallProcess', 'process_delete', out3, [], 'target', 1);
end
clear out2 out3 out4 sStudy x

% step2: statistics: one sample
% Process: t-test paired
bst_process('CallProcess', 'process_test_parametric1', file, [],...
    'timewindow', [], 'scoutsel', {}, 'scoutfunc', 1, ...
    'isnorm', 0, 'avgtime', 0, 'Comment', ['param ttest ' savename3], ...
    'test_type', 'ttest_onesample', 'tail', 'two'); % One-sample Student's t-test

% step3: delete individual files in Group_analyses
bst_process('CallProcess', 'process_delete', file, [], 'target', 1);
clear out file

%% Import and compute Blob topography for decoding in Bst


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Build figure of decoding for each reg by band
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Script generated by Brainstorm (03-Feb-2020)
reg = 2; % 1-> periodic | 2-> syncope | 3 -> Groove

if reg == 2
        name_reg = 'syncope';
elseif reg == 3
        name_reg = 'groove';
end


fth = logspace(0,2,100);
name_freq = {'low-delta' '2Hz' '4Hz' 'theta' 'alpha_low-beta' 'beta',};
% Input files
% % **param ttest Decod_6reg_Source_PS_az_tot_9reg_k50_6bands_chanFDR005_ (30s)_02 | alpha=0.001 (FDR:1,3)
% sFiles = {'Group_analysis/@intra/timefreq_fft_pthresh_210312_1253.mat'};

% |||||param ttest Decod_6reg_Source_PS_az_tot_9reg_k50_6bands_chanFDR005_ (30s) | alpha=0.005 (FDR:1)
sFiles = {'Group_analysis/@intra/timefreq_fft_pthresh_210811_1330.mat'};
RawFiles = {'C:\Users\zalta\Desktop\Blob_Groove\blob reg by freq band'};

% Start a new report
bst_report('Start', sFiles);

for i0 = 1:6
    % Process: Export to SPM8/SPM12 (volume)
    sFiles = bst_process('CallProcess', 'process_export_spmvol', sFiles, [], ...
        'outputdir',      {RawFiles{1}, 'Nifti1'}, ...
        'filetag',        name_reg, ...
        'isconcat',       0, ...
        'timewindow',     [reg, reg], ...
        'timedownsample', 0, ...
        'timemethod',     2, ...  % Keep time dimension (4D volume)
        'freq_export',    i0, ...
        'voldownsample',  0, ...
        'iscut',          0);
    
    % Process: Export to SPM8/SPM12 (volume)
    sFiles = bst_process('CallProcess', 'process_export_spmvol', sFiles, [], ...
        'outputdir',      {RawFiles{1}, 'Nifti1'}, ...
        'filetag',        name_reg, ...
        'isconcat',       0, ...
        'timewindow',     [reg, reg], ...
        'timedownsample', 0, ...
        'timemethod',     2, ...  % Keep time dimension (4D volume)
        'freq_export',    i0, ...
        'voldownsample',  0, ...
        'iscut',          1);
    
    
    % show size of the file for selection
    s = dir([RawFiles{1} '\' name_reg '_' name_freq{i0} '_02.nii']);
    filesize = s.bytes;
    
    if filesize > 1e4
        delete([RawFiles{1} '\' name_reg '_' name_freq{i0} '_02.nii']);
    else
        delete([RawFiles{1} '\' name_reg '_' name_freq{i0} '.nii']);
        delete([RawFiles{1} '\' name_reg '_' name_freq{i0} '_02.nii']);
    end

end

% % % % Then import nii files as surfaces in bst 
% % % % Then merge the blob by frequency bands
% % % 
% % % right click in the font of the figure -> "Figure" -> "Change Backgroung color" -> "RVB"
% % % Image Font color -> [0 0 0] = Black
% % % 
% % % cortex_15002V
% % %     Sulci on
% % %     color        -> [240 240 240]
% % %     Transparency -> 20%
% % %     Smooth       -> 10%
% % %     
% % % Syncop
% % %     color        -> [75 75 204]
% % %     Transparency -> 20%
% % %     Smooth       -> 40%
% % %     
% % % Groove
% % %     color        -> [204 75 75]
% % %     Transparency -> 20%
% % %     Smooth       -> 40%
% % %     
% % % Orientation
% % %     Left         -> view(172,5)
% % %     Right        -> view(8,5)
% % %     Add internal views -> remove the controlateral hemisphere and corresponding blobs
% % %     
% % %     
% % % Orientation for each scout
% % % 
% % %     'lPAC'      -> view(172,5)
% % %     'rPAC'      -> view(8,5)
% % %     'SMA'       -> view(270,90)
% % %     'rMotor'    -> view(8,5)
% % %     'lParietal' -> view(172,5)
% % % 
% % % saveas(gcf,'C:\Users\zalta\Desktop\Blob_Groove\blob reg by freq band\images blob\lPAC.png')