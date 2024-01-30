clear
clc 

X = load('E:\GrooveS\data\Group_analysis\@intra\timefreq_fft_210803_1751.mat');
load('E:\GROOVE_A_CONSERVER\files_for_decoding\scout_Volume_1673__Frites_Atlas.mat');   % scouts commun atlas | 7 scouts

t = [];
t = X.TF(:,:,1);
t = find(t(:,:,:) >0);

% Variable
kclust = 10;
nscout = [2 4];
name = {'rPAC' 'lPAC' 'SMA' 'rMotor' 'lParietal'};
Name = 'Volume 1673: Frites_Atlas_Contrast';

a = []; ind = [];
for ivox = 1:length(t)  
    iloc = X.GridLoc(t(ivox),:);
    idist = sqrt( sum( bsxfun(@minus,X.GridLoc, iloc) .^ 2, 2) );
    [~,id_cur] = mink(idist, kclust); % n closest voxels
    
    Scouts(ivox).Seed = t(ivox);
    Scouts(ivox).Vertices = id_cur;
    a(:, ivox) =  abs(Scouts(ivox).Vertices - Scouts(ivox).Seed); 

end

[~, ind] = mink(a,kclust/2);

% merge similar scouts in one (lPAC and SMA)
for iscout = nscout
Scouts(iscout).Vertices = unique([Scouts(iscout).Vertices(ind(:,iscout)); Scouts(iscout+1).Vertices(ind(:,iscout+1))]); ;
end

Scouts(nscout+1) = [];
[Scouts.Label] = name{:};

save('C:\Users\zalta\Desktop\scout_Frites_Atlas_Contrast.mat', 'Name','Scouts','TessNbVertices');   % scouts commun atlas | 7 scouts


clear
clc
SUJETS  = {'002','003','004','005','006','007','008','009','010','012'...
           '014','015','016','017','018','019','020','021','022','023','024','025' ...
           '026','027','028','029','030','031','032'}; % (001&013out) without 011

 [~, hname] = system('hostname');
if strcmp(hname(1:2),'bn')  == 1
    Dir = '/Users/Bn/';
elseif strcmp(hname(1),'R') == 1
    Dir = 'C:\Users\arnau\';
elseif strcmp(hname(1),'C') == 1
    Dir = 'C:\Users\zalta\';
end         
      
DropDir = [Dir 'Dropbox\Groove\'];

addpath(genpath('/Applications/Work/Softwares/Toolbox/'))
addpath([Dir 'Dropbox\tools\'])
addpath([Dir 'Dropbox\Groove\6_Arnaud\glmfst'])


if     strcmp(hname(1),'R') == 1
    % load MEG data and regressors
    load([DriveName('Maxtor') ':\GrooveS_downsample\MNE_MVPA\Chan_data_7regr.mat'], 'X'); % X = behav; Y = Meg; fth = freqs
    
    % load bad channels & bad trials
    load([DriveName('Maxtor') ':/GrooveS_downsample/bad.mat']) % chan
    load([DriveName('Maxtor') ':/GrooveS_downsample/bad_bn_process.mat']); % trials
    
    % load atlas of scouts
    load([DriveName('Maxtor') ':\Groove_script\scout_Frites_Atlas_Contrast.mat']);
else
    % load MEG data and regressors
    load([Dir 'Documents\GrooveS_downsample\MNE_MVPA\Chan_data_7regr.mat'], 'X'); % X = behav; Y = Meg; fth = freqs
    
    % load bad channels & bad trials
    load([Dir 'Documents/GrooveS_downsample/bad.mat']) % chan
    load([Dir 'Documents/GrooveS_downsample/bad_bn_process.mat']); % trials
    
    % load atlas of scouts
    load([ Dir 'Desktop\scout_Frites_Atlas_Contrast.mat']);   % 5 scouts

end
% parameters
ntrial = 144; nvox = 1673;
fth = logspace(0,2, 100);
OUT = zeros(size(X,1),size({Scouts.Label}, 2), 3, length(fth)); % 3 regressors


%% Decoding
for isuj = 1:size(X,1)
    disp(['participant ' num2str(isuj) '/' num2str(size(X,1))]);
    
%%%%%%%%%%%%%%%%%%%%%POUR RAPPEL%%%%%%%%%%%%%%%%%%%%%%%
%x_cur0 -> low(1) vs. med|high (0) -> Periodic
%x_cur1 -> sync
%x_cur2 -> bhv 
%x_cur6 -> sync orthogonalized with bhv
%x_cur3 -> bhv orthogonalized with sync
%x_cur4 -> sync orthogonalized with periodic (xcur0)
%x_cur5 -> bhv othogonalized with synco orthogonalized with periodic (xcur4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 1. Load TF data (Ps_az)
    load([DriveName('Maxtor') ':\Groove_script\Results\groove_' SUJETS{isuj} '\source_Ps_az.mat']); %!

    y_cur = TF;
    y_cur = permute(TF, [2 3 1]);
    y_cur = reshape(y_cur, 3, nvox, length(fth), []);
    
    % regressors
    x_cur = squeeze(X(isuj,:,:));
    x_cur(:,1) = x_cur(:,1)*-1; %!

    % remove bad trials
    goodtrial = badtrial_bnprocess(isuj,:)';
    goodtrial = find(goodtrial == 1);
    x_cur = x_cur(goodtrial,:); % include/exclude bads
    x_cur = x_cur(:,1:3); % 3 regressors
    
    % multivariate analysis subject by subject
    nfold = 10 ; %size(x_cur,1); % 10
    alpha = 2; % 1 ridge parameter (to adjust?)
    r2z   = @(r)0.5*log((1+r)./(1-r));
    
    for z0 = 1:size({S.Label}, 2)
            disp(['scout ' num2str(z0) '/' num2str(size({Scouts.Label}, 2))]);

        y = y_cur(:,Scouts(z0).Vertices,:,goodtrial); 
        y = reshape(y, 3*size(Scouts(z0).Vertices,1), length(fth), length(goodtrial));
        
        for ifq = 1:length(fth)
            for i0 = 1:size(x_cur,2)
                out_dec = decode_ridg(x_cur(:,i0), squeeze(y(:,ifq,:))', nfold, alpha); % use ridge if only one regressor (lambda or e (3rd parameter) = 1)
                OUT(isuj,z0,i0,ifq) = r2z(out_dec.rdec); % regressors|fq
                clear out_dec
            end
        end
    end
end
clear iStudy R y_cur x_cur r2z
save([Dir 'Dropbox\Groove\6_Arnaud\TheBestOf\Results\Decoding_Scouts_Ps_az_Frites_Atlas_Contrast.mat'], 'OUT', 'fth');

clear nfold alpha
disp('TERMINARES');


%% create blob to show the scouts in the cerebral volume

load('C:\Users\zalta\Desktop\scout_Frites_Atlas_Contrast.mat')


cd C:\Users\zalta\Downloads\brainstorm3
brainstorm

% IN Brainstorm 

% export 'Frites_Atlas_anat_scouts' in matlab as W
W.TF = zeros(1673, 7, 2);

for i0 = 1:5
W.TF(Scouts(i0).Vertices,i0,:) = 1
end

% import W in bst as 'Frites_Atlas_anat_scouts'

% Frites_Atlas_anat_scouts: in bst
sFiles = {'Group_analysis/@intra/timefreq_fft_210726_1550.mat'}; 
RawFiles = {['C:\Users\zalta\Desktop\Blob_Groove\blob scouts\Contrast_Frites_Atlas']};

% Start a new report
bst_report('Start', sFiles);

for i0 = 1:5

name_reg = Scouts(i0).Label;
    % Process: Export to SPM8/SPM12 (volume)
    sFiles = bst_process('CallProcess', 'process_export_spmvol', sFiles, [], ...
        'outputdir',      {RawFiles{1}, 'Nifti1'}, ...
        'filetag',        name_reg, ...
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


% right click in the font of the figure -> "Figure" -> "Change Backgroung color" -> "RVB"
% Image Font color -> [0 0 0] = Black
% 
% cortex_15002V
%     Sulci on
%     color        -> [240 240 240]
%     Transparency -> 20%
%     Smooth       -> 10%
%     
% The scouts
%     color        -> [100 100 100]
%     Transparency -> 20%
%     Smooth       -> 75%
%     
% Orientation for each scout
% 
%     'lPAC'      -> view(172,5)
%     'rPAC'      -> view(8,5)
%     'SMA'       -> view(270,90)
%     'rMotor'    -> view(8,5)
%     'lParietal' -> view(172,5)
