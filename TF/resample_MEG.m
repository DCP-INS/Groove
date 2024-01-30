%% Downsample Groove data (task and resting)

clearvars
close all
clc

% set path of raw MEG data
origDir = 'xxxx/GrooveS/';

SUJETS  = {'002','003','004','005','006','007','008','009','010','011','012'...
           '014','015','016','017','018','019','020','021','022','023','024','025' ...
           '026','027','028','029','030','031','032'}; % (001&013out)
  
RUNS    = {'TRIGGER_16', 'TRIGGER_32', 'TRIGGER_64'};
Rname   = {'Low', 'Medium', 'High'};

% Select task (istask = true) or resting data (istask = false)
istask = true; 

% computation parameters
nchan   = 299;
ntrial  = 48;
fs = 508.6275; % resampling
badtrial = ones(length(SUJETS),length(RUNS),ntrial);
badchan = ones(length(SUJETS),nchan);
s = [];

for isuj = 1:length(SUJETS)
    
    if logical(istask)
        subject = sprintf('groove_%s',SUJETS{isuj});
    else
        subject = sprintf('rest_groove_%s',SUJETS{isuj});
    end
    
    for idata = 1:length(RUNS)
        run = RUNS{idata};

        % Process: Select files
        files = bst_process('CallProcess', 'process_select_files_data', [], [], ...
            'subjectname', subject, 'condition', run, ...
            'includebad', 1, 'includeintra', 0, 'includecommon', 0);

       % Control number of trials
        if idata == 4; files = files(1:ntrial); 
        elseif length(files) ~=ntrial; fprintf('error n trial files\n'); pause; end

         for i0 = 1:ntrial       
            try

            % load channel data
            Data = in_bst_data(files(i0).FileName);
            Data.ChannelFlag = badchan;

            % Process: Interpolate bad electrodes
        %     sfiles = files(i0).FileName;
            Data = bst_process('CallProcess', 'process_eeg_interpbad', Data, [], ...
            'maxdist',     5, ...
            'sensortypes', 'MEG', ...
            'overwrite',   0);


            % select bad trials : if good -> 1 ; if bad -> 0
            if isempty(find(strcmp({Data.Events.label}, 'BAD')==1,1))  
            else
                badtrial(isuj,idata,i0) = 0;
            end

            % select badchan
            Data.ChannelFlag(Data.ChannelFlag==-1) = 0;
            badchan(isuj,:) = badchan(isuj,:).*Data.ChannelFlag';

            % resample channel data      
            [Data, t] = process_resample('Compute', Data.F, Data.Time, fs);
            s(idata, i0,:,:) = Data;
            if ~exist('tim', 'var')
                tim = t;
            elseif ~eq(mean(tim), mean(t))
                sprintf('error'); pause;
            end

            catch ME
                if isempty (regexp(ME.message,'Cannot load recordings file','once')) == 1
                else
                    fprintf([subject '_' run '_' num2str(i0) '\n']);
                    s(idata, i0,:,:) = s(idata, i0-1,:,:);
                end
            end

            clear Data t tf deb fin ME
        end
    end    
    
    % save data
    save([ 'XXXX\GrooveS_downsample\signal_brut_downsample_' subject '.mat'],'s','tim','-v7.3');
        clear s tim files run 
end

badchan(badchan==0) = -1;        
save('XXX/GrooveS_downsample/bad.mat','badtrial','badchan','iMeg');
clear badchan badtrial
    