clearvars
close all
clc

cd('C:/Users/arnau/Desktop/Groove_final')

% neural decoding: size-effect and statistics
load('./Data/stimuli/Mvpa_chanPs_accEd.mat') % suj|fq_neural|fq_acc

OUT(10,:,:) = []; % remove outlier participant

fig = figure; set(gcf,'color','w');

for i0 = 2
%     subplot(2,1,i0)
    
    if i0 == 1 % avg suj
        x = squeeze(mean(OUT,1));
        x = x'; 
        s = pcolor(fth, fq_acc, x); % neural in x, acoustic in y
        caxis([0, .2]); 
        h = colorbar;
        h.Label.String = 'Accuracy (r)';
        h.Ticks = [0 .1 .2];
        
    else % stats (FDR-corrected)
        pval = 0.05;
        [~,p] = ttest(OUT,0, pval, 'both');
        p = mafdr(p(:));
        p = reshape(p, numel(fth), numel(fq_acc));
        x = -log10(p);
        x = x'; 
        
        s = imagesc(fth, fq_acc, x); % neural in x, acoustic in y

        caxis([-log10(0.05), -log10(0.001)])
        
         h = colorbar(gca,'Position',...
        [0.78 0.2 0.0078 0.33],...
        'Ticks',[-log10(0.05) -log10(0.01) -log10(0.001)] ,...
        'TickLabels',{'0.05','0.01','0.001'});
    
        h.Label.String = 'q-value'; % FDR corrected
    end

    xlabel('Neural (Hz)')
    ylabel('Acoustic (Hz)')
    
    xticks([2 4 8 16 30 60 100])
    yticks([2 8 16])

    set(gca, 'XScale', 'log');
    set(gca,'Ydir', 'normal', 'Layer','top', 'Box','on', 'TickLength',[.01 .01], 'FontSize', 11, 'FontName', 'Calibri')    

helper_ILL(gcf, 2.4, 3,false);
print('./Figures/Fig1f','-painters', '-dpdf');
end