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
    fig.Position = fig.Position./sqrt(4);

%     % image of the scout blob
%     img = imread(['./Data/Blob_Groove/blob PAC/images blob\ROIs/' scout_Name{i0} '.png']);
%     image(img, 'XData', [10 45], 'YData', [.026 .01])
    
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

    helper_ILL(gcf, 1.1, 5);
    print('./Figures/Fig5e','-painters', '-dpdf');

%     saveas(gcf, ['C:\Users\zalta\Dropbox\Groove\6_Arnaud\Groove_au_propre\Figures\Fig5_PAC_' scout_Name{i0} '_under_45.png'])
end