clearvars
close all
clc

cd('C:\Users\arnau\Desktop')

% Plot decoding scouts
load('.\Groove_final\Data\scout_Frites_Atlas_Contrast.mat');   % 5 scouts
load('.\Groove_final\Data\Decoding_Scouts_Ps_az_Frites_Atlas_Contrast.mat')

OUT(10,:,:,:) = []; % remove 011 if not removed during computation
size(OUT) % must be 29 participants

% Names of scouts
titles = {'right auditory cortex' 'left auditory cortex'...
          'supplementary motor area' 'right motor cortex' 'left parietal cortex'};

for i0 = 1:length(Scouts) 
    figure; set(gcf,'color','w'); colormap([.5 .5 .5; 1 1 1]); hold on
    
    subplot(2,2,1)
    
    s1 = squeeze(OUT(:,i0,1,:)); % Periodic 
    s2 = squeeze(OUT(:,i0,2,:)); % syncop
    s3 = squeeze(OUT(:,i0,3,:)); % Groove

    %Stats
    pval = 0.005;
    for z0 = 1:length(fth)
        [h, p1] = ttest(squeeze(s1(:,z0)) ,0, pval, 'both');
        p_tot1(z0) = p1;

        [h, p2] = ttest(squeeze(s2(:,z0)) ,0, pval, 'both');
        p_tot2(z0) = p2; 

        [h, p3] = ttest(squeeze(s3(:,z0)),0, pval,  'both');
        p_tot3(z0) = p3; 
        clear h p
    end

    [h] = mafdr(p_tot1);%, 'BHFDR', 'true');
    z1 = NaN(1, length(fth)); z1(h<pval) = 0; fth(z1 == 0);

    [h] = mafdr(p_tot2);%, 'BHFDR', 'true');
    z2 = NaN(1, length(fth)); z2(h<pval) = 0; fth(z2 == 0);

    [h] = mafdr(p_tot3);%, 'BHFDR', 'true');
    z3 = NaN(1, length(fth)); z3(h<pval) = 0; fth(z3 == 0);

% image of the scout blob
    hold on
    img = imread(['.\Groove_final\Data\Blob_Groove\blob scouts\Contrast_Frites_Atlas\Image scouts\' Scouts(i0).Label '.png']);
    image(img, 'XData', [5 100], 'YData', [.4 .20])  %  'XData', [10 100], 'YData', [.35 .20])

    %figure
    Ip = 100;
    xn = linspace(1, length(fth), length(fth) * Ip);

    % yn = spline(fth, s1, xn);
    %     plot(xn, mean(yn), 'Color',[.3 .3 .8], 'Linewidth',2) % color in RGB -> [76.5 76.5 204]
    %     boundedline(xn, mean(yn), std(yn)/sqrt(size(yn,1)), ...
    %         'b-', 'transparency', 'alpha', .25);
    % p =  plot(fth, z1-0.1, 'Color', [.3 .3 .8], 'Linewidth', 3, 'DisplayName','Periodic');


    yn = spline(fth, s2, xn);
        plot(xn, mean(yn), 'Color',[.3 .3 .8], 'Linewidth',2)
        boundedline(xn, mean(yn), std(yn)/sqrt(size(yn,1)), ...
            'cmap',[.3 .3 .8], 'transparency', 'alpha', .25);
    s = plot(fth, z2-0.08, 'Color', [.3 .3 .8], 'Linewidth', 2.5, 'DisplayName','Syncope');
        plot(fth, z2-0.08, ':.', 'Color', [.3 .3 .8], 'Linewidth', 2.5)


    yn = spline(fth, s3, xn);
        plot(xn, mean(yn), 'Color',[.8 .3 .3], 'Linewidth',2)
        boundedline(xn, mean(yn), std(yn)/sqrt(size(yn,1)), ...
            'cmap',[.8 .3 .3], 'transparency', 'alpha', .25);
    g = plot(fth, z3-0.065, 'Color', [.8 .3 .3], 'Linewidth', 2.5, 'DisplayName','Groove');
        plot(fth, z3-0.065, ':.', 'Color', [.8 .3 .3], 'Linewidth', 2.5)


    % stat contrast
    pval = 0.05; %!
    [~, p] = ttest(s2, s3, pval, 'both');
    [h] = mafdr(p); %,'BHFDR', 'true');
    q = NaN(1, length(fth)); q(h<pval & ( z2==0 | z3==0 ) ) = 0; % contrast significant & 1 main effect significant !

    plot(fth, q-0.05, ':.', 'Color', 'k', 'Linewidth', 2.5)
    plot(fth, q-0.05,       'Color', 'k', 'Linewidth', 2.5)



    yline(0,':');
    set(gca, 'FontSize', 11, 'FontName', 'Calibri')
    nbTicks = 6 ;

    xticks([2 5 10 20 50 100])
    xticklabels([2 5 10 20 50 100])
    yticks([-.1 0 .1 .2 .3  .4])
    
    xlabel( 'Neural frequency (Hz)' )
    ylabel('Coding precision (z)')
    
    set(gca,'Layer','top','Box','off','TickLength',[.001 .001])
    ylim([-.09 .4])
    title(titles{i0}, 'Interpreter' ,'none')
    set(gca,'XScale', 'log')
    
    helper_ILL(gcf, 1.29, 6.3*2);
    print(['C:\Users\arnau\Desktop\Groove_final\Figures\FigS3_' Scouts(i0).Label],'-painters', '-dpdf');

%     saveas(gcf, ['.\Groove_final\Figures\FigS3_' Scouts(i0).Label '.png']); 
end
%     legend([s g], 'EdgeColor', 'w')
    