%% navigate to the folder and download DS
addpath(genpath('~/Develop/matlab'));
addpath(genpath('~/Google Drive/'));
cd ~/Develop/Mendoza__ReplicationEvolution/Data/
load('DS_Michi.mat');

%% generate a figure, for each chromosome/arm - a fraction of unreplicated (y axis) VS distance to the end (x axis)
% green line - pMet3cdc20; red line - dbf2-2
chr_nums = [1:16];
% dist to the end of the chromosome, KB
K = 75;

unq_mutant = {'cdc20' , 'dbf2'};
clrs1 = summer(12);
clrs2 = hot(12);
clrs_set = [clrs1(3,:) ; clrs2(3,:)];
figure; 
for L = 1:length(chr_nums)
    subplot(8,4,2*L-1); hold on; grid on; set(gca , 'FontSize' , 12);
    for I = 1:length(unq_mutant)
        idx = find( strcmp(DS.MutantID , unq_mutant{I}) & DS.chr_num == chr_nums(L) & DS.rptbl_flg == 1);
        data = NaN( length(idx) , K);
        for J = 1:length(idx)
            y = DS.unreplicated_cells_not_trimmed{idx(J)};
            start_point_kb = DS.start_point_kb{idx(J)};
            for Z = 1:K
                idx_current_kb = find(start_point_kb == Z);
                data(J , Z) = nanmedian(y(idx_current_kb));
                
            end
        end

        x = [1:K];
        mean_y = nanmedian(data);
        std_y = nanstd(data);
        plot(x , mean_y , 'LineWidth' , 3 , 'color' , clrs_set(I,:) ,...
            'Display' , strcat( 'mean,' , unq_mutant{I} ) );

        h = fill([x';flipud(x')],[mean_y'-std_y';flipud(mean_y'+std_y')], clrs_set(I,:) ,...
            'linestyle','none' ,  'Display' , strcat( '2xstd,' , unq_mutant{I} ) );
        set(h,'facealpha',.3);
    end
    xlim([0 K]);
    ylim([-5 60]);
    set(gca , 'Xtick' , [0 25 50 75]);
    title(strcat ( DS.chr{idx(1)} , ', left'));
end

% each chromosome separetely, Right arm, dbf2&cdc20 mutants together


for L = 1:length(chr_nums)
    subplot(8,4,2*L); hold on; grid on; set(gca , 'FontSize' , 12);
    for I = 1:2
        idx = find( strcmp(DS.MutantID , unq_mutant{I}) & DS.chr_num == chr_nums(L) & DS.rptbl_flg == 1);
        data = NaN( length(idx) , K);
        for J = 1:length(idx)
            y = DS.unreplicated_cells_not_trimmed{idx(J)};
            start_point_kb = DS.start_point_kb{idx(J)};
            for Z = 1:K
                
                idx_current_kb = find(start_point_kb == nanmax(start_point_kb) - Z + 1);
                data(J , Z) = nanmedian(y(idx_current_kb));
            end
        end

        x = [1:K];
        mean_y = nanmedian(data);
        std_y = nanstd(data);
        plot(x , mean_y , 'LineWidth' , 3 , 'color' , clrs_set(I,:) ,...
            'Display' , strcat( 'mean,' , unq_mutant{I} ) );

        h = fill([x';flipud(x')],[mean_y'-std_y';flipud(mean_y'+std_y')], clrs_set(I,:) ,...
            'linestyle','none' ,  'Display' , strcat( '2xstd,' , unq_mutant{I} ) );
        set(h,'facealpha',.3);
    end
    xlim([0 K]);
    ylim([-5 60]);
    set(gca , 'Xtick' , [0 25 50 75]);
    title(strcat ( DS.chr{idx(1)} , ', right'));
   
end

set(gcf , 'PaperPosition' , [0 0 15 30]);
%print('-dpsc2' , 'cdc20_dbf2__underreplication__VS__dist.eps');

%% colors for the legend (once mutant names are settled):

% cdc20 - summer(12 , 3);
% dbf2 - hot(12 , 3);














