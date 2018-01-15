%% navigate to the folder and download DS
addpath(genpath('~/Develop/matlab'));
addpath(genpath('~/Google Drive/'));
cd ~/Develop/Mendoza__ReplicationEvolution/Data/
load('DS_Michi.mat');

%%
chr_nums = [1:16];
my_legend = cell(2*length(chr_nums) , 1);
mean_data_dbf2 = NaN(2*length(chr_nums) , 1);
mean_data_cdc20 = NaN(2*length(chr_nums) , 1);
std_data_dbf2 = NaN(2*length(chr_nums) , 1);
std_data_cdc20 = NaN(2*length(chr_nums) , 1);
for L = 1:length(chr_nums)
    idx_dbf2 = find( strcmp(DS.MutantID , 'dbf2') & DS.rptbl_flg == 1 & DS.chr_num == chr_nums(L) );
    idx_cdc20 = find( strcmp(DS.MutantID , 'cdc20') & DS.rptbl_flg == 1 & DS.chr_num == chr_nums(L) );
    
    data_dbf2 = DS.percent_unreplicated_left_not_trimmed(idx_dbf2);
    data_cdc20 = DS.percent_unreplicated_left_not_trimmed(idx_cdc20);
    mean_data_dbf2(2*L-1) = mean(data_dbf2);
    mean_data_cdc20(2*L-1) = mean(data_cdc20);
    std_data_dbf2(2*L-1) = std(data_dbf2);
    std_data_cdc20(2*L-1) = std(data_cdc20);
    my_legend{2*L-1} = strcat(DS.chr{idx_dbf2(1)} , '- left');
    
    data_dbf2 = DS.percent_unreplicated_right_not_trimmed(idx_dbf2);
    data_cdc20 = DS.percent_unreplicated_right_not_trimmed(idx_cdc20);
    mean_data_dbf2(2*L) = mean(data_dbf2);
    mean_data_cdc20(2*L) = mean(data_cdc20);
    std_data_dbf2(2*L) = std(data_dbf2);
    std_data_cdc20(2*L) = std(data_cdc20);
    my_legend{2*L} = strcat(DS.chr{idx_dbf2(1)} , '- right');
end
%% y-axis - cdc20; x-axis - dbf2; every dot - percent of unreplicated at the first 1KB
clrs1 = parula(12);
clrs2 = hot(12);
clrs3 = summer(12);
clrs4 = lines(6);
clrs5 = spring(12);
clrs_set = [clrs1(1,:) ; clrs2(3,:) ; .3 .3 .3 ; clrs3(3,:) ; clrs4(4,:) ; clrs1(10,:) ; clrs1(6,:) ; clrs5(3,:) ; ];
marker_set = {'d' , 'o' , 's' , 'x'};
N = length(chr_nums) * 2;
figure; hold on; grid on; set(gca , 'FontSize' , 16);
for I = 1:N
    clrs_idx = (I-1 - rem(I-1 , 4))/4 + 1;
	marker_idx = rem(I-1,4)+1;
    
    plot([mean_data_dbf2(I) mean_data_dbf2(I)], ...
            [mean_data_cdc20(I)-std_data_cdc20(I)/2 mean_data_cdc20(I)+std_data_cdc20(I)/2], ...
            'LineWidth' , 1 , 'color' , clrs_set(clrs_idx,:) );
        
    plot([mean_data_dbf2(I)-std_data_dbf2(I)/2 mean_data_dbf2(I)+std_data_dbf2(I)/2], ...
            [mean_data_cdc20(I) mean_data_cdc20(I)], ...
            'LineWidth' , 1 , 'color' , clrs_set(clrs_idx,:) );
end
for I = 1:N
    clrs_idx = (I-1 - rem(I-1 , 4))/4 + 1;
	marker_idx = rem(I-1,4)+1;
    if round(I/4) == I/4
        scatter(mean_data_dbf2(I) , mean_data_cdc20(I) , 150 , clrs_set(clrs_idx,:) , marker_set{marker_idx} ,...
            'LineWidth' , 4);
        
    else
        scatter(mean_data_dbf2(I) , mean_data_cdc20(I) , 150 , clrs_set(clrs_idx,:) , 'filled' , marker_set{marker_idx} ,...
            'LineWidth' , 4);
    end
end


plot([20 90] , [20 90] , 'LineWidth' , 2 , 'LineStyle' , '--' , 'color' , [.2 .2 .2]);
set(gcf , 'PaperPosition' , [0 0 20 20]);
%print('-dpsc2' , 'fig__scatter__PercentUnderreplication.eps');
%% legend


clrs1 = parula(12);
clrs2 = hot(12);
clrs3 = summer(12);
clrs4 = lines(6);
clrs5 = spring(12);
clrs_set = [clrs1(1,:) ; clrs2(3,:) ; .3 .3 .3 ; clrs3(3,:) ; clrs4(4,:) ; clrs1(10,:) ; clrs1(6,:) ; clrs5(3,:) ; ];
marker_set = {'d' , 'o' , 's' , 'x'};
N = length(chr_nums) * 2;
figure; hold on; grid on; set(gca , 'FontSize' , 16);
for I = 1:N
    clrs_idx = (I-1 - rem(I-1 , 4))/4 + 1;
	marker_idx = rem(I-1,4)+1;
    if round(I/4) == I/4
        scatter(mean_data_dbf2(I) , mean_data_cdc20(I) , 150 , clrs_set(clrs_idx,:) , marker_set{marker_idx} ,...
            'LineWidth' , 4 , 'Display' , my_legend{I});
        
    else
        scatter(mean_data_dbf2(I) , mean_data_cdc20(I) , 150 , clrs_set(clrs_idx,:) , 'filled' , marker_set{marker_idx} ,...
            'LineWidth' , 4 , 'Display' , my_legend{I});
    end
end
legend('Location' , 'EastOutside');

