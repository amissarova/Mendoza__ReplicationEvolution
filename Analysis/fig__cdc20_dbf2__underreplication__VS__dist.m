%% navigate to the folder and download DS
addpath(genpath('~/Develop/matlab'));
addpath(genpath('~/Google Drive/'));
cd ~/Develop/Mendoza__ReplicationEvolution/Data/
load('DS_Michi__200bp.mat');

%% generate a figure, for each chromosome/arm - a fraction of unreplicated (y axis) VS distance to the end (x axis)
% green line - pMet3cdc20; red line - dbf2-2
DS = G;
chr_nums = [1:16];
% dist to the end of the chromosome, KB
K = 33;

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
            y = DS.percent_unreplicated_cells_not_trimmed{idx(J)};
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
            y = DS.percent_unreplicated_cells_not_trimmed{idx(J)};
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

% green line - pMet3cdc20; red line - dbf2-2
DS = G;
chr_nums = [1 11];
% dist to the end of the chromosome, KB
K = 100;

unq_mutant = {'cdc20' , 'dbf2'};
clrs1 = summer(12);
clrs2 = hot(12);
clrs_set = [clrs1(3,:) ; clrs2(3,:)];
figure; 
for L = 1:length(chr_nums)
    subplot(length(chr_nums),2,2*L-1); hold on; grid on; set(gca , 'FontSize' , 11);
    for I = 1:length(unq_mutant)
        idx = find( strcmp(DS.MutantID , unq_mutant{I}) & DS.chr_num == chr_nums(L) & DS.rptbl_flg == 1);
        start_point_kb = DS.start_point_kb{idx(1)};
        start_point_kb = round(start_point_kb);
        N = nanmax(start_point_kb);
        data = NaN( length(idx) , K );
        for J = 1:length(idx)
            y = DS.percent_unreplicated_cells_not_trimmed{idx(J)};
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
    ylim([-20 60]);
    set(gca , 'Xtick' , [0 25 50 75 100]);
    title(strcat ( DS.chr{idx(1)} , ', left'));
    xlabel('Position on chromosome, KB');
    ylabel('Degree of underreplication');
    
    subplot(length(chr_nums),2,2*L); hold on; grid on; set(gca , 'FontSize' , 11);
    for I = 1:length(unq_mutant)
        idx = find( strcmp(DS.MutantID , unq_mutant{I}) & DS.chr_num == chr_nums(L) & DS.rptbl_flg == 1);
        start_point_kb = DS.start_point_kb{idx(1)};
        start_point_kb = round(start_point_kb);
        N = nanmax(start_point_kb);
        data = NaN( length(idx) , K );
        for J = 1:length(idx)
            y = DS.percent_unreplicated_cells_not_trimmed{idx(J)};
            for Z = 1:K
                idx_current_kb = find(start_point_kb == N-Z+1);
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
    ylim([-20 60]);
    set(gca , 'Xtick' , [0 25 50 75 100]);
    title(strcat ( DS.chr{idx(1)} , ', right'));
    xlabel('Position on chromosome, KB');
    ylabel('Degree of underreplication');
end
set(gcf , 'PaperPosition' , [0 0 20 20]);
print('-dpsc2' , 'cdc20_dbf2__underreplication__VS__dist__example.eps');

%%
%load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');
%
idx = find(~isnan(DS.percent_unreplicated_not_trimmed_cdc20) & ...
    DS.percent_unreplicated_not_trimmed_cdc20 < 100000 & ...
    DS.percent_unreplicated_not_trimmed_cdc20 > -100000 & ...
    ~isnan(DS.percent_unreplicated_not_trimmed_dbf2) & ...
    DS.percent_unreplicated_not_trimmed_dbf2 < 100000 & ...
    DS.percent_unreplicated_not_trimmed_dbf2 > -100000 );
D = DS(idx , :);
D = sortrows(D , {'dist_to_the_end_kb' , 'percent_unreplicated_not_trimmed_cdc20'});
figure; hold on; grid on; set(gca , 'FontSize' , 12);
clrs = hot(12);
scatter(D.dist_to_the_end_kb , D.percent_unreplicated_not_trimmed_cdc20 , 20 , [.2 .2 .2] , 'filled');
ylim([-120 100]);
%y = smooth(D.dist_to_the_end_kb , D.percent_unreplicated_not_trimmed_cdc20 , 160 , 'rlowess');
plot(D.dist_to_the_end_kb , y , 'LineWidth' , 2 , 'color' , clrs(3,:));    
xlabel('Distance to the end');
ylabel('Degree of underreplication');
print('-dpsc2' , 'cdc20__all__dist__VS__underreplication_w_spline.eps');

%%
idx = find(~isnan(DS.percent_unreplicated_not_trimmed_cdc20) & ...
    DS.percent_unreplicated_not_trimmed_cdc20 < 100000 & ...
    DS.percent_unreplicated_not_trimmed_cdc20 > -100000);
D = DS(idx , :);
R_coeff = NaN(150 , 3);
for I = 1:150
    idx1 = find(D.dist_to_the_end_kb <= I);
    idx2 = find(D.dist_to_the_end_kb > I);
    C = corrcoef(D.dist_to_the_end_kb(idx1) , D.percent_unreplicated_not_trimmed_cdc20(idx1) );
    R_coeff(I,2) = C(1,2);
    C = corrcoef(D.dist_to_the_end_kb(idx2) , D.percent_unreplicated_not_trimmed_cdc20(idx2) );
    R_coeff(I,3) = C(1,2);
    R_coeff(I,1) = I;
end

%% Get the magnitude of underreplication (degree of underreplication at the last 1kb of underreplication)
%load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_Michi__200bp.mat');
data_cdc20 = NaN(32,1);
data_dbf2 = NaN(32,1);
for I = 1:16
    idx = find( strcmp(G.MutantID , 'dbf2') & G.chr_num == I & G.rptbl_flg == 1);
    
    data = NaN( length(idx) , 1);
	for J = 1:length(idx)
        y = G.percent_unreplicated_cells_not_trimmed{idx(J)};
        data(J,:) = nanmedian(y(1:5));
    end
    data_dbf2(2*I-1) = nanmedian(data);
    
    data = NaN( length(idx) , 1);
	for J = 1:length(idx)
        y = G.percent_unreplicated_cells_not_trimmed{idx(J)};
        data(J,:) = nanmedian(y(end-4:end));
    end
    data_dbf2(2*I) = nanmedian(data);
end


%%
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');

T = grpstats(DS , {'chr', 'chr_num', 'arm'});
T = T(: , {'chr' , 'chr_num' , 'arm'});
T = T(~strcmp(T.arm , 'center' ) , :);

%%
T.length_underreplication = NaN(length(T) , 1);
for I = 1:length(T)
    idx = find(DS.chr_num == T.chr_num(I) & strcmp(DS.arm , T.arm{I}));
	data = DS.percent_unreplicated_not_trimmed_cdc20(idx);
    x = DS.dist_to_the_end_kb(idx);
    p = 0;
	K = 20;
	count = 0;
	N = length(data)+1;
    if strcmp(T.arm{I} , 'left')
        while p < .05
            count = count + 1;
            [~,p] = ttest2(data , data(count:count+K));
        end
        T.length_underreplication(I) = x(count + K/2);
    else
        while p < .05
            count = count + 1;
            [~,p] = ttest2(data , data(N-count-K:N-count));
        end
        T.length_underreplication(I) = x(N - count - K/2);
    end
end
        
    
    
    



























