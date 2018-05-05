%% navigate to the folder and download DS
addpath(genpath('~/Develop/matlab'));
addpath(genpath('~/Google Drive/'));
cd ~/Develop/Mendoza__ReplicationEvolution/Data/
load('MutantChr200.mat');
DS = G;
%% Fig 2A: example of chromosome 5, left arm - % unreplicated against VS distance to the end, kbp
chr_num = 5; K = 75; unq_mutant = {'cdc20' , 'dbf2'};
clrs1 = summer(12); clrs2 = hot(12); clrs_set = [clrs1(3,:) ; clrs2(3,:)];
figure('units','centimeters','position',[5 5 8 8]);
hold on; grid on; set(gca , 'FontSize' , 12);
for I = 1:length(unq_mutant)
	idx = find( strcmp(DS.MutantID , unq_mutant{I}) & DS.chr_num == chr_num & DS.rptbl_flg == 1);
	data = NaN( length(idx) , K);
    for J = 1:length(idx)
        y = DS.percent_underreplicated{idx(J)}*100;
        start_point_kb = DS.start_point_kb{idx(J)};
        for Z = 1:K
            idx_current_kb = find(start_point_kb == Z);
            data(J , Z) = nanmedian(y(idx_current_kb));
        end
    end

	x = [1:K];
	mean_y = nanmedian(data); std_y = nanstd(data);
	plot(x , mean_y , 'LineWidth' , 3 , 'color' , clrs_set(I,:) ,...
            'Display' , strcat( 'mean,' , unq_mutant{I} ) );
    h = fill([x';flipud(x')],[mean_y'-std_y';flipud(mean_y'+std_y')], clrs_set(I,:) ,...
            'linestyle','none' ,  'Display' , strcat( '2xstd,' , unq_mutant{I} ) );
    set(h,'facealpha',.3);
end
xlim([0 K]);
ylim([0 60]);
set(gca , 'Xtick' , [0 25 50 75]);
set(gca , 'Ytick' , [0 20 40 60]);
title(strcat ( DS.chr{idx(1)} , ', left'));
%set(gcf , 'PaperPosition' , [0 0 20 10]);
xlabel('Distance to the end, kbp');
ylabel('% unreplicated cells');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2A' , '-r300');

%%
load('DS_Tsveti.mat');
DS = G;
chr_num = 5; K = 75; unq_mutant = {'M' , 'S'};
clrs1 = summer(12); clrs2 = hot(12); clrs_set = [clrs1(3,:) ; clrs2(3,:)];
figure('units','centimeters','position',[5 5 8 8]);
hold on; grid on; set(gca , 'FontSize' , 12);
for I = 1:length(unq_mutant)
	idx = find( strcmp(DS.MutantID , unq_mutant{I}) & DS.chr_num == chr_num);
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
	mean_y = nanmedian(data); std_y = nanstd(data);
	plot(x , mean_y , 'LineWidth' , 3 , 'color' , clrs_set(I,:) ,...
            'Display' , strcat( 'mean,' , unq_mutant{I} ) );
    h = fill([x';flipud(x')],[mean_y'-std_y';flipud(mean_y'+std_y')], clrs_set(I,:) ,...
            'linestyle','none' ,  'Display' , strcat( '2xstd,' , unq_mutant{I} ) );
    set(h,'facealpha',.3);
end
xlim([0 K]);
ylim([-10 60]);
set(gca , 'Xtick' , [0 25 50 75]);
set(gca , 'Ytick' , [0 20 40 60]);
title(strcat ( DS.chr{idx(1)} , ', left'));
%set(gcf , 'PaperPosition' , [0 0 20 10]);
xlabel('Distance to the end, kbp');
ylabel('% unreplicated cells');

%%
load('DS_stat__200bp.mat');
DS = sortrows(DS , {'chr_num' , 'start_point'});
load('DS_Tsveti.mat');
DS.percent_underrep_Tsveti_M = NaN(length(DS) , 1);
DS.percent_underrep_Tsveti_S = NaN(length(DS) , 1);
for I = 1:length(DS)
    
    idx_M = find(DS.chr_num(I) == G.chr_num & strcmp(G.MutantID , 'M'));
    x = G.start_point{idx_M(1)};
    idx = find(x == DS.start_point(I));
    data = NaN(length(idx_M) , 1);
    for J = 1:length(idx_M)
        y = G.unreplicated_cells_not_trimmed{idx_M(J)};
        data(J) = y(idx);
    end
    DS.percent_underrep_Tsveti_M(I) = nanmedian(data);
        
    idx_S = find(DS.chr_num(I) == G.chr_num & strcmp(G.MutantID , 'S'));
    x = G.start_point{idx_S(1)};
    idx = find(x == DS.start_point(I));
    data = NaN(length(idx_S) , 1);
    for J = 1:length(idx_S)
        y = G.unreplicated_cells_not_trimmed{idx_S(J)};
        data(J) = y(idx);
    end
    DS.percent_underrep_Tsveti_S(I) = nanmedian(data);   
end
        
%%
DS = sortrows(DS , 'percent_underrep_Tsveti_M');

clrs = parula(12); clrs2 = hot(12);
idx = find(DS.percent_underrep_Tsveti_M >= 0);
D = DS(idx , :);
figure('units','centimeters','position',[5 5 8 8]); hold on; grid on; set(gca , 'FontSize' , 10);
scatter(D.percent_underrep_Tsveti_M , D.percent_underrep_Tsveti_S , 30 , clrs(5,:) , 'filled' , ...
    'MarkerFaceAlpha' , .2);
%S = CalcDM_Newman06( D.percent_underrep_Tsveti_M , D.percent_underrep_Tsveti_S  , 200, 0 );
%plot(D.percent_underrep_Tsveti_M , S.DM_ypred , 'LineWidth' , 2 , 'color' , clrs2(3,:));
xlabel('Metaphase');
ylabel('Metaphase + CDK inhibition');
R = corrcoef(D.percent_underrep_Tsveti_M , D.percent_underrep_Tsveti_S );
[~,p] = ttest2(D.percent_underrep_Tsveti_M , D.percent_underrep_Tsveti_S );
title(sprintf('R = %.2f. p < 0.001.' , R(1,2)));
plot([0 50 ] , [0 50] , 'LineWidth' , 2 , 'LineStyle' , '--' , 'color' , [.2 .2 .2]);
xlim([0 50]);
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig4/4A' , '-r300');
%% Fig 2BC -- scatter: x-axis - length underreplication, y-axis - % unreplicated cells, one dot -- one chromosomal end
figure; hold on; grid on; set(gca , 'FontSize' , 10);
markers_left = {'d' , 'o' , 's'}; markers_right = {'x' , '+' , '*'};
clrs1 = hot(12); clrs2 = lines(6); clrs3 = parula(12); clrs4 = summer(12); clrs5 = spring(12);
clrs_set = [clrs1(3,:) ; clrs2(4,:) ; clrs3(2,:) ; clrs3(10,:) ; clrs4(2,:) ; clrs5(4,:) ; clrs3(8,:)];
unq_mutant = {'cdc20' }; chr_nums = [1:16];
for I1 = 1:length(unq_mutant)
    for I2 = 1:length(chr_nums)
        idx = find(strcmp(DS.MutantID , unq_mutant{I1}) & DS.chr_num == chr_nums(I2));
        scatter(nanmedian(DS.length_underreplicated_left(idx)) , ...
            nanmedian(DS.percent_underreplicated_left(idx))*100 , ...
            100 , clrs_set(rem(I2,7)+1,:),'filled' , markers_left{rem(I2,3)+1} , 'Display' , strcat(DS.chr{idx(1)} , '-left'));
        scatter(nanmedian(DS.length_underreplicated_right(idx)) , ...
            nanmedian(DS.percent_underreplicated_right(idx))*100 , ...
            100 , clrs_set(rem(I2,7)+1,:) , markers_right{rem(I2,3)+1} , 'LineWidth' , 3 , 'Display' , strcat(DS.chr{idx(1)} , '-right'));
    end
end
xlabel('Length of under-replicated subtelomeric region, kbp');
ylabel('% unreplicated cells at the last kbp');
legend('location' , 'EastOutside');
set(gcf , 'PaperPosition' , [0 0 20 20]);
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2BC' , '-r300');

%% Fig 2B -- hist option (only cdc20)
figure('units','centimeters','position',[5 5 8 8]); hold on; grid on; set(gca , 'FontSize' , 10);
clrs1 = summer(12); 
data = []; 
unq_mutant = {'cdc20'}; chr_nums = [1:16];
for I1 = 1:length(unq_mutant)
    for I2 = 1:length(chr_nums)
        idx = find(strcmp(DS.MutantID , unq_mutant{I1}) & DS.chr_num == chr_nums(I2));
        data = [data nanmedian(DS.length_underreplicated_left(idx)) nanmedian(DS.length_underreplicated_right(idx)) ];
    end
end
h = histogram(data , [0:5:40] );
%h.EdgeColor = 'w';
set(h , 'EdgeColor' , 'w' , 'FaceColor' , clrs1(3,:));
xlabel('Length of under-replicated subtelomeric region, kbp');
ylabel('Number of chromosomal ends');
%set(gcf , 'PaperPosition' , [0 0 20 20]);
title(sprintf('Mean = %.2f.' , nanmean(data)));
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2B_hist' , '-r300');

%% Fig 2C -- hist option (only cdc20)
figure('units','centimeters','position',[5 5 8 8]); hold on; grid on; set(gca , 'FontSize' , 10);
clrs1 = summer(12); 
data = []; 
unq_mutant = {'cdc20'}; chr_nums = [1:16];
for I1 = 1:length(unq_mutant)
    for I2 = 1:length(chr_nums)
        idx = find(strcmp(DS.MutantID , unq_mutant{I1}) & DS.chr_num == chr_nums(I2));
        data = [data nanmedian(DS.percent_underreplicated_left(idx))*100 nanmedian(DS.percent_underreplicated_right(idx))*100 ];
    end
end
h = histogram(data , [40:5:70] );
set(h , 'EdgeColor' , 'w' , 'FaceColor' , clrs1(3,:));
xlabel('% unreplicated cells at the last kbp');
ylabel('Number of chromosomal ends');
%set(gcf , 'PaperPosition' , [0 0 20 20]);
title(sprintf('Mean = %.2f.' , nanmean(data)));
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2C_hist' , '-r300');


%% Fig 2B -- bar option w/ errors (only cdc20)
figure; hold on; grid on; set(gca , 'FontSize' , 10);
clrs1 = parula(12); clrs2 = lines(6);
data_median = []; 
data_std = [];
unq_mutant = {'cdc20'}; chr_nums = [1:16];
for I1 = 1:length(unq_mutant)
    for I2 = 1:length(chr_nums)
        idx = find(strcmp(DS.MutantID , unq_mutant{I1}) & DS.chr_num == chr_nums(I2));
        data_median = [data_median nanmedian(DS.length_underreplicated_left(idx)) nanmedian(DS.length_underreplicated_right(idx)) ];
        data_std = [data_std nanstd(DS.length_underreplicated_left(idx)) nanstd(DS.length_underreplicated_right(idx))];
    end
end
for I = 1:length(data_median)/2
    bar(2*I-1, data_median(2*I-1) , 'FaceColor' , clrs1(10,:) , 'EdgeColor' , 'w');
    bar(2*I, data_median(2*I) , 'FaceColor' , clrs2(4,:) , 'EdgeColor' , 'w');
end
% for I = 1:length(data_median)/2
%     errorbar(2*I-1, data_median(2*I-1) , data_std(2*I-1) , 'LineWidth' , 2 , 'color' , [.2 .2 .2]);
%     errorbar(2*I, data_median(2*I) , data_std(2*I) , 'LineWidth' , 2 , 'color' , [.2 .2 .2]);
% end
xlabel('Chromosome');
ylabel('Length of under-replicated subtelomeric region, kbp');
set(gca , 'Xtick' , [1.5:2:31.5] , 'XtickLabel' , {'I' , 'II' , 'III' , 'IV' , ...
    'V' , 'VI' , 'VII' , 'VIII' , ...
    'IX' , 'X' , 'XI' , 'XII' , ...
    'XIII' , 'XIV' , 'XV' , 'XVI'} , 'FontSize' , 8);
legend('Left arm' , 'Right arm' , 'location' , 'se');
set(gcf , 'PaperPosition' , [0 0 10 8]);
title(sprintf('Mean = %.2f.' , nanmean(data_median)));
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2B_bar' , '-r300');

%% Fig 2C -- bar option w/ errors (only cdc20)
figure; hold on; grid on; set(gca , 'FontSize' , 10);
clrs1 = parula(12); clrs2 = lines(6);
data_median = []; 
data_std = [];
unq_mutant = {'cdc20'}; chr_nums = [1:16];
for I1 = 1:length(unq_mutant)
    for I2 = 1:length(chr_nums)
        idx = find(strcmp(DS.MutantID , unq_mutant{I1}) & DS.chr_num == chr_nums(I2));
        data_median = [data_median nanmedian(DS.percent_underreplicated_left(idx)) nanmedian(DS.percent_underreplicated_right(idx)) ];
        data_std = [data_std nanstd(DS.percent_underreplicated_left(idx)) nanstd(DS.percent_underreplicated_right(idx))];
    end
end
data_median = data_median*100;
for I = 1:length(data_median)/2
    bar(2*I-1, data_median(2*I-1) , 'FaceColor' , clrs1(10,:) , 'EdgeColor' , 'w');
    bar(2*I, data_median(2*I) , 'FaceColor' , clrs2(4,:) , 'EdgeColor' , 'w');
end
% for I = 1:length(data_median)/2
%     errorbar(2*I-1, data_median(2*I-1) , data_std(2*I-1) , 'LineWidth' , 2 , 'color' , [.2 .2 .2]);
%     errorbar(2*I, data_median(2*I) , data_std(2*I) , 'LineWidth' , 2 , 'color' , [.2 .2 .2]);
% end
xlabel('Chromosome');
ylabel('% unreplicated cells at the last kbp');
set(gca , 'Xtick' , [1.5:2:31.5] , 'XtickLabel' , {'I' , 'II' , 'III' , 'IV' , ...
    'V' , 'VI' , 'VII' , 'VIII' , ...
    'IX' , 'X' , 'XI' , 'XII' , ...
    'XIII' , 'XIV' , 'XV' , 'XVI'} , 'FontSize' , 8);
legend('Left arm' , 'Right arm' , 'location' , 'ne');
set(gcf , 'PaperPosition' , [0 0 10 8]);
title(sprintf('Mean = %.2f.' , nanmean(data_median)));
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2C_bar' , '-r300');




%% Fig Sup 2A_1: generate a figure, for each chromosome/arm - a fraction of unreplicated (y axis) VS distance to the end (x axis)
% green line - pMet3cdc20; red line - dbf2-2
chr_nums = [1:16];
% dist to the end of the chromosome, KB
K = 30;

unq_mutant = {'cdc20' , 'dbf2'};
clrs1 = summer(12);
clrs2 = hot(12);
clrs_set = [clrs1(3,:) ; clrs2(3,:)];
figure; 
for L = 1:length(chr_nums)
    subplot(8,4,2*L-1); hold on; grid on; set(gca , 'FontSize' , 9);
    for I = 1:length(unq_mutant)
        idx = find( strcmp(DS.MutantID , unq_mutant{I}) & DS.chr_num == chr_nums(L) & DS.rptbl_flg == 1);
        data = NaN( length(idx) , K);
        for J = 1:length(idx)
            y = DS.percent_underreplicated{idx(J)}*100;
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
    ylim([0 60]);
    set(gca , 'Xtick' , [0 10 20 30]);
    set(gca , 'Ytick' , [0 20 40 60]);
    title(strcat ( DS.chr{idx(1)} , ', left'));
end

% each chromosome separetely, Right arm, dbf2&cdc20 mutants together

%
for L = 1:length(chr_nums)
    subplot(8,4,2*L); hold on; grid on; set(gca , 'FontSize' , 9);
    for I = 1:2
        idx = find( strcmp(DS.MutantID , unq_mutant{I}) & DS.chr_num == chr_nums(L) & DS.rptbl_flg == 1);
        data = NaN( length(idx) , K);
        for J = 1:length(idx)
            y = DS.percent_underreplicated{idx(J)}*100;
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
    ylim([0 60]);
    set(gca , 'Xtick' , [0 10 20 30]);
    set(gca , 'Ytick' , [0 20 40 60]);
    title(strcat ( DS.chr{idx(1)} , ', right'));
end

set(gcf , 'PaperPosition' , [0 0 20 30]);
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/Sup2A_1' , '-r300');

%% Fig Sup 2A_2: each rep separately, all gDNA (only cdc20)
unq_mutant = {'cdc20'}; chr_nums = [1:16];
clrs1 = hot(12); clrs2 = lines(6); clrs3 = parula(12); clrs4 = summer(12); clrs5 = spring(12);
clrs_set = [clrs1(3,:) ; clrs2(4,:) ; clrs3(2,:) ; clrs3(10,:) ; clrs4(2,:) ; clrs3(6,:)];
for I = 1:length(unq_mutant)
    figure('units','centimeters','position',[5 5 8 8]);
    for J = 1:length(chr_nums)
        subplot(4,4,J); hold on; grid on; set(gca , 'FontSize' , 6);
        idx = find( strcmp(DS.MutantID , unq_mutant{I}) & DS.chr_num == chr_nums(J) & DS.rptbl_flg == 1);
        for Z = 1:length(idx)
            x = DS.start_point{idx(Z)}/1000;
            y = DS.percent_underreplicated{idx(Z)}*100;
            plot(x,y,'color' , clrs_set(Z,:));
        end
        title(DS.chr{idx(1)});
        ylim([-10 70]);
        set(gca , 'Ytick' , [0 25 50]);
    end
end

print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/Sup2A_2' , '-r300');











