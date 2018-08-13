%% Pipeline to generate panels for fig2 and potential Supplements for fig 2
addpath(genpath('~/Develop/matlab'));
%cd ~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/
%cd 20170531
%% Fig 2A: example for one chromosome, under-replication VS distance to the end (pMet3-cdc20 and dbf2-2)
load('~/Develop/Mendoza__ReplicationEvolution/Data/MutantChr200_Michi.mat');
DS = G;
chr_num = 5; K = 75; 
unq_mutant = {'cdc20' , 'dbf2'};
%unq_mutant = {'M'};
clrs1 = summer(12); clrs2 = hot(12); clrs_set = [clrs1(3,:) ; clrs2(3,:)];
figure('units','centimeters','position',[5 5 15 10]);
hold on; grid on; set(gca , 'FontSize' , 12);
for I = 1:length(unq_mutant)
	idx = find( strcmp(DS.MutantID , unq_mutant{I}) & DS.chr_num == chr_num );
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
xlim([0 K]); ylim([-9 60]);
set(gca , 'Xtick' , [0 25 50 75] , 'Ytick' , [0 20 40 60]);
title(strcat ( DS.chr{idx(1)} , ', left'));
xlabel('Distance to the end, kbp'); ylabel('% unreplicated cells');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2A' , '-r600');


%% Fig 2C: scatter, all data, dist to the end VS under-rep with highlughting particularly underreplicated
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp_new.mat');
figure('units','centimeters','position',[5 5 15 10]);
hold on; grid on; set(gca , 'FontSize' , 13);
clrs = winter(12);
% all 200-bp windows
h1 = scatter(DS.dist_to_the_end_kb , DS.percent_underreplicated_cdc20*100 , ...
    13 , [.65 .65 .65] , 'filled' , 'MarkerFaceAlpha',.1 , 'Display' , 'all 200-bp windows');
thresh_underrep = 21.0745;
idx = find(DS.percent_underreplicated_cdc20*100 > thresh_underrep);
h2 = scatter(DS.dist_to_the_end_kb(idx) , DS.percent_underreplicated_cdc20(idx)*100 , ...
    20 , clrs(9,:) , 'filled' , 'MarkerFaceAlpha',.4 , 'Display' , 'all 200-bp windows');
ylim([-9 60]);
legend('Whole genome' , 'Under-replicated' , 'location' , 'nw');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2B_scatter' , '-r600');

%% Optional Fig 2C-scatter for DBF2
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp_new.mat');
figure('units','centimeters','position',[5 5 15 10]);
hold on; grid on; set(gca , 'FontSize' , 13);
clrs = winter(12);
% all 200-bp windows
h1 = scatter(DS.dist_to_the_end_kb , DS.percent_underreplicated_dbf2*100 , ...
    13 , [.65 .65 .65] , 'filled' , 'MarkerFaceAlpha',.1 , 'Display' , 'all 200-bp windows');
thresh_underrep = 21.0745;
idx = find(DS.percent_underreplicated_dbf2*100 > thresh_underrep);
h2 = scatter(DS.dist_to_the_end_kb(idx) , DS.percent_underreplicated_dbf2(idx)*100 , ...
    20 , clrs(9,:) , 'filled' , 'MarkerFaceAlpha',.4 , 'Display' , 'all 200-bp windows');
ylim([-25 60]);
legend('whole genome' , 'Under-replicated' , 'location' , 'ne');
%print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/20170531/2B_scatter_dbf2' , '-r300');


%% Fig 2B: hist with all the data and highlighting which are under-replicated
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp_new.mat');
data = DS.percent_underreplicated_cdc20*100; data = data(data > -100); data = sort(data);
s = skewness(data); K = length(data); eps = 0;
while s > eps
    K = K - 1;
    s = skewness(data(1:K));
end
clrs = winter(12);
figure('units','centimeters','position',[5 5 5 7]); hold on; grid on;
[y,x] = ksdensity(data , [-30:1:60]);
thresh_underrep = data(K); 
plot([thresh_underrep thresh_underrep] , [0 .0042] , 'LineWidth' , 4 , 'color' , clrs(9,:));
for I = 1:length(x)
    if x(I) > thresh_underrep
        plot([x(I) x(I)] , [0 y(I)] , 'LineWidth' , 2 , 'color' , clrs(9,:));
    end
end
plot(x , y/sum(y) , 'LineWidth' , 2.5 , 'color' , [.25 .25 .25]);      
xlim([-20 60]);        
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2B_hist' , '-r600');

%% Fig 2D: G4, transposons, fragile sites against all the data points
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp_new.mat');
data = NaN(length(DS) , 4);
%
data(:,1) = DS.percent_underreplicated_cdc20*100;
%
thresh_G4 = 80; idx = find(DS.G4 > thresh_G4);
data(1:length(idx) , 4) = DS.percent_underreplicated_cdc20(idx)*100;
%
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features_new.mat');
idx = find(strcmp(DS.TYPE , 'transposable_element_gene'));
data(1:length(idx) , 2) = DS.percent_underreplicated_cdc20(idx)*100;
%
idx = find(strcmp(DS.TYPE , 'interstitial_del_dup') | strcmp(DS.TYPE , 'terminal_del_dup'));
data(1:length(idx) , 3) = DS.percent_underreplicated_cdc20(idx)*100;

clrs1 = summer(12); clrs2 = hot(12); clrs3 = parula(12); clrs4 = lines(6); clrs5 = spring(12);
clrs_set = [.65 .65 .65 ; clrs1(3,:) ; clrs2(3,:) ; clrs3(5,:) ; ...
    clrs4(4,:) ; clrs3(10,:) ; clrs3(3,:) ; clrs5(3,:)];
legend_titles = {'whole genome' ,'transposable elements' , ...
    'fragile sites' , 'G4-rich regions'};
figure('units','centimeters','position',[5 5 8 10]); hold on; grid on; set(gca , 'FontSize' , 10);
h1 = boxplot(data , 'color' , [.2 .2 .2] , 'symbol' , '');
set(h1 , 'LineWidth' , 1.5);
h = findobj(gca,'Tag','Box'); 
N = length(h);
for j=1:N
    patch(get(h(N-j+1),'XData'),get(h(N-j+1),'YData'), clrs_set(j,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
legend('location' , 'SouthOutside');
ylim([-15 80]);
set(gca , 'Xtick' , '');
%ylabel('Under-replication, %');
%print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/20170531/2C' , '-r300');

%% Fig 2E: underrep VS dist to the end for gene ontologies 
thresh_nderrep = 21.0745;
load('~/Develop/Mendoza__ReplicationEvolution/Data/YeastGOwUnderrep.mat');
figure('units','centimeters','position',[5 5 20 13]); hold on; grid on; set(gca , 'FontSize' , 16)
idx = find(T.SampleSize > 5);
scatter(T.median_dist_to_the_end(idx) , T.median_underreplication(idx) , 30 , [.5 .5 .5] , 'filled' , ...
    'MarkerFaceAlpha' , .5);
idx = find(T.SampleSize > 5 & T.median_underreplication > 21.0745);
clrs1 = parula(12); clrs2 = summer(12); clrs3 = lines(6); clrs4 = spring(12); clrs5 = pink(12); clrs6 = hot(12);
clrs_set = [clrs1(3,:) ; clrs3(4,:) ; clrs1(10,:) ; clrs2(4,:) ; clrs4(5,:) ; ...
    clrs1(7,:) ; clrs5(5,:) ; clrs6(3,:) ; clrs2(1,:)];
markers = {'d' , 's'};
for I = 1:length(idx)
%     errorbar(T.mean_dist_to_the_end(idx(I)) , T.mean_underreplication(idx(I)) , ...
%         T.cv_underreplication(idx(I)) , 'LineWidth' , 1 , 'color' , [.3 .3 .3]);
%     herrorbar(T.mean_dist_to_the_end(idx(I)) , T.mean_underreplication(idx(I)) , ...
%         T.cv_dist_to_the_end(idx(I)) );
    %T.median_underreplication(idx(I)) 
    if T.median_dist_to_the_end(idx(I)) < 50
        scatter(T.median_dist_to_the_end(idx(I)) , T.median_underreplication(idx(I)) ,100 , clrs3(4,:) , 'filled', 'd' , 'MarkerEdgeColor' , 'w' , 'Display' , T.FamilyName{idx(I)} );
     T.FamilyName{idx(I)}
    else
        scatter(T.median_dist_to_the_end(idx(I)) , T.median_underreplication(idx(I)) ,100 , clrs4(5,:) , 'filled', 'd' ,'Display' , T.FamilyName{idx(I)} );
        text(T.median_dist_to_the_end(idx(I))+5 , T.median_underreplication(idx(I))+2 , T.FamilyName{idx(I)} , ...
            'FontSize' , 14);
    end
 end
%legend('location' , 'EastOutside');
ylim([-5 50]); xlim([0 500]);
%xlabel('Distance to telomeric start, kb');
%ylabel('Under-replication, %');
%print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/20170531/2D_median' , '-r300');

% legend 2E
load('~/Develop/Mendoza__ReplicationEvolution/Data/YeastGOwUnderrep.mat');
figure('units','centimeters','position',[5 5 20 13]); hold on; grid on; set(gca , 'FontSize' , 12)
idx = find(T.SampleSize > 7 & T.mean_underreplication > 21.0745);
clrs1 = parula(12); clrs2 = summer(12); clrs3 = lines(6); clrs4 = spring(12); clrs5 = pink(12); clrs6 = hot(12);
clrs_set = [clrs1(3,:) ; clrs3(4,:) ; clrs1(10,:) ; clrs2(4,:) ; clrs4(5,:) ; ...
    clrs1(7,:) ; clrs5(5,:) ; clrs6(3,:) ; clrs2(1,:)];
for I = 1:length(idx)
    scatter(T.mean_dist_to_the_end(idx(I)) , T.mean_underreplication(idx(I)) , 100 , clrs_set(I,:) , 'filled','d' ,...
        'Display' , strcat( T.FamilyName{idx(I)} , sprintf(', med. %% under-replication = %.2f' , T.median_underreplication(idx(I)))) );
end
legend('location' , 'best');
ylim([-5 50]); xlim([0 500]);
xlabel('Distance to telomeric start, kb');
ylabel('Under-replication, %');
%print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/20170531/2D_LEGEND' , '-r300');

%% Fig 2F: boxplots for SNPs and InDel frequency (for 200bp-windows) for different under-replication bins & gene loss

load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp_new.mat');
DS.underrep_bin = NaN(length(DS) , 1);
for I = 1:length(DS)
    if DS.percent_underreplicated_cdc20(I)*100 <= 20
        DS.underrep_bin(I) = 1;
    elseif DS.percent_underreplicated_cdc20(I)*100 <= 30
        DS.underrep_bin(I) = 2;
    elseif DS.percent_underreplicated_cdc20(I)*100 <= 40
        DS.underrep_bin(I) = 3;
    else
        DS.underrep_bin(I) = 4;
    end
end
clrs1 = winter(8);
figure('units','centimeters','position',[5 5 4 9]);
legend_titles = {'<5' , '5-10' , '10-20' , '>20'};
hold on; grid on; set(gca , 'FontSize' , 12);
h1 = boxplot(DS.freq_SNP , DS.underrep_bin , 'color' , [.2 .2 .2] , 'symbol' , '');
ylim([0 .175]);
set(h1 , 'LineWidth' , 1.5);
h = findobj(gca,'Tag','Box'); set(h1 , 'LineWidth' , 1.5);
for j=1:4
    patch(get(h(4-j+1),'XData'),get(h(4-j+1),'YData'), clrs1(2*j,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
set(gca , 'Xtick' , []); 
%xlabel('% unreplicated cells');
%ylabel('SNP frequency');
%legend('location' , 'SouthOutside');
%print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2I_SNP' , '-r300');

figure('units','centimeters','position',[5 5 4 9]); hold on; grid on; set(gca , 'FontSize' , 12);
h1 = boxplot(DS.freq_indel , DS.underrep_bin , 'color' , [.2 .2 .2] , 'symbol' , '');
ylim([0 .175]);
set(h1 , 'LineWidth' , 1.5);
h = findobj(gca,'Tag','Box'); set(h1 , 'LineWidth' , 1.5);
for j=1:4
    patch(get(h(4-j+1),'XData'),get(h(4-j+1),'YData'), clrs1(2*j,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
set(gca , 'Xtick' , []); 
%xlabel('% unreplicated cells');
%ylabel('InDel frequency');
%legend('location' , 'SouthOutside');
%print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2I_Indel' , '-r300');

% Fig 2J: boxplot for gene preservation (for ORFs) for different under-replication bins

load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features_new.mat');
idx = find(~isnan(DS.frac_conservation)); DS = DS(idx , :);
DS.underrep_bin = NaN(length(DS) , 1);
for I = 1:length(DS)
    if DS.percent_underreplicated_cdc20(I)*100 <= 20
        DS.underrep_bin(I) = 1;
    elseif DS.percent_underreplicated_cdc20(I)*100 <= 30
        DS.underrep_bin(I) = 2;
    elseif DS.percent_underreplicated_cdc20(I)*100 <= 40
        DS.underrep_bin(I) = 3;
    else
        DS.underrep_bin(I) = 4;
    end
end
clrs1 = winter(8);
figure('units','centimeters','position',[5 5 4 9]);
%legend_titles = {'<20' , '20-30' , '30-40' , '>40'};
%legend_titles = {'<5' , '5-10' , '10-20' , '>20'};
hold on; grid on; set(gca , 'FontSize' , 12);
h1 = boxplot(DS.frac_conservation*100 , DS.underrep_bin , 'color' , [.2 .2 .2] , 'symbol' , '');
ylim([85 100]);
set(h1 , 'LineWidth' , 1.5);
h = findobj(gca,'Tag','Box'); set(h1 , 'LineWidth' , 1.5);
for j=1:4
    patch(get(h(4-j+1),'XData'),get(h(4-j+1),'YData'), clrs1(2*j,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
set(gca , 'Xtick' , []); 
%xlabel('% unreplicated cells');
%ylabel('% strains with preserved gene');
%legend('location' , 'SouthOutside');
%print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2I_preservation' , '-r300');


%% Fig 2G: under-replication against relative fitness from Baryshnikova
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features_new.mat');
idx = find(strcmp(DS.TYPE , 'ORF')); D = DS(idx , :);
D = D(: , {'ORF' , 'dist_to_the_end_kb' , 'percent_underreplicated_cdc20' , ...
    'percent_underreplicated_Tsveti_M' , 'percent_underreplicated_dbf2' , 'Phenotype_Observable'});
T = dataset('file' , '~/Develop/Mendoza__ReplicationEvolution/Data/ExternalData/Baryshnikova10.tab');
T = T(: , {'ORF' , 'fitness'}); T = unique(T);
D = join(D , T , 'Type' , 'left' , 'Keys' , 'ORF' , 'MergeKeys' , true);
idx = find(strcmp(D.Phenotype_Observable , 'Inviable'));
for I = 1:length(idx)
    D.fitness(idx(I)) = 0;
end
idx = find(~isnan(D.fitness)); D = D(idx , :);
clrs = winter(12);
figure('units','centimeters','position',[5 5 7 9]); hold on; grid on;
%scatter(D.percent_underreplicated_cdc20*100 , D.fitness , 30 , [.5 .5 .5] , 'filled' , 'MarkerFaceAlpha' , .1);
dscatter(D.percent_underreplicated_cdc20*100 , D.fitness);
colorbar;
colormap('parula');
xlim([-20 60]);
xlabel('Under-replication , %');
ylabel('Relative to WT fitness');
D = sortrows(D , 'percent_underreplicated_cdc20');  
R = corrcoef(D.percent_underreplicated_cdc20*100 , D.fitness); 
title(sprintf('R = %.1f.' , R(1,2)));
%print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2H_dscatter' , '-r300');

%%
% ---Supplementary figures for fig 2 below


%% S4: all chromosomal arms, figure type Fig 2A
load('~/Develop/Mendoza__ReplicationEvolution/Data/MutantChr200_Michi.mat');
% left arms
DS = G;
chr_num = 5; K = 75; 
unq_mutant = {'cdc20' , 'dbf2'};
%unq_mutant = {'M'};
clrs1 = summer(12); clrs2 = hot(12); clrs_set = [clrs1(3,:) ; clrs2(3,:)];
figure('units','centimeters','position',[5 5 15 20]);
for L = 1:16
    subplot(4,4,L); hold on; grid on; set(gca , 'FontSize' , 10);
    for I = 1:length(unq_mutant)
        idx = find( strcmp(DS.MutantID , unq_mutant{I}) & DS.chr_num == L);
        data = NaN( length(idx) , K);
        for J = 1:length(idx)
            y = DS.percent_underreplicated{idx(J)}*100;
            start_point_kb = DS.start_point_kb{idx(J)};
            start_point_kb = round(start_point_kb);
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
    xlim([0 K]); ylim([-9 60]);
    set(gca , 'Xtick' , [0 25 50 75] , 'Ytick' , [0 20 40 60]);
    title(strcat ( DS.chr{idx(1)} ));
    %xlabel('Distance to the end, kbp'); ylabel('% unreplicated cells');
end
%print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/20170531/Sup_2A_left' , '-r300');

% right arms
load('~/Develop/Mendoza__ReplicationEvolution/Data/MutantChr200_Michi.mat');
DS = G;
chr_num = 5; K = 75; 
unq_mutant = {'cdc20' , 'dbf2'};
%unq_mutant = {'M'};
clrs1 = summer(12); clrs2 = hot(12); clrs_set = [clrs1(3,:) ; clrs2(3,:)];
figure('units','centimeters','position',[5 5 15 20]);
for L = 1:16
    subplot(4,4,L); hold on; grid on; set(gca , 'FontSize' , 10);
    for I = 1:length(unq_mutant)
        idx = find( strcmp(DS.MutantID , unq_mutant{I}) & DS.chr_num == L);
        data = NaN( length(idx) , K);
        for J = 1:length(idx)
            y = DS.percent_underreplicated{idx(J)}*100;
            start_point_kb = DS.start_point_kb{idx(J)};
            start_point_kb = round(start_point_kb);
            for Z = 1:K
                idx_current_kb = find(start_point_kb == max(start_point_kb)-Z+1);
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
    xlim([0 K]); ylim([-9 60]);
    set(gca , 'Xtick' , [0 25 50 75] , 'Ytick' , [0 20 40 60]);
    title(strcat ( DS.chr{idx(1)} ));
    %xlabel('Distance to the end, kbp'); ylabel('% unreplicated cells');
end
%print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/20170531/Sup_2A_right' , '-r300');

%% S5: histogrmas,across all chromosmal arms, length of underreplication & percent of unreplicated cells

load('~/Develop/Mendoza__ReplicationEvolution/Data/MutantChr200_Michi.mat');
DS = G;
figure('units','centimeters','position',[5 5 11 11]);
clrs1 = summer(12); unq_mutant = {'cdc20'}; chr_nums = [1:16];

% Length of under-replication
subplot(2,1,1); hold on; grid on; set(gca , 'FontSize' , 8);
data = [];
for I1 = 1:length(unq_mutant)
    for I2 = 1:length(chr_nums)
        idx = find(strcmp(DS.MutantID , unq_mutant{I1}) & DS.chr_num == chr_nums(I2));
        data = [data nanmedian(DS.length_underreplicated_left(idx)) nanmedian(DS.length_underreplicated_right(idx)) ];
    end
end
h = histogram(data , [30:2:70] );
set(h , 'FaceColor' , clrs1(4,:) , 'EdgeColor' , [.2 .2 .2]);
xlabel('Length of under-replicated subtelomeric region, kbp');
ylabel('# of chromosomal ends');
title(sprintf('Mean = %.2f.' , nanmean(data)));

% Percent of unreplicated cells at the last kbp
subplot(2,1,2); hold on; grid on; set(gca , 'FontSize' , 8);
data = []; 
unq_mutant = {'cdc20'}; chr_nums = [1:16];
for I1 = 1:length(unq_mutant)
    for I2 = 1:length(chr_nums)
        idx = find(strcmp(DS.MutantID , unq_mutant{I1}) & DS.chr_num == chr_nums(I2));
        data = [data nanmedian(DS.percent_underreplicated_left(idx))*100 nanmedian(DS.percent_underreplicated_right(idx))*100 ];
    end
end
h = histogram(data , [40:2:70] );
set(h , 'FaceColor' , clrs1(4,:) , 'EdgeColor' , [.2 .2 .2]);
xlabel('% unreplicated cells at the last kbp');
ylabel('# of chromosomal ends');
title(sprintf('Mean = %.2f.' , nanmean(data)));

%print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/SupAlsu2' , '-r300');


%% S6B: histograms of under-replication, binned by replication timing

load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp_new.mat');
clrs1 = summer(12); clrs2 = pink(12); clrs3 = lines(6); clrs4 = winter(12);
clrs_set = [clrs1(6,:) ; clrs4(7,:) ; clrs3(4,:) ; clrs4(7,:)];
figure('units','centimeters','position',[5 5 10 10]); hold on; grid on; set(gca , 'FontSize' , 15);
data = DS.Trep_spline; 
trep_thresh = [0 30 45 65];
legend_names = {'<30 minutes' , '30-45 minutes' , '>45 minutes'};

for I = 1:length(trep_thresh)-1
    idx = find(DS.Trep_spline >= trep_thresh(I) & DS.Trep_spline < trep_thresh(I+1) );
    if ~isempty(idx)
        [y,x] = ksdensity(DS.percent_underreplicated_cdc20(idx)*100 , [-20:1:70]);
        plot(x,y/sum(y) , 'LineWidth' , 2.5 , 'color' , clrs_set(I,:));
    end
end
xlim([-30 70]);
legend(legend_names , 'location' , 'ne');
% ylabel('Probability density'); xlabel('Under-replication, %');
%print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2B' , '-r300');

%% S6A: violin plot for under-replication, binned by replication timing
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp_new.mat');
data = NaN(length(DS) , 3); 
trep_thresh = [0 30 45 65];
legend_names = {'<30' , '30-45' , '>45'};

for I = 1:length(trep_thresh)-1
    idx = find(DS.Trep_spline >= trep_thresh(I) & DS.Trep_spline < trep_thresh(I+1) & ...
        DS.percent_underreplicated_cdc20 > -5);
    data(1:length(idx) , I) = DS.percent_underreplicated_cdc20(idx)*100;
end

clrs1 = summer(12); clrs2 = pink(12); clrs3 = lines(6); clrs4 = winter(12);
clrs_set = [clrs1(6,:) ; clrs4(7,:) ; clrs3(4,:) ; clrs4(7,:)];
figure('units','centimeters','position', [5 5 10 10]); hold on; grid on; set(gca , 'FontSize' , 15);
violin(data , 'facecolor',clrs_set,'edgecolor','none','bw',5,'mc',[] , 'medc' , []);
ylim([-20 70]);
set(gca , 'Xtick' , [1:3] , 'XtickLabel' , legend_names)
%print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/20170531/repTime__underrep__violin_Opt3' , '-r300');


%% S7 : like 2D, but for dbf2 (whole genome, transposons, fragile sites, G4)
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp_new.mat');
data = NaN(length(DS) , 4);
%
data(:,1) = DS.percent_underreplicated_dbf2*100;
%
thresh_G4 = 80; idx = find(DS.G4 > thresh_G4);
data(1:length(idx) , 4) = DS.percent_underreplicated_dbf2(idx)*100;
%
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features_new.mat');
idx = find(strcmp(DS.TYPE , 'transposable_element_gene'));
data(1:length(idx) , 2) = DS.percent_underreplicated_dbf2(idx)*100;
%
idx = find(strcmp(DS.TYPE , 'interstitial_del_dup') | strcmp(DS.TYPE , 'terminal_del_dup'));
data(1:length(idx) , 3) = DS.percent_underreplicated_dbf2(idx)*100;

clrs1 = summer(12); clrs2 = hot(12); clrs3 = parula(12); clrs4 = lines(6); clrs5 = spring(12);
clrs_set = [.65 .65 .65 ; clrs1(3,:) ; clrs2(3,:) ; clrs3(5,:) ; ...
    clrs4(4,:) ; clrs3(10,:) ; clrs3(3,:) ; clrs5(3,:)];
legend_titles = {'whole genome' ,'transposable elements' , ...
    'fragile sites' , 'G4-rich regions'};
figure('units','centimeters','position',[5 5 8 10]); hold on; grid on; set(gca , 'FontSize' , 10);
h1 = boxplot(data , 'color' , [.2 .2 .2] , 'symbol' , '');
set(h1 , 'LineWidth' , 1.5);
h = findobj(gca,'Tag','Box'); 
N = length(h);
for j=1:N
    patch(get(h(N-j+1),'XData'),get(h(N-j+1),'YData'), clrs_set(j,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
legend('location' , 'SouthOutside');
ylim([-25 70]);
set(gca , 'Xtick' , '');
%ylabel('Under-replication, %');
%print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/20170531/2C' , '-r300');


