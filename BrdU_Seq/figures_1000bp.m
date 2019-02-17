%%%%% script to generate figures regardibg BrUSeq

%% Fig1. Check gDNA coverage distribution for all the samples
load('~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000.mat');
idx = find(~strcmp(DS.chr , 'chrMito'));
data = NaN(length(idx) , 6+5+3);
data(:,1) = DS.cov_MG1_1(idx); data(:,2) = DS.cov_MG1_2(idx); data(:,3) = DS.cov_MG1_3(idx); data(:,4) = DS.cov_MG1_4(idx); data(:,5) = DS.cov_MG1_5(idx); data(:,6) = DS.cov_MG1_6(idx); 
data(:,7) = DS.cov_MM_1(idx); data(:,8) = DS.cov_MM_2(idx); data(:,9) = DS.cov_MM_3(idx); data(:,10) = DS.cov_MM_4(idx); data(:,11) = DS.cov_MM_5(idx);
data(:,12) = DS.cov_MG1noBrdU_1(idx); data(:,13) = DS.cov_MG1noBrdU_2(idx); data(:,14) = DS.cov_MG1noBrdU_3(idx); 

figure('units','centimeters','position',[5 5 45 20]); hold on; grid on; set(gca , 'FontSize' , 12);
h = boxplot(data , 'color' , [.2 .2 .2] , 'symbol' , '');
ylim([0 600]); ylabel('Coverage'); title('gDNA');
set(gca , 'Xtick' , [1:14] , 'XtickLabel' , {'G1-1','G1-2','G1-3','G1-4','G1-5','G1-6',...
    'M-1','M-2','M-3','M-4','M-5' , 'noBr-1' , 'noBr-2' , 'noBr-3'});
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/figures/boxplot__coverage__gDNA' , '-r300');

%% Fig2. average mtDNA coveage/average coverage -- average coverage
load('~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000.mat');
data = NaN(length(DS) , 6+5+3);
data(:,1) = DS.cov_MG1_1; data(:,2) = DS.cov_MG1_2; data(:,3) = DS.cov_MG1_3; data(:,4) = DS.cov_MG1_4; data(:,5) = DS.cov_MG1_5; data(:,6) = DS.cov_MG1_6; 
data(:,7) = DS.cov_MM_1; data(:,8) = DS.cov_MM_2; data(:,9) = DS.cov_MM_3; data(:,10) = DS.cov_MM_4; data(:,11) = DS.cov_MM_5;
data(:,12) = DS.cov_MG1noBrdU_1; data(:,13) = DS.cov_MG1noBrdU_2; data(:,14) = DS.cov_MG1noBrdU_3; 
idx_gDNA = find(~strcmp(DS.chr , 'chrMito'));
idx_mtDNA = find(strcmp(DS.chr , 'chrMito'));
labels = {'G1-1','G1-2','G1-3','G1-4','G1-5','G1-6',...
    'M-1','M-2','M-3','M-4','M-5' , 'noBr-1' , 'noBr-2' , 'noBr-3'};

figure('units','centimeters','position',[5 5 45 20]); 
clrs1 = winter(12); clrs2 = summer(12); clrs3 = hot(12);
clrs_set = [repmat(clrs1(7,:) , 6 , 1); repmat(clrs2(8,:) , 5 , 1); repmat(clrs3(3,:) , 3 , 1);];

subplot(1,2,2);hold on; grid on; set(gca , 'FontSize' , 12);
for I = 1:14
    scatter(nanmedian(data(idx_gDNA,I)) , nanmedian(data(idx_mtDNA,I))/nanmedian(data(idx_gDNA,I)) , 150 , clrs_set(I,:) , 'filled');
    text(nanmedian(data(idx_gDNA,I))+1 , nanmedian(data(idx_mtDNA,I))/nanmedian(data(idx_gDNA,I))+1 , labels{I} , 'FontSize' , 12);
end
xlabel('median cov. gDNA'); ylabel('median cov. mtDNA/ median cov. gDNA');

subplot(1,2,1);hold on; grid on; set(gca , 'FontSize' , 12);
for I = 1:14
    scatter(nanmedian(data(idx_gDNA,I)) , nanmedian(data(idx_mtDNA,I)) , 150 , clrs_set(I,:) , 'filled');
    text(nanmedian(data(idx_gDNA,I))+1 , nanmedian(data(idx_mtDNA,I))+1 , labels{I} , 'FontSize' , 12);
end
xlabel('median cov. gDNA'); ylabel('median cov. mtDNA');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/figures/scatter__coverage__mt_and_g' , '-r300');

%% Fig 3. Comapring coverage with other features

load('~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000.mat');
D = dataset('file' , '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/BrdU_Stat.tab');
data = NaN(length(DS) , 6+5+3);
data(:,1) = DS.cov_MG1_1; data(:,2) = DS.cov_MG1_2; data(:,3) = DS.cov_MG1_3; data(:,4) = DS.cov_MG1_4; data(:,5) = DS.cov_MG1_5; data(:,6) = DS.cov_MG1_6; 
data(:,7) = DS.cov_MM_1; data(:,8) = DS.cov_MM_2; data(:,9) = DS.cov_MM_3; data(:,10) = DS.cov_MM_4; data(:,11) = DS.cov_MM_5;
data(:,12) = DS.cov_MG1noBrdU_1; data(:,13) = DS.cov_MG1noBrdU_2; data(:,14) = DS.cov_MG1noBrdU_3; 
idx_gDNA = find(~strcmp(DS.chr , 'chrMito'));
idx_mtDNA = find(strcmp(DS.chr , 'chrMito'));

% 4 features: number of mapped reads; fraction of mapped reads; synchrony;
% molarity after DNA prep;

data_features = NaN(14 , 4);
data_features(:,1) = D.NumReadsMapped(1:14);
data_features(:,2) = D.NumReadsMapped(1:14)./D.NumReads(1:14);
data_features(:,3) = D.Synchrony(1:14);
data_features(:,4) = D.MolarityAfterLibraryPrep__nM(1:14);

clrs1 = winter(12); clrs2 = summer(12); clrs3 = hot(12);
clrs_set = [repmat(clrs1(7,:) , 6 , 1); repmat(clrs2(8,:) , 5 , 1); repmat(clrs3(3,:) , 3 , 1);];

labels = {'G1-1','G1-2','G1-3','G1-4','G1-5','G1-6',...
    'M-1','M-2','M-3','M-4','M-5' , 'noBr-1' , 'noBr-2' , 'noBr-3'};

feature_titles = {'Number of mapped reads' , 'Frac. of mapped reads' , 'Synchrony' , 'Molarity'};
save_names  = {'numMappedReads' , 'fracMappedReads' , 'synchrony' , 'molarity'};

for J = 1:4
    figure('units','centimeters','position',[5 5 45 20]); 
    subplot(1,2,1); hold on; grid on; set(gca , 'FontSize' , 12);
    for I = 1:14
        scatter(nanmedian(data(idx_gDNA,I)) , data_features(I,J) , 100 , clrs_set(I,:) , 'filled');
        text(nanmedian(data(idx_gDNA,I)) , data_features(I,J) , labels{I} , 'FontSize' , 12);
    end
    xlabel('cov. gDNA'); ylabel(feature_titles{J});
    subplot(1,2,2); hold on; grid on; set(gca , 'FontSize' , 12);
    for I = 1:14
        scatter(nanmedian(data(idx_mtDNA,I))/nanmedian(data(idx_gDNA,I)) , data_features(I,J) , 100 , clrs_set(I,:) , 'filled');
        text(nanmedian(data(idx_mtDNA,I))/nanmedian(data(idx_gDNA,I)) , data_features(I,J) , labels{I} , 'FontSize' , 12);
    end
    xlabel('cov. mtDNA / cov. gDNA'); ylabel(feature_titles{J});
    save_name = strcat('~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/figures/coverage__VS__' , save_names{J} );
    print('-dpng' , save_name , '-r300');
end


%% Fig 4. Coverage and CN -- VS -- distance to the end for all samples in M-G1, M-M and M-G1noBrdU

load('~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000.mat');
idx = find(~strcmp(DS.chr , 'chrMito'));
data = NaN(length(idx) , 6+5+3);
data(:,1) = DS.cov_MG1_1(idx); data(:,2) = DS.cov_MG1_2(idx); data(:,3) = DS.cov_MG1_3(idx); data(:,4) = DS.cov_MG1_4(idx); data(:,5) = DS.cov_MG1_5(idx); data(:,6) = DS.cov_MG1_6(idx); 
data(:,7) = DS.cov_MM_1(idx); data(:,8) = DS.cov_MM_2(idx); data(:,9) = DS.cov_MM_3(idx); data(:,10) = DS.cov_MM_4(idx); data(:,11) = DS.cov_MM_5(idx);
data(:,12) = DS.cov_MG1noBrdU_1(idx); data(:,13) = DS.cov_MG1noBrdU_2(idx); data(:,14) = DS.cov_MG1noBrdU_3(idx); 

title_names = {'G1-1','G1-2','G1-3','G1-4','G1-5','G1-6','M-1','M-2','M-3','M-4','M-5','noBr-1' , 'noBr-2' , 'noBr-3'};
ylims = [700 150 300 700 300 300 300 300 30 700 300 200 300 350];
figure('units','centimeters','position',[5 5 45 20]); 
clrs = winter(12);
for I = 1:14
    subplot(4,4,I);hold on; grid on; set(gca , 'FontSize' , 12);
    scatter(DS.dist_to_the_end_kb(idx) , data(:,I) , 15 , clrs(7,:) , 'filled');
    xlabel('Distance to the end'); ylabel('Coverage'); title(title_names{I});
    ylim([0 ylims(I)]);
end
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/figures/scatter__coverage__VS__dist' , '-r300');

data = NaN(length(idx) , 6+5+3);
data(:,1) = DS.CN_MG1_1(idx); data(:,2) = DS.CN_MG1_2(idx); data(:,3) = DS.CN_MG1_3(idx); data(:,4) = DS.CN_MG1_4(idx); data(:,5) = DS.CN_MG1_5(idx); data(:,6) = DS.CN_MG1_6(idx); 
data(:,7) = DS.CN_MM_1(idx); data(:,8) = DS.CN_MM_2(idx); data(:,9) = DS.CN_MM_3(idx); data(:,10) = DS.CN_MM_4(idx); data(:,11) = DS.CN_MM_5(idx);
data(:,12) = DS.CN_MG1noBrdU_1(idx); data(:,13) = DS.CN_MG1noBrdU_2(idx); data(:,14) = DS.CN_MG1noBrdU_3(idx); 

title_names = {'G1-1','G1-2','G1-3','G1-4','G1-5','G1-6','M-1','M-2','M-3','M-4','M-5','noBr-1' , 'noBr-2' , 'noBr-3'};
figure('units','centimeters','position',[5 5 45 20]); 
clrs = winter(12);
for I = 1:14
    subplot(4,4,I);hold on; grid on; set(gca , 'FontSize' , 12);
    scatter(DS.dist_to_the_end_kb(idx) , data(:,I) , 15 , clrs(7,:) , 'filled');
    xlabel('Distance to the end'); ylabel('Coverage'); title(title_names{I});
    ylim([0 3]);
end
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/figures/scatter__CN__VS__dist' , '-r300');

%% Fig 5. Functional analysis of loci that have high CN: 
% add if in all the replicates across the same stage: CN > 1.5 
load('~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000.mat');
DS.G1_highCN = NaN(length(DS) , 1);
DS.M_highCN = NaN(length(DS) , 1);
for I = 1:length(DS)
    if ~strcmp(DS.chr{I} , 'chrMito')
        DS.G1_highCN(I) = sum(DS.CN_MG1_1(I) > 1.5) + sum(DS.CN_MG1_4(I) > 1.5) ;
        DS.M_highCN(I) = sum(DS.CN_MM_4(I) > 1.5);
    end
end

% 5A: boxplots for under-replication, distance to the end, RT, distance to
% the closest ARS

feature_titles = {'Under-replication, %' , 'Dist. to the end' , 'RT' , 'Dist. to ARS'};
ylims = [-0.3 0.7 ; 0 600 ; 0 65 ; 0 20];
idx = find(~strcmp(DS.chr , 'chrMito'));
data = [DS.percent_underreplicated_cdc20(idx) DS.dist_to_the_end_kb(idx) DS.Trep_spline(idx) DS.dist_to_ARS(idx)];
figure('units','centimeters','position',[5 5 45 20]); 
for J = 1:4
    subplot(2,4,J); hold on; grid on; set(gca , 'FontSize' , 10);
    boxplot( data(:,J) , DS.G1_highCN(idx) , 'color' , [.2 .2 .2] , 'symbol' , '');
    set(gca , 'Xtick' , [1:3] , 'XtickLabel' , {'0 CN > 1.5' , '1 CN > 1.5' , '2 CN > 1.5'});
    xlabel('BrdU incorporated'); ylabel(feature_titles{J}); title('G1');
    ylim(ylims(J,:));
    
    subplot(2,4,J+4); hold on; grid on; set(gca , 'FontSize' , 10);
    boxplot(data(:,J) , DS.M_highCN(idx), 'color' , [.2 .2 .2] , 'symbol' , '');
    set(gca , 'Xtick' , [1:2] , 'XtickLabel' , {'CN < 1.5' , 'CN > 1.5'});
    xlabel('BrdU incorporated'); ylabel(feature_titles{J}); title('M');
    ylim(ylims(J,:));
end
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/figures/features__VS__highCN' , '-r300');

% 5B. Barplots for fraction of regins being a transposon for regions that
% incorporated BrdU and not
clrs = winter(4*3);
figure('units','centimeters','position',[5 5 45 20]);
subtel_thresh = [0 50];
titles_subtel = {'all gDNA' , 'not subtelomeric'};

for J = 1:2
    subplot(2,2,2*J-1); hold on; grid on; set(gca , 'FontSize' , 12);
    for I = 1:3
        bar(I , sum(DS.transposable_bool == 1 & DS.G1_highCN == I-1 & DS.dist_to_the_end_kb > subtel_thresh(J))/sum(DS.G1_highCN == I-1 & DS.dist_to_the_end_kb > subtel_thresh(J)) , 'FaceColor' , clrs(4*I-2,:));
    end
    set(gca , 'Xtick' , [1:3] , 'XtickLabel' , {'0 CN > 1.5' , '1 CN > 1.5' , '2 CN > 1.5'});
    ylabel('Fraction of transposable elements'); title(strcat('G1,' , titles_subtel{J}));

    subplot(2,2,2*J); hold on; grid on; set(gca , 'FontSize' , 12);
    for I = 1:2
        bar(I , sum(DS.transposable_bool == 1 & DS.M_highCN == I-1 & DS.dist_to_the_end_kb > subtel_thresh(J))/sum(DS.M_highCN == I-1 & DS.dist_to_the_end_kb > subtel_thresh(J)), 'FaceColor' , clrs(4*I-2,:));
    end
    set(gca , 'Xtick' , [1:2] , 'XtickLabel' , {'CN < 1.5' , 'CN > 1.5'});
    ylabel('Fraction of transposable elements'); title(strcat('M,' , titles_subtel{J}));
end

print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/figures/transposableEffect' , '-r300');

%% Fig 6. For each region: comparing CN between G1-1; G1-4 and M-4

load('~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000.mat');
clrs = winter(12);
figure; 

subplot(2,2,1); hold on; grid on; set(gca , 'FontSize' , 12);
idx = find(~strcmp(DS.chr , 'chrMito') & DS.rDNA_bool == 0);
scatter(DS.CN_MG1_1(idx) , DS.CN_MG1_4(idx) , 15 , clrs(6,:) , 'filled');
R = corrcoef(DS.CN_MG1_1(idx) , DS.CN_MG1_4(idx));
title(sprintf('R = %.2f.',R(1,2)));
plot([0 3] , [0 3] , 'LineWidth' , 1.5 , 'color' , [.2 .2 .2]);
xlabel('G1-1'); ylabel('G1-4');

subplot(2,2,2); hold on; grid on; set(gca , 'FontSize' , 12);
idx = find(~strcmp(DS.chr , 'chrMito') & DS.rDNA_bool == 0);
scatter(DS.CN_MG1_4(idx) , DS.CN_MM_4(idx) , 15 , clrs(6,:) , 'filled');
plot([0 3] , [0 3] , 'LineWidth' , 1.5 , 'color' , [.2 .2 .2]);
R = corrcoef(DS.CN_MM_4(idx) , DS.CN_MG1_4(idx));
title(sprintf('R = %.2f.',R(1,2)));
xlabel('G1-4'); ylabel('M-4');

subplot(2,2,3); hold on; grid on; set(gca , 'FontSize' , 12);
idx = find(~strcmp(DS.chr , 'chrMito') & DS.rDNA_bool == 0);
scatter(DS.CN_MM_4(idx) , DS.CN_MG1_1(idx) , 15 , clrs(6,:) , 'filled');
plot([0 3] , [0 3] , 'LineWidth' , 1.5 , 'color' , [.2 .2 .2]);
xlabel('M-4'); ylabel('G1-1');
plot([0 3] , [0 3] , 'LineWidth' , 1.5 , 'color' , [.2 .2 .2]);
R = corrcoef(DS.CN_MM_4(idx) , DS.CN_MG1_1(idx));
title(sprintf('R = %.2f.',R(1,2)));

print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/figures/comparison__CN__M_and_G1' , '-r300');
%% Fig 7. rDNA
load('~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000.mat');
DS.vec_CNR_MG1_MM = cell(length(DS) , 1);
DS.median_CNR_MG1_MM = NaN(length(DS) , 1);
DS.mean_CNR_MG1_MM = NaN(length(DS) , 1);
for I = 1:length(DS)
    DS.vec_CNR_MG1_MM{I} = [DS.CN_MG1_3(I)/DS.CN_MM_1(I) DS.CN_MG1_4(I)/DS.CN_MM_2(I) ];
    DS.median_CNR_MG1_MM(I) = nanmedian(DS.vec_CNR_MG1_MM{I});
    DS.mean_CNR_MG1_MM(I) = nanmean(DS.vec_CNR_MG1_MM{I});
end
figure; hold on; grid on; set(gca , 'FontSize' , 12);
idx = find(~strcmp(DS.chr , 'chrMito'));
boxplot(DS.median_CNR_MG1_MM(idx) , DS.rDNA_bool(idx) , 'color' , [.2 .2 .2] , 'symbol' , '');
xlabel('rDNA');
ylabel('CN-ratio (G1/M)');
ylim([.5 2]);
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/figures/rDNA' , '-r300');


%% Analyses 2. Add CN-ratio across replicates
% 4 options -- 
% a) only account samples with median coverage > 300 (G1-1, G1-4, M-4 );
% b) only account samples with median coverage > 100 (G1-1, G1-4, G1-5, G1-6, M-1, M-4, M-5, noBr-2 , noBr-3)
% c) account everything, except M3 which has average coverage ~5
% d) account only matched up by date(i.e. coming from the same sample)
% replicates; for G1-M: G1-3 -- M-1; G1-4 -- M-2; exclude G1-5 -- M-3 since M-3 has
% too low coverage; for G1-G1noBrdU: G1-4 -- noBr1; G1-6 -- noBr2
% e) discard samples with low molarity  or samples for which Nuria observed
% weird profile

% A: G1-1, G1-4; M-4
load('~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000.mat');
% add ratio for CN
DS.median_CN_MM = NaN(length(DS) , 1);
DS.mean_CN_MM = NaN(length(DS) , 1);
DS.median_CN_MG1 = NaN(length(DS) , 1);
DS.mean_CN_MG1 = NaN(length(DS) , 1);
DS.median_CNR_MG1_MM = NaN(length(DS) , 1);
DS.mean_CNR_MG1_MM = NaN(length(DS) , 1);
for I = 1:length(DS)
    DS.median_CN_MM(I) = DS.CN_MM_4(I);
    DS.median_CN_MG1(I) = nanmedian([DS.CN_MG1_1(I) DS.CN_MG1_4(I)]);
    DS.median_CNR_MG1_MM(I) = DS.median_CN_MG1(I)/DS.median_CN_MM(I);
    
    DS.mean_CN_MM(I) = DS.CN_MM_4(I);
    DS.mean_CN_MG1(I) = nanmean([DS.CN_MG1_1(I) DS.CN_MG1_4(I)]);
    DS.mean_CNR_MG1_MM(I) = DS.mean_CN_MG1(I)/DS.mean_CN_MM(I);
end
save('~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000_A.mat' , 'DS');

% B: G1-1, G1-4, G1-5, G1-6; M-1, M-4, M-5
load('~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000.mat');
% add ratio for CN
DS.median_CN_MM = NaN(length(DS) , 1);
DS.mean_CN_MM = NaN(length(DS) , 1);
DS.median_CN_MG1 = NaN(length(DS) , 1);
DS.mean_CN_MG1 = NaN(length(DS) , 1);
DS.median_CNR_MG1_MM = NaN(length(DS) , 1);
DS.mean_CNR_MG1_MM = NaN(length(DS) , 1);

DS.median_CN_MG1noBrdU = NaN(length(DS) , 1);
DS.mean_CN_MG1noBrdU = NaN(length(DS) , 1);
DS.median_CNR_MG1_MG1noBrdU = NaN(length(DS) , 1);
DS.mean_CNR_MG1_MG1noBrdU = NaN(length(DS) , 1);
for I = 1:length(DS)
    DS.median_CN_MM(I) = nanmedian([DS.CN_MM_1(I) DS.CN_MM_4(I) DS.CN_MM_5(I)]);
    DS.median_CN_MG1(I) = nanmedian([DS.CN_MG1_1(I) DS.CN_MG1_4(I) DS.CN_MG1_5(I) DS.CN_MG1_6(I)]);
    DS.median_CNR_MG1_MM(I) = DS.median_CN_MG1(I)/DS.median_CN_MM(I);
    DS.mean_CN_MM(I) = nanmean([DS.CN_MM_1(I) DS.CN_MM_4(I) DS.CN_MM_5(I)]);
    DS.mean_CN_MG1(I) = nanmean([DS.CN_MG1_1(I) DS.CN_MG1_4(I) DS.CN_MG1_5(I) DS.CN_MG1_6(I)]);
    DS.mean_CNR_MG1_MM(I) = DS.mean_CN_MG1(I)/DS.mean_CN_MM(I);
    
    DS.median_CN_MG1noBrdU(I) = nanmedian([DS.CN_MG1noBrdU_2(I) DS.CN_MG1noBrdU_3(I)]);
    DS.mean_CN_MG1noBrdU(I) = nanmean([DS.CN_MG1noBrdU_2(I) DS.CN_MG1noBrdU_3(I)]);
    DS.median_CNR_MG1_MG1noBrdU(I) = DS.median_CN_MG1(I)/DS.median_CN_MG1noBrdU(I);
    DS.mean_CNR_MG1_MG1noBrdU(I) = DS.mean_CN_MG1(I)/DS.mean_CN_MG1noBrdU(I);
end
save('~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000_B.mat' , 'DS');

% C: everything except M3
load('~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000.mat');
% add ratio for CN
DS.median_CN_MM = NaN(length(DS) , 1);
DS.mean_CN_MM = NaN(length(DS) , 1);
DS.median_CN_MG1 = NaN(length(DS) , 1);
DS.mean_CN_MG1 = NaN(length(DS) , 1);
DS.median_CNR_MG1_MM = NaN(length(DS) , 1);
DS.mean_CNR_MG1_MM = NaN(length(DS) , 1);

DS.median_CN_MG1noBrdU = NaN(length(DS) , 1);
DS.mean_CN_MG1noBrdU = NaN(length(DS) , 1);
DS.median_CNR_MG1_MG1noBrdU = NaN(length(DS) , 1);
DS.mean_CNR_MG1_MG1noBrdU = NaN(length(DS) , 1);

for I = 1:length(DS)
    DS.median_CN_MM(I) = nanmedian([DS.CN_MM_1(I) DS.CN_MM_2(I) DS.CN_MM_4(I) DS.CN_MM_5(I)]);
    DS.median_CN_MG1(I) = nanmedian([DS.CN_MG1_1(I) DS.CN_MG1_2(I) DS.CN_MG1_3(I) DS.CN_MG1_4(I) DS.CN_MG1_5(I) DS.CN_MG1_6(I)]);
    DS.median_CNR_MG1_MM(I) = DS.median_CN_MG1(I)/DS.median_CN_MM(I);
    DS.mean_CN_MM(I) = nanmean([DS.CN_MM_1(I) DS.CN_MM_2(I) DS.CN_MM_4(I) DS.CN_MM_5(I)]);
    DS.mean_CN_MG1(I) = nanmean([DS.CN_MG1_1(I) DS.CN_MG1_2(I) DS.CN_MG1_3(I) DS.CN_MG1_4(I) DS.CN_MG1_5(I) DS.CN_MG1_6(I)]);
    DS.mean_CNR_MG1_MM(I) = DS.mean_CN_MG1(I)/DS.mean_CN_MM(I);
    
    DS.median_CN_MG1noBrdU(I) = nanmedian([DS.CN_MG1noBrdU_1(I) DS.CN_MG1noBrdU_2(I) DS.CN_MG1noBrdU_3(I)]);
    DS.mean_CN_MG1noBrdU(I) = nanmean([DS.CN_MG1noBrdU_1(I) DS.CN_MG1noBrdU_2(I) DS.CN_MG1noBrdU_3(I)]);
    DS.median_CNR_MG1_MG1noBrdU(I) = DS.median_CN_MG1(I)/DS.median_CN_MG1noBrdU(I);
    DS.mean_CNR_MG1_MG1noBrdU(I) = DS.mean_CN_MG1(I)/DS.mean_CN_MG1noBrdU(I);
end
save('~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000_C.mat' , 'DS');

% D: matched up by date(i.e. coming from the same sample) replicates: G1-M: G1-3
% -- M-1; G1-4 -- M-2; G1-G1noBrdU: G1-4 -- noBr1; G1-6 -- noBr2
load('~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000.mat');
% add ratio for CN (vector and median/mean)
DS.vec_CNR_MG1_MM = cell(length(DS) , 1);
DS.median_CNR_MG1_MM = NaN(length(DS) , 1);
DS.mean_CNR_MG1_MM = NaN(length(DS) , 1);
DS.vec_CNR_MG1_MG1noBrdU = cell(length(DS) , 1);
DS.median_CNR_MG1_MG1noBrdU = NaN(length(DS) , 1);
DS.mean_CNR_MG1_MG1noBrdU = NaN(length(DS) , 1);
for I = 1:length(DS)
    DS.vec_CNR_MG1_MM{I} = [DS.CN_MG1_3(I)/DS.CN_MM_1(I) DS.CN_MG1_4(I)/DS.CN_MM_2(I) ];
    DS.median_CNR_MG1_MM(I) = nanmedian(DS.vec_CNR_MG1_MM{I});
    DS.mean_CNR_MG1_MM(I) = nanmean(DS.vec_CNR_MG1_MM{I});
    
    DS.vec_CNR_MG1_MG1noBrdU{I} = [DS.CN_MG1_4(I)/DS.CN_MG1noBrdU_1(I) DS.CN_MG1_6(I)/DS.CN_MG1noBrdU_2(I) ];
    DS.median_CNR_MG1_MG1noBrdU(I) = nanmedian(DS.vec_CNR_MG1_MG1noBrdU{I});
    DS.mean_CNR_MG1_MG1noBrdU(I) = nanmean(DS.vec_CNR_MG1_MG1noBrdU{I});
end
save('~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000_D.mat' , 'DS');

% E: throw away samples with low concentration: G1-5, G1-3 and G1-4
load('~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000.mat');
% add ratio for CN
DS.median_CN_MM = NaN(length(DS) , 1);
DS.mean_CN_MM = NaN(length(DS) , 1);
DS.median_CN_MG1 = NaN(length(DS) , 1);
DS.mean_CN_MG1 = NaN(length(DS) , 1);
DS.median_CNR_MG1_MM = NaN(length(DS) , 1);
DS.mean_CNR_MG1_MM = NaN(length(DS) , 1);

DS.median_CN_MG1noBrdU = NaN(length(DS) , 1);
DS.mean_CN_MG1noBrdU = NaN(length(DS) , 1);
DS.median_CNR_MG1_MG1noBrdU = NaN(length(DS) , 1);
DS.mean_CNR_MG1_MG1noBrdU = NaN(length(DS) , 1);

for I = 1:length(DS)
    DS.median_CN_MM(I) = nanmedian([DS.CN_MM_1(I) DS.CN_MM_2(I) DS.CN_MM_5(I)]);
    DS.median_CN_MG1(I) = nanmedian([DS.CN_MG1_1(I) DS.CN_MG1_2(I) DS.CN_MG1_3(I) DS.CN_MG1_4(I) DS.CN_MG1_6(I)]);
    DS.median_CNR_MG1_MM(I) = DS.median_CN_MG1(I)/DS.median_CN_MM(I);
    DS.mean_CN_MM(I) = nanmean([DS.CN_MM_1(I) DS.CN_MM_2(I) DS.CN_MM_5(I)]);
    DS.mean_CN_MG1(I) = nanmean([DS.CN_MG1_1(I) DS.CN_MG1_2(I) DS.CN_MG1_3(I) DS.CN_MG1_4(I) DS.CN_MG1_6(I)]);
    DS.mean_CNR_MG1_MM(I) = DS.mean_CN_MG1(I)/DS.mean_CN_MM(I);
    
    DS.median_CN_MG1noBrdU(I) = nanmedian([DS.CN_MG1noBrdU_1(I) DS.CN_MG1noBrdU_2(I) DS.CN_MG1noBrdU_3(I)]);
    DS.mean_CN_MG1noBrdU(I) = nanmean([DS.CN_MG1noBrdU_1(I) DS.CN_MG1noBrdU_2(I) DS.CN_MG1noBrdU_3(I)]);
    DS.median_CNR_MG1_MG1noBrdU(I) = DS.median_CN_MG1(I)/DS.median_CN_MG1noBrdU(I);
    DS.mean_CNR_MG1_MG1noBrdU(I) = DS.mean_CN_MG1(I)/DS.mean_CN_MG1noBrdU(I);
end
save('~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000_E.mat' , 'DS');

%% Fig8. scatter of CNR against distance to the end for datasets
clrs = winter(12);
DS_names = {'~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000_A.mat' , ...
    '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000_B.mat' , ...
    '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000_C.mat' , ...
    '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000_D.mat',...
    '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000_E.mat'};
title_names = {'A' , 'B' , 'C' , 'D' , 'E'};
% Fig3A: G1-M
figure('units','centimeters','position',[5 5 45 20]); 
for I = 1:length(DS_names)
    DS_name = DS_names{I}; title_name = title_names{I};
    load(DS_name);

    subplot(length(DS_names),2,2*I-1);hold on; grid on;
    scatter(DS.dist_to_the_end_kb , DS.mean_CNR_MG1_MM , 15 , clrs(7,:) , 'filled');
    xlabel('Dist. to the end'); ylabel('mean CNR'); title(title_name);

    subplot(length(DS_names),2,2*I);hold on; grid on;
    scatter(DS.percent_underreplicated_cdc20 , DS.mean_CNR_MG1_MM , 15 , clrs(7,:) , 'filled');
	idx = find(~isnan(DS.mean_CNR_MG1_MM) & DS.mean_CNR_MG1_MM < 1000000);
    R = corrcoef(DS.percent_underreplicated_cdc20(idx) , DS.mean_CNR_MG1_MM(idx)); R = R(1,2);
    title(strcat( title_name , sprintf(': R = %.2f.' , R)));
    xlabel('Under-replication, %'); ylabel('mean CNR'); xlim([-0.3 0.6]);
end
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/figures/scatter__CNR__VS__dist__G1_and_M' , '-r300');

% Fig3B: G1-G1noBrdU
figure('units','centimeters','position',[5 5 45 20]); 
for I = 2:length(DS_names)
    DS_name = DS_names{I}; title_name = title_names{I};
    load(DS_name);

    subplot(length(DS_names),2,2*I-1);hold on; grid on;
    scatter(DS.dist_to_the_end_kb , DS.mean_CNR_MG1_MG1noBrdU , 15 , clrs(7,:) , 'filled');
    xlabel('Dist. to the end'); ylabel('mean CNR'); title(title_name);

    subplot(length(DS_names),2,2*I);hold on; grid on;
    scatter(DS.percent_underreplicated_cdc20 , DS.mean_CNR_MG1_MG1noBrdU , 15 , clrs(7,:) , 'filled');
    idx = find(~isnan(DS.mean_CNR_MG1_MG1noBrdU) & DS.mean_CNR_MG1_MG1noBrdU < 1000000);
    R = corrcoef(DS.percent_underreplicated_cdc20(idx) , DS.mean_CNR_MG1_MG1noBrdU(idx)); R = R(1,2);
    title(strcat( title_name , sprintf(': R = %.2f.' , R)));
    xlabel('Under-replication, %'); ylabel('mean CNR'); xlim([-0.3 0.6]);
end
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/figures/scatter__CNR__VS__dist__G1_and_G1noBrdU' , '-r300');

%% Fig9. boxplots of CNR by under-replication bins
DS_names = {'~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000_A.mat' , ...
    '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000_B.mat' , ...
    '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000_C.mat' , ...
    '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000_D.mat',...
    '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000_E.mat'};
title_names = {'A' , 'B' , 'C' , 'D' , 'E'};
ylims = [1.5 1.5 1.5 1.8];
% Fig4A: G1-M
figure('units','centimeters','position',[5 5 45 15]); 
for Z = 1:length(DS_names)
    DS_name = DS_names{Z}; title_name = title_names{Z};
    load(DS_name);
    DS.underrep_bin = NaN(length(DS) , 1);
    for I = 1:length(DS)
        if DS.percent_underreplicated_cdc20(I) < 0.2
            DS.underrep_bin(I) = 1;
        elseif DS.percent_underreplicated_cdc20(I) < 0.3
            DS.underrep_bin(I) = 2;
        elseif DS.percent_underreplicated_cdc20(I) < 0.4
            DS.underrep_bin(I) = 3;
        elseif DS.percent_underreplicated_cdc20(I) < 0.5
            DS.underrep_bin(I) = 4;
        else
            DS.underrep_bin(I) = 5;
        end
    end
    subplot(1 , length(DS_names) , Z); hold on; grid on; set(gca , 'FontSize' , 12);
    boxplot(DS.mean_CNR_MG1_MM , DS.underrep_bin , 'color' , [.2 .2 .2] , 'symbol' , '');
    ylabel('mean CNR'); set(gca , 'Xtick' , [1:5] , 'XtickLabel' , {'<20' , '20-30', '30-40' , '40-50' , '>50'});
    title(title_name);
    ylim([0.6 ylims(Z)]);
end
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/figures/boxplot__CNR__VS__dist__G1_and_M' , '-r300');

% Fig4B: G1-M
figure('units','centimeters','position',[5 5 45 15]); 
for Z = 2:length(DS_names)
    DS_name = DS_names{Z}; title_name = title_names{Z};
    load(DS_name);
    DS.underrep_bin = NaN(length(DS) , 1);
    for I = 1:length(DS)
        if DS.percent_underreplicated_cdc20(I) < 0.2
            DS.underrep_bin(I) = 1;
        elseif DS.percent_underreplicated_cdc20(I) < 0.3
            DS.underrep_bin(I) = 2;
        elseif DS.percent_underreplicated_cdc20(I) < 0.4
            DS.underrep_bin(I) = 3;
        elseif DS.percent_underreplicated_cdc20(I) < 0.5
            DS.underrep_bin(I) = 4;
        else
            DS.underrep_bin(I) = 5;
        end
    end
    subplot(1 , length(DS_names) , Z); hold on; grid on; set(gca , 'FontSize' , 12);
    boxplot(DS.mean_CNR_MG1_MG1noBrdU , DS.underrep_bin , 'color' , [.2 .2 .2] , 'symbol' , '');
    ylabel('mean CNR'); set(gca , 'Xtick' , [1:5] , 'XtickLabel' , {'<20' , '20-30', '30-40' , '40-50' , '>50'});
    title(title_name);
    ylim([0.6 ylims(Z)]);
end
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/figures/boxplot__CNR__VS__dist__G1_and_G1noBrdU' , '-r300');
    
    
    
    
    
    
    
    
    
    
    

