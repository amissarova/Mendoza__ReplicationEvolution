%% generate DS with coverage info for all samples
% run once
cd ~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq
addpath(genpath('~/Develop/matlab'));
DS = dataset();

folder_names = {'~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/coverage_data_1000/G1_HU',...
    '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/coverage_data_1000/M_G1',...
    '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/coverage_data_1000/M_M',...
    '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/coverage_data_1000/M_G1__noBrdU'};
folder_IDs = {'G1HU' , 'MG1' , 'MM' , 'MG1noBrdU'};
for Z = 1:length(folder_names)
    folder_name = folder_names{Z}; folder_ID = folder_IDs{Z};
    samples = dir(folder_name); samples = samples(3:end);
    cd(folder_name)
    for I = 1:length(samples)
        D = dataset('file' , samples(I).name);
        D = D(: , [1:4]);
        D.CN = NaN(length(D) , 1);
        D.pct = NaN(length(D) , 1);
        m = modefit( D.cov , 0 , [-5:5:1000] );
        if m > 0
            D.CN = D.cov/m;
            for J = 1:length(D)
                D.pct(J) = nanmean(D.cov < D.cov(J));
            end
        end
        
        D.Properties.VarNames{4} = strcat('cov_' , folder_ID , '_' , sprintf('%d',I));
        D.Properties.VarNames{5} = strcat('CN_' , folder_ID , '_' , sprintf('%d',I));
        D.Properties.VarNames{6} = strcat('pct_' , folder_ID , '_' , sprintf('%d',I));
        if isempty(DS)
            DS = D;
        else
            DS = join(DS , D , 'Type' , 'Left' , 'Keys' , {'chr' , 'start_point' , 'end_point'} , 'MergeKeys' , true);
        end
    end
end
% keep only gDNA chromosomes
idx = [];
for I = 1:length(DS)
    if findstr(DS.chr{I} , 'chr') & ~strcmp(DS.chr{I} , 'chrMito')
        idx = [idx I];
    end
end
DS = DS(idx,:);
% add info from underreplication dataset regarding each window (dist to the end, RT, dist to ARS)
D = DS;
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp_new.mat');
DS = DS(: , {'chr' , 'start_point' , 'end_point' , ...
    'Trep_spline' , 'dist_to_ARS' , 'dist_to_the_end_kb' , ...
    'percent_underreplicated_cdc20' , 'percent_underreplicated_dbf2'});
D.dist_to_the_end_kb = NaN(length(D) , 1);
D.Trep_spline = NaN(length(D) , 1);
D.dist_to_ARS = NaN(length(D) , 1);
D.percent_underreplicated_cdc20 = NaN(length(D) , 1);
D.percent_underreplicated_dbf2 = NaN(length(D) , 1);

for I = 1:length(D)
    idx = find(DS.start_point >= D.start_point(I) & DS.start_point < D.end_point(I) & strcmp(DS.chr , D.chr{I}));
    D.dist_to_the_end_kb(I) = nanmean(DS.dist_to_the_end_kb(idx));
    D.Trep_spline(I) = nanmean(DS.Trep_spline(idx));
    D.dist_to_ARS(I) = nanmean(DS.dist_to_ARS(idx));
    D.percent_underreplicated_cdc20(I) = nanmean(DS.percent_underreplicated_cdc20(idx));
    D.percent_underreplicated_dbf2(I) = nanmean(DS.percent_underreplicated_dbf2(idx));
end
DS = D;
% save DS
cd ~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq
save('~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000.mat' , 'DS');

%% Fig1. Check the coverage distribution for all the samples
load('~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000.mat');
data = NaN(length(DS) , 6+5+3);
data(:,1) = DS.cov_MG1_1; data(:,2) = DS.cov_MG1_2; data(:,3) = DS.cov_MG1_3; data(:,4) = DS.cov_MG1_4; data(:,5) = DS.cov_MG1_5; data(:,6) = DS.cov_MG1_6; 
data(:,7) = DS.cov_MM_1; data(:,8) = DS.cov_MM_2; data(:,9) = DS.cov_MM_3; data(:,10) = DS.cov_MM_4; data(:,11) = DS.cov_MM_5;
data(:,12) = DS.cov_MG1noBrdU_1; data(:,13) = DS.cov_MG1noBrdU_2; data(:,14) = DS.cov_MG1noBrdU_3; 

figure('units','centimeters','position',[5 5 45 20]); hold on; grid on; set(gca , 'FontSize' , 12);
h = boxplot(data , 'color' , [.2 .2 .2] , 'symbol' , '');
ylim([0 600]); ylabel('Coverage');
set(gca , 'Xtick' , [1:14] , 'XtickLabel' , {'G1-1','G1-2','G1-3','G1-4','G1-5','G1-6',...
    'M-1','M-2','M-3','M-4','M-5' , 'noBr-1' , 'noBr-2' , 'noBr-3'});
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/figures/boxplot__coverage__G1_and_M' , '-r300');

%% Fig 2. Coverage across distance to the end for all samples in M-G1, M-M and M-G1noBrdU
load('~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000.mat');
data = NaN(length(DS) , 6+5+3);
data(:,1) = DS.cov_MG1_1; data(:,2) = DS.cov_MG1_2; data(:,3) = DS.cov_MG1_3; data(:,4) = DS.cov_MG1_4; data(:,5) = DS.cov_MG1_5; data(:,6) = DS.cov_MG1_6; 
data(:,7) = DS.cov_MM_1; data(:,8) = DS.cov_MM_2; data(:,9) = DS.cov_MM_3; data(:,10) = DS.cov_MM_4; data(:,11) = DS.cov_MM_5;
data(:,12) = DS.cov_MG1noBrdU_1; data(:,13) = DS.cov_MG1noBrdU_2; data(:,14) = DS.cov_MG1noBrdU_3; 

title_names = {'G1-1','G1-2','G1-3','G1-4','G1-5','G1-6','M-1','M-2','M-3','M-4','M-5','noBr-1' , 'noBr-2' , 'noBr-3'};
ylims = [700 150 300 700 300 300 300 300 30 700 300 200 300 350];
figure('units','centimeters','position',[5 5 45 20]); 
clrs = winter(12);
for I = 1:14
    subplot(4,4,I);hold on; grid on; set(gca , 'FontSize' , 12);
    scatter(DS.dist_to_the_end_kb , data(:,I) , 15 , clrs(7,:) , 'filled');
    xlabel('Distance to the end'); ylabel('Coverage'); title(title_names{I});
    ylim([0 ylims(I)]);
end
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/figures/scatter__coverage__VS__dist__G1_and_M' , '-r300');

%% Analyses 2. Add CN-ratio across replicates
% 4 options -- 
% a) only account samples with median coverage > 300 (G1-1, G1-4, M-4 );
% b) only account samples with median coverage > 100 (G1-1, G1-4, G1-5, G1-6, M-1, M-4, M-5, noBr-2 , noBr-3)
% c) account everything, except M3 which has average coverage ~5
% d) account only matched up by date(i.e. coming from the same sample)
% replicates; for G1-M: G1-3 -- M-1; G1-4 -- M-2; exclude G1-5 -- M-3 since M-3 has
% too low coverage; for G1-G1noBrdU: G1-4 -- noBr1; G1-6 -- noBr2

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


%% Fig3. scatter of CNR against distance to the end for datasets
clrs = winter(12);
DS_names = {'~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000_A.mat' , ...
    '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000_B.mat' , ...
    '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000_C.mat' , ...
    '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000_D.mat'};
title_names = {'A' , 'B' , 'C' , 'D'};
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

%% Fig4. boxplots of CNR by under-replication bins
DS_names = {'~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000_A.mat' , ...
    '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000_B.mat' , ...
    '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000_C.mat' , ...
    '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000_D.mat'};
title_names = {'A' , 'B' , 'C' , 'D'};
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
    
    
    
    
    
    
    
    
    
    
    

