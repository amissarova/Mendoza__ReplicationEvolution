%% Pipeline to generate panels for fig2 and potential Supplements for fig 2

%% Fig 2A: example for one chromosome, under-replication VS distance to the end (pMet3-cdc20 and dbf2-2)
load('~/Develop/Mendoza__ReplicationEvolution/Data/MutantChr200.mat');
DS = G;
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
xlim([0 K]); ylim([-9 60]);
set(gca , 'Xtick' , [0 25 50 75] , 'Ytick' , [0 20 40 60]);
title(strcat ( DS.chr{idx(1)} , ', left'));
xlabel('Distance to the end, kbp'); ylabel('% unreplicated cells');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2A' , '-r300');

%% Fig 2B: all 200bp-genomic windows: under-replication VS dist to the end + running median and visualization of DM
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');
idx = find(DS.dist_to_the_end_kb <= 250);
DS = DS(idx , :);
clrs2 = hot(12); clrs3 = lines(6); clrs4 = parula(12); clrs5 = summer(12);
figure('units','centimeters','position',[5 5 8 8]); hold on; grid on;

% all 200-bp windows
h1 = scatter(DS.dist_to_the_end_kb , DS.percent_unreplicated_not_trimmed_cdc20_smooth , ...
    5 , clrs5(7,:) , 'filled' , 'MarkerFaceAlpha',.1 , 'Display' , 'all 200-bp windows');

% running median
h2 = plot(DS.dist_to_the_end_kb , DS.median_values_DM , 'LineWidth' , 3.5 , 'color' , clrs5(1,:) , 'Display' , 'running median');

% example for one point above running median w/ high absolute under-rep but
% low DM
idx = find(DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman < 10 & ...
    DS.percent_unreplicated_not_trimmed_cdc20_smooth > 30 & ...
    DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman > 5, 15); idx = idx(end);
% line to the running median itself 
h3 = plot([DS.dist_to_the_end_kb(idx) DS.dist_to_the_end_kb(idx)] , ...
    [DS.median_values_DM(idx) DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx)] , ...
    'LineWidth' , 1.5 , 'color' , clrs2(3,:), 'Display' , sprintf('DM = %.2f' , DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman(idx)));
% a highlighted point
h5 = scatter(DS.dist_to_the_end_kb(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx),...
    30 , clrs2(3,:) , 'filled');
% tick to the line to emphasize that we are measuring length
h6 = plot([DS.dist_to_the_end_kb(idx)-6 DS.dist_to_the_end_kb(idx)+6] , ...
    [DS.median_values_DM(idx) DS.median_values_DM(idx)] , ...
    'LineWidth' , 1.5 , 'color' , clrs2(3,:));
% tick to the line to emphasize that we are measuring length
h7 = plot([DS.dist_to_the_end_kb(idx)-6 DS.dist_to_the_end_kb(idx)+6] , ...
    [DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx)-1 DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx)-1] , ...
    'LineWidth' , 1.5 , 'color' , clrs2(3,:));

% example for one point below running median to illustrate negative DM
idx = find(DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman < -6 & ...
    DS.percent_unreplicated_not_trimmed_cdc20_smooth > -6, 105); idx = idx(end);
% line to the running median itself 
h4 = plot([DS.dist_to_the_end_kb(idx) DS.dist_to_the_end_kb(idx)] , ...
    [DS.median_values_DM(idx) DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx)] , ...
    'LineWidth' , 1.5 , 'color' , clrs3(4,:) ,'Display' , sprintf('DM = %.2f' , DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman(idx)));
% a highlighted point
h8 = scatter(DS.dist_to_the_end_kb(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx),...
    30 , clrs3(4,:) , 'filled');
% tick to the line to emphasize that we are measuring length
h9 = plot([DS.dist_to_the_end_kb(idx)-6 DS.dist_to_the_end_kb(idx)+6] , ...
    [DS.median_values_DM(idx) DS.median_values_DM(idx)] , ...
    'LineWidth' , 1.5 , 'color' , clrs3(4,:));
% tick to the line to emphasize that we are measuring length
h10 = plot([DS.dist_to_the_end_kb(idx)-6 DS.dist_to_the_end_kb(idx)+6] , ...
    [DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx)+1 DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx)+1] , ...
    'LineWidth' , 1.5 , 'color' , clrs3(4,:));
% legend for now is commneted but if you want -- uncomment and specify
% which lines/dots you want to mention in the legend
%legend([h1 h2 h3 h4] ,'location' , 'ne');

ylim([-9 60]); xlim([0 250]);
set(gca , 'Ytick' , [0 20 40 60]);
xlabel('Distance to the end, kbp'); ylabel('% unreplicated cells');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2B' , '-r300');

%% fig 2C (potentially out of the main figure or REDO): R^2 for linear models for :
% dist to the end; + T-rep ; + Ttranscription rate
% +G4, +GC

% this part to generate random split in 20 folds
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');
T = DS;
rand_vec = rand(length(T) , 1);
T.random_num = NaN(length(T) , 1);
for I = 1:length(T)
    T.random_num(I) = rand_vec(I);
end
T = sortrows(T , 'random_num');
K = 20;
N = round(length(T)/K);
all_idx = [1:length(T)];
idx = cell(K , 1);
for I = 1:K-1
    temp_idx = all_idx(1:N);
    all_idx = setdiff(all_idx , temp_idx);
    idx{I} = temp_idx;
end 
idx{K} = all_idx;

% this part to generate linear models and calculate R^2 for test and trains
% samples
% Variable 23 -- log2(distance); Variable 11 -- splined replication timing;
% Variable 37 -- mean Transcription rate; Variable 38 -- G4 value;
% Variable 27 -- GC
% Variable 18 -- under-replication (smoothed by spline)
parameters_idx = cell(5,1);
parameters_idx{1} = [23]; parameters_idx{2} = [23 11]; 
parameters_idx{3} = [23 11 37]; parameters_idx{4} = [23 11 37 38]; parameters_idx{5} = [23 11 37 38 27];
R_train = cell(length(parameters_idx),1);
R_test = cell(length(parameters_idx),1);
for I = 1:length(parameters_idx)
    X = double( T(:,parameters_idx{I}));
    Y = double( T(:, 18)) ;
    temp_train = NaN(K,1);
    temp_test = NaN(K,1);
    for J = 1:K
        idx_train = setdiff([1:length(T)] , idx{J});
        idx_test = idx{J};
        Mdl = fitglm(X(idx_train , :) , Y(idx_train , :) );
        Y_train = predict(Mdl , X(idx_train , :));
        Y_test = predict(Mdl , X(idx_test , :));
        
        mean_y = nanmean(Y(idx_train ));
        SS_tot = sum((Y(idx_train) - mean_y).^2);
        SS_res = sum((Y(idx_train) - Y_train).^2);
        temp_train(J) = 1 - (SS_res/SS_tot);
        
        mean_y = nanmean(Y(idx_test));
        SS_tot = sum((Y(idx_test) - mean_y).^2);
        SS_res = sum((Y(idx_test) - Y_test).^2);
        temp_test(J) = 1 - (SS_res/SS_tot);
    end
    R_train{I} = temp_train;
    R_test{I} = temp_test;
end
% figure itself
legend_titles = {'log (distance to the end)' , '+Replication time' , '+Transcription rate' , '+G4' , '+GC'};
clrs = parula(12); clrs1 = hot(12); clrs2 = summer(12); clrs3 = lines(6); 
clrs_set = [clrs(7,:) ; clrs(10,:) ; clrs(3,:) ; clrs1(3,:) ; clrs2(2,:) ; clrs3(4,:)];
figure('units','centimeters','position',[5 5 8 8]); hold on; grid on; set(gca , 'FontSize' , 10 ); 
data = [R_test{1} R_test{2} R_test{3} R_test{4} R_test{5}];
h1 = boxplot(data , 'color' , [.2 .2 .2], 'symbol','');
h = findobj(gca,'Tag','Box');set(h1 , 'LineWidth' , 1.5);
N = 5;
for j=1:length(h)
    patch(get(h(N-j+1),'XData'),get(h(N-j+1),'YData'), clrs_set(N-j+1,:) ,'FaceAlpha' , .5);
end
legend(legend_titles,'location' , 'SouthOutside');
set(gca , 'Xtick' , []);
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2C' , '-r300');

%% fig 2D: for all 200bp-windows, G4-quadruplexes, absolute under-rep and DM for poor and rich regions

load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');
legend_titles = {'G4-poor' , 'G4-rich'};
clrs1 = lines(6); clrs2 = parula(12); clrs_set = [clrs1(4,:) ; clrs2(6,:)];
fh = figure('units','centimeters','position',[5 5 8 8 ]); 
idx = DS.G4 > 100 ; 
subplot(1,2,2); hold on; grid on; set(gca , 'FontSize' , 10);
h1 = boxplot(DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman , idx ,'symbol','','color' , [.2 .2 .2]);
h = findobj(gca,'Tag','Box');set(h1 , 'LineWidth' , 1.5);
for j=1:length(h)
    patch(get(h(2-j+1),'XData'),get(h(2-j+1),'YData'), clrs_set(2-j+1,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
ylim([-15 20]);
ylabel('under-replication DM');
[~,p] = ttest2(DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman(DS.G4 <= 100) , ...
    DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman(DS.G4 > 100) )
%title(sprintf('P = %.4f.' , p));
set(gca,'xtick',[]);
legend('location', 'SouthOutside');
subplot(1,2,1); hold on; grid on; set(gca , 'FontSize' , 10);
h1 = boxplot(DS.percent_unreplicated_not_trimmed_cdc20_smooth , idx ,'symbol','','color' , [.2 .2 .2]);
h = findobj(gca,'Tag','Box');set(h1 , 'LineWidth' , 1.9);
for j=1:length(h)
    patch(get(h(2-j+1),'XData'),get(h(2-j+1),'YData'), clrs_set(2-j+1,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
[~,p] = ttest2(DS.percent_unreplicated_not_trimmed_cdc20_smooth(DS.G4 <= 100) , ...
    DS.percent_unreplicated_not_trimmed_cdc20_smooth(DS.G4 > 100) )
%title(sprintf('P = %.4f.' , p));
ylim([-15 70]);
ylabel('% unreplicated cells');
set(gca,'xtick',[]);
legend('location', 'SouthOutside');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2D' , '-r300');

%% fig 2E: for all ORFs, transcription rate, absolute under-rep and DM, for low-medium and highly transcribed ORFs

load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features.mat');
idx = find(~isnan(DS.mean_PROseq) & strcmp(DS.TYPE , 'ORF') );
DS = DS(idx , :);
DS.PROseq_bin = NaN(length(DS) , 1);
for I = 1:length(DS)
    if DS.mean_PROseq(I) < 3.2
        DS.PROseq_bin(I) = 1;
    else
        DS.PROseq_bin(I) = 2;
    end
end
figure('units','centimeters','position',[5 5 10 10]);
clrs1 = winter(40);
legend_titles = {'low-medium','high'};
subplot(1,2,1); hold on; grid on; set(gca , 'FontSize' , 10);
h1 = boxplot(DS.percent_unreplicated_not_trimmed_cdc20_smooth , DS.PROseq_bin , 'color' , [.2 .2 .2] , 'symbol' , '');
set(h1 , 'LineWidth' , 1.5);
h = findobj(gca,'Tag','Box'); 
for j=1:2
    patch(get(h(2-j+1),'XData'),get(h(2-j+1),'YData'), clrs1(12*j+4,:) ,...
        'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
set(gca , 'Xtick' , []);
ylabel('% unreplicated cells');
ylim([-15 20]);
[h,p] = ttest2(DS.percent_unreplicated_not_trimmed_cdc20_smooth(DS.PROseq_bin == 1) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(DS.PROseq_bin == 2));
%title(sprintf('P = %.4f.' , p));
legend('location' , 'SouthOutside');

subplot(1,2,2); hold on; grid on; set(gca , 'FontSize' , 10);
h1 = boxplot(DS.percent_unreplicated_not_trimmed_cdc20_DM , DS.PROseq_bin , 'color' , [.2 .2 .2] , 'symbol' , '');
set(h1 , 'LineWidth' , 1.5);
h = findobj(gca,'Tag','Box'); set(h1 , 'LineWidth' , 1.5);
for j=1:2
    patch(get(h(2-j+1),'XData'),get(h(2-j+1),'YData'), clrs1(12*j+4,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
set(gca , 'Xtick' , []); 
ylabel('under-replication DM');
ylim([-15 15]);
legend('location' , 'SouthOutside');
[h,p] = ttest2(DS.percent_unreplicated_not_trimmed_cdc20_DM(DS.PROseq_bin == 1) , DS.percent_unreplicated_not_trimmed_cdc20_DM(DS.PROseq_bin == 2));
%title(sprintf('P = %.4f.' , p));
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2E' , '-r300');

%% Fig 2F: under-replication (absolute and DM): ORF VS tRNA, transposons and fragile sites

load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features.mat');
DS.type_general = cell(length(DS) , 1);
for I = 1:length(DS)
    temp_type = DS.TYPE{I};
    if ~strcmp(temp_type , 'interstitial_del_dup') & ~strcmp(temp_type , 'terminal_del_dup')
        DS.type_general{I} = temp_type;
    else
        DS.type_general{I} = 'del_dup';
    end
end
type_of_interest = {'ORF' , 'tRNA' , 'transposable_element_gene', 'del_dup'};
data = NaN(length(DS) , length(type_of_interest));
data_DM = NaN(length(DS) , length(type_of_interest));
for I = 1:length(type_of_interest)
    idx = find(strcmp(DS.type_general , type_of_interest{I}));
    data(1:length(idx) , I) = DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx);
    data_DM(1:length(idx) , I) = DS.percent_unreplicated_not_trimmed_cdc20_DM(idx);
end
clrs1 = summer(12); clrs2 = hot(12); clrs3 = parula(12); clrs4 = lines(6); clrs5 = spring(12);
clrs_set = [.65 .65 .65 ; clrs1(3,:) ; clrs2(3,:) ; clrs3(5,:) ; ...
    clrs4(4,:) ; clrs3(10,:) ; clrs3(3,:) ; clrs5(3,:)];
legend_titles = {'ORFs' , 'tRNAs' , 'transposable elements' , ...
    'fragile sites'};
figure('units','centimeters','position',[5 5 10 10]); 
subplot(1,2,1); hold on; grid on; set(gca , 'FontSize' , 10);
h1 = boxplot(data , 'color' , [.2 .2 .2] , 'symbol' , '');
set(h1 , 'LineWidth' , 1.5);
h = findobj(gca,'Tag','Box'); 
for j=1:4
    patch(get(h(4-j+1),'XData'),get(h(4-j+1),'YData'), clrs_set(j,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
legend('location' , 'SouthOutside');
ylim([-15 50]);
set(gca , 'Xtick' , '');
ylabel('% unreplicated cells');

subplot(1,2,2); hold on; grid on; set(gca , 'FontSize' , 10);
h1 = boxplot(data_DM , 'color' , [.2 .2 .2] , 'symbol' , '');
set(h1 , 'LineWidth' , 1.5);
h = findobj(gca,'Tag','Box'); 
for j=1:4
    patch(get(h(4-j+1),'XData'),get(h(4-j+1),'YData'), clrs_set(j,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
legend('location' , 'SouthOutside');
ylim([-20 40]);
set(gca , 'Xtick' , '');
ylabel('under-replication DM');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2F' , '-r300');

%% Fig 2G: boxplot for SNPs and InDel frequency (for 200bp-windows) for different under-replication bins

load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');
DS.underrep_bin = NaN(length(DS) , 1);
for I = 1:length(DS)
    if DS.percent_unreplicated_not_trimmed_cdc20_smooth(I) <= 20
        DS.underrep_bin(I) = 1;
    elseif DS.percent_unreplicated_not_trimmed_cdc20_smooth(I) <= 30
        DS.underrep_bin(I) = 2;
    elseif DS.percent_unreplicated_not_trimmed_cdc20_smooth(I) <= 40
        DS.underrep_bin(I) = 3;
    else
        DS.underrep_bin(I) = 4;
    end
end
clrs1 = winter(8);
figure('units','centimeters','position',[5 5 9 9]);
legend_titles = {'<20' , '20-30' , '30-40' , '>40'};
subplot(1,2,1); hold on; grid on; set(gca , 'FontSize' , 12);
h1 = boxplot(DS.freq_SNP , DS.underrep_bin , 'color' , [.2 .2 .2] , 'symbol' , '');
ylim([0 .175]);
set(h1 , 'LineWidth' , 1.5);
h = findobj(gca,'Tag','Box'); set(h1 , 'LineWidth' , 1.5);
for j=1:4
    patch(get(h(4-j+1),'XData'),get(h(4-j+1),'YData'), clrs1(2*j,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
set(gca , 'Xtick' , []); 
%xlabel('% unreplicated cells');
ylabel('SNP frequency');
legend('location' , 'SouthOutside');

subplot(1,2,2); hold on; grid on; set(gca , 'FontSize' , 12);
h1 = boxplot(DS.freq_indel , DS.underrep_bin , 'color' , [.2 .2 .2] , 'symbol' , '');
ylim([0 .175]);
set(h1 , 'LineWidth' , 1.5);
h = findobj(gca,'Tag','Box'); set(h1 , 'LineWidth' , 1.5);
for j=1:4
    patch(get(h(4-j+1),'XData'),get(h(4-j+1),'YData'), clrs1(2*j,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
set(gca , 'Xtick' , []); 
%xlabel('% unreplicated cells');
ylabel('InDel frequency');
legend('location' , 'SouthOutside');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2G' , '-r300');

%% Fig 2H: boxplot for gene preservation (for ORFs) for different under-replication bins

load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features.mat');
idx = find(~isnan(DS.frac_conservation)); DS = DS(idx , :);
DS.underrep_bin = NaN(length(DS) , 1);
for I = 1:length(DS)
    if DS.percent_unreplicated_not_trimmed_cdc20_smooth(I) <= 20
        DS.underrep_bin(I) = 1;
    elseif DS.percent_unreplicated_not_trimmed_cdc20_smooth(I) <= 30
        DS.underrep_bin(I) = 2;
    elseif DS.percent_unreplicated_not_trimmed_cdc20_smooth(I) <= 40
        DS.underrep_bin(I) = 3;
    else
        DS.underrep_bin(I) = 4;
    end
end
clrs1 = winter(8);
figure('units','centimeters','position',[5 5 9 9]);
legend_titles = {'<20' , '20-30' , '30-40' , '>40'};
hold on; grid on; set(gca , 'FontSize' , 12);
h1 = boxplot(DS.frac_conservation*100 , DS.underrep_bin , 'color' , [.2 .2 .2] , 'symbol' , '');
ylim([80 100]);
set(h1 , 'LineWidth' , 1.5);
h = findobj(gca,'Tag','Box'); set(h1 , 'LineWidth' , 1.5);
for j=1:4
    patch(get(h(4-j+1),'XData'),get(h(4-j+1),'YData'), clrs1(2*j,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
set(gca , 'Xtick' , []); 
%xlabel('% unreplicated cells');
%ylabel('% strains with preserved gene');
legend('location' , 'SouthOutside');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2H' , '-r300');

%% FigSupAlsu1: all chromosomes, all repeatable replicates: under-replication VS position on genome

load('~/Develop/Mendoza__ReplicationEvolution/Data/MutantChr200.mat');
DS = G;
unq_mutant = {'cdc20'}; chr_nums = [1:16];
clrs1 = hot(12); clrs2 = lines(6); clrs3 = parula(12); clrs4 = summer(12); clrs5 = spring(12);
clrs_set = [clrs1(3,:) ; clrs2(4,:) ; clrs3(2,:) ; clrs3(10,:) ; clrs4(2,:) ; clrs3(6,:)];
for I = 1:length(unq_mutant)
    figure('units','centimeters','position',[5 5 15 15]);
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
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/SupAlsu1' , '-r300');


%% FigSupAlsu2: only for pMet3cdc20 -- histogram of length of under-replication and % of unreplicated at the last kb for all chromosomal parts

load('~/Develop/Mendoza__ReplicationEvolution/Data/MutantChr200.mat');
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
h = histogram(data , [0:2:40] );
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

print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/SupAlsu2' , '-r300');

%% Fig SupAlsu3: boxplot of under-replication and under-replication DM by replication timing bins

load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');

% replication bins
trep_grid = [0 40 45 50 70];
legend_titles = {'<= 40 min' , '40-45 min' , '45-50 min' , '> 50 min'};

figure ('units','centimeters','position',[5 5 12 12]);
subplot(1,2,2); hold on; grid on; set(gca , 'FontSize' , 10);
data = NaN(length(DS) , length(trep_grid)-1);
clrs = winter(4*length(trep_grid)-4);
for I = 1:length(trep_grid)-1
    idx = find(DS.Trep_spline > trep_grid(I) & DS.Trep_spline < trep_grid(I+1) );
    data(1:length(idx) , I) = DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman(idx);
end
h1 = boxplot(data , 'color' , [.2 .2 .2]  , 'symbol' , '');
h = findobj(gca,'Tag','Box'); set(h1 , 'LineWidth' , 1.5);
for j=1:4
    patch(get(h(4-j+1),'XData'),get(h(4-j+1),'YData'), clrs(4*j-2,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end

ylim([-20 20]);
set(gca , 'Xtick' , []);

ylabel('under-replication DM');
legend('location' , 'SouthOutside');

subplot(1,2,1); hold on; grid on; set(gca , 'FontSize' , 10);
data = NaN(length(DS) , length(trep_grid)-1);
clrs = winter(4*length(trep_grid)-4);
for I = 1:length(trep_grid)-1
    idx = find(DS.Trep_spline > trep_grid(I) & DS.Trep_spline < trep_grid(I+1) );
    data(1:length(idx) , I) = DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx);
end
h1 = boxplot(data , 'color' , [.2 .2 .2]  , 'symbol' , '');
h = findobj(gca,'Tag','Box'); set(h1 , 'LineWidth' , 1.5);
for j=1:4
    patch(get(h(4-j+1),'XData'),get(h(4-j+1),'YData'), clrs(4*j-2,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
set(h , 'LineWidth' , 1.5);
ylim([-20 70]);
set(gca , 'Xtick' , []);
legend('location' , 'SouthOutside');
ylabel('% unreplicated cells');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/SupAlsu3' , '-r300');


%% Fig SupAlsu4: scatterplot, x-axis -- Replication timing, y-axis -- under-replication

load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');
clrs = parula(12);
figure('units','centimeters','position',[5 5 10 10]); hold on; grid on;
scatter(DS.Trep_spline , DS.percent_unreplicated_not_trimmed_cdc20_smooth , 20 , clrs(5,:) , 'filled' , 'MarkerFaceAlpha' , .1)
xlim([10 65]);
ylim([-20 70]);
R = corrcoef(DS.Trep_spline , DS.percent_unreplicated_not_trimmed_cdc20_smooth);
title(sprintf('R = %.1f.' , R(1,2)));
xlabel('Replication time, minutes');
ylabel('% unreplicated cells');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/SupAlsu4' , '-r300');

%% Fig SupAlsu5: boxplot for SNP and Indel frequency (for 200bp-windows) for different under-replication DM bins

load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');
DS.underrep_DM_bin = NaN(length(DS) , 1);
for I = 1:length(DS)
    if DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman(I) <= 0
        DS.underrep_DM_bin(I) = 1;
    else 
        DS.underrep_DM_bin(I) = 2;
    end
end
clrs1 = autumn(8);
figure('units','centimeters','position',[5 5 9 9]);
legend_titles = {'<0' , '>0'};
subplot(1,2,1); hold on; grid on; set(gca , 'FontSize' , 12);
h1 = boxplot(DS.freq_SNP , DS.underrep_DM_bin , 'color' , [.2 .2 .2] , 'symbol' , '');
ylim([0 .1]);
set(h1 , 'LineWidth' , 1.5);
h = findobj(gca,'Tag','Box'); set(h1 , 'LineWidth' , 1.5);
for j=1:2
    patch(get(h(2-j+1),'XData'),get(h(2-j+1),'YData'), clrs1(3*j,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
set(gca , 'Xtick' , []); 
%xlabel('under-replication DM');
ylabel('SNP frequency');
legend('location' , 'SouthOutside');

subplot(1,2,2); hold on; grid on; set(gca , 'FontSize' , 12);
h1 = boxplot(DS.freq_indel , DS.underrep_DM_bin , 'color' , [.2 .2 .2] , 'symbol' , '');
ylim([0 .06]);
set(h1 , 'LineWidth' , 1.5);
h = findobj(gca,'Tag','Box'); set(h1 , 'LineWidth' , 1.5);
for j=1:2
    patch(get(h(2-j+1),'XData'),get(h(2-j+1),'YData'), clrs1(3*j,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
set(gca , 'Xtick' , []); 
%xlabel('under-replication DM');
ylabel('InDel frequency');
legend('location' , 'SouthOutside');

print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/SupAlsu5' , '-r300');



















