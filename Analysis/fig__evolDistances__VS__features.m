%%

load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');

idx = find(~isnan(DS.percent_unreplicated_not_trimmed_cdc20) & ~isnan(DS.percent_unreplicated_not_trimmed_dbf2) );
T = DS(idx , :);
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

%% Freq SNPs: dist, Trep, dist+Trep; underrep
DS.dist_log = log2(DS.dist_to_the_end_kb);
T.dist_log = log2(T.dist_to_the_end_kb);
%%
parameters_idx = cell(4,1);
parameters_idx{1} = [42]; parameters_idx{2} = [13]; parameters_idx{3} = [42 13]; parameters_idx{4} = [25]; 
R_train = cell(length(parameters_idx),1);
R_test = cell(length(parameters_idx),1);
for I = 1:length(parameters_idx)
    X = double( T(:,parameters_idx{I}));
    Y = double( T(:, 21)) ;
    temp_train = NaN(K,1);
    temp_test = NaN(K,1);
    for J = 1:K
        %aa = [I J]
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

figure; hold on; grid on; set(gca , 'FontSize' , 20); 
data = [R_test{1} R_test{2} R_test{3} R_test{4} ];
boxplot(data , 'notch' , 'on' , 'color' , [.2 .2 .2], 'symbol','');

clrs = parula(16);
%clrs_set = [clrs2(2,:); clrs2(10,:); clrs1(3,:) ; clrs4(2,:); clrs3(4,:); clrs5(3,:); clrs2(6,:); clrs1(3,:)];
%labels = {'dist' , 'dist+Trep' , 'dist+chr', 'dist+chr+Trep' , 'dist+ARS' , 'dist+expression+Trep+chr' , 'all'};
for I = 1:length(parameters_idx)
    scatter(repmat(I,K,1)+rand(K,1)*0.2-0.1 , R_test{I} , 100 , clrs(4*I-2,:) , 'filled' , 'd' , 'MarkerEdgeColor' , [.2 .2 .2]);
%     [~ , p] = ttest2(R_test{I} , R_train{I});
%     if p < 0.05
%         scatter(3*I-1.5 , 1 , 100 , 'k'  , 'x');
%     end
end
xticklabels = {'1.distance to telomer, log' , '2.Replication Timing' , '3. 1+2' , '4.% unreplicated in M' };
%set(gca , 'Xtick' , [1:8] , 'XtickLabel' , xticklabels);

ylabel('R^{2}');
%xtickangle(-45);
set(gcf , 'PaperPosition' , [0 0 20 20]);
title('% of strains, where gene is preserved'); 
print('-dpsc2' , 'frac_conservation__R2.eps');





%% Freq SNPs
parameters_idx = cell(6,1);
parameters_idx{1} = [13]; parameters_idx{2} = [11]; parameters_idx{3} = [9]; parameters_idx{4} = [13 11 9]; 
parameters_idx{5} = [19]; parameters_idx{6} = [18]; 
%
R_train = cell(length(parameters_idx),1);
R_test = cell(length(parameters_idx),1);
for I = 1:length(parameters_idx)
    X = double( T(:,parameters_idx{I}));
    Y = double( T(:, 16 )) ;
    temp_train = NaN(K,1);
    temp_test = NaN(K,1);
    for J = 1:K
        %aa = [I J]
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

figure; hold on; grid on; set(gca , 'FontSize' , 20); 
data = [R_test{1} R_test{2} R_test{3} R_test{4} R_test{5} R_test{6}];
boxplot(data , 'notch' , 'on' , 'color' , [.2 .2 .2]);

clrs = parula(24);
%clrs_set = [clrs2(2,:); clrs2(10,:); clrs1(3,:) ; clrs4(2,:); clrs3(4,:); clrs5(3,:); clrs2(6,:); clrs1(3,:)];
%labels = {'dist' , 'dist+Trep' , 'dist+chr', 'dist+chr+Trep' , 'dist+ARS' , 'dist+expression+Trep+chr' , 'all'};
for I = 1:length(parameters_idx)
    scatter(repmat(I,K,1)+rand(K,1)*0.2-0.1 , R_test{I} , 70 , clrs(4*I-2,:) , 'filled' , 'd');
%     [~ , p] = ttest2(R_test{I} , R_train{I});
%     if p < 0.05
%         scatter(3*I-1.5 , 1 , 100 , 'k'  , 'x');
%     end
end
%xticklabels = {'1.dist to ARS' , '1.Replication Timing' , '2.Test' , '2.Train' , '3.Test' , '3.Train' , '4.Test' , '4.Train'};
set(gca , 'Xtick' , [1:8] , 'XtickLabel' , xticklabels);

ylabel('R^{2}');
%xtickangle(-45);
set(gcf , 'PaperPosition' , [0 0 20 20]);
title('SNP frequency');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Manuscript/Fig11A');

%% Freq INDELs
parameters_idx = cell(6,1);
parameters_idx{1} = [13]; parameters_idx{2} = [11]; parameters_idx{3} = [9]; parameters_idx{4} = [13 11 9]; 
parameters_idx{5} = [19]; parameters_idx{6} = [18]; 
%
R_train = cell(length(parameters_idx),1);
R_test = cell(length(parameters_idx),1);
for I = 1:length(parameters_idx)
    X = double( T(:,parameters_idx{I}));
    Y = double( T(:, 22 )) ;
    temp_train = NaN(K,1);
    temp_test = NaN(K,1);
    for J = 1:K
        %aa = [I J]
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

figure; hold on; grid on; set(gca , 'FontSize' , 20); 
data = [R_test{1} R_test{2} R_test{3} R_test{4} R_test{5} R_test{6}];
boxplot(data , 'notch' , 'on' , 'color' , [.2 .2 .2]);

clrs = parula(24);
for I = 1:length(parameters_idx)
    scatter(repmat(I,K,1)+rand(K,1)*0.2-0.1 , R_test{I} , 70 , clrs(4*I-2,:) , 'filled' , 'd');
end
set(gca , 'Xtick' , [1:8] , 'XtickLabel' , xticklabels);

ylabel('R^{2}');
set(gcf , 'PaperPosition' , [0 0 20 20]);
title('InDel frequency');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Manuscript/Fig11B');

%%

load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features.mat');
idx = find(~isnan(DS.percent_unreplicated_not_trimmed_cdc20) & ~isnan(DS.percent_unreplicated_not_trimmed_dbf2) & ...
    ~isnan(DS.frac_conservation));
T = DS(idx , :);
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

%% Preservation frac
parameters_idx = cell(6,1);
parameters_idx{1} = [23]; parameters_idx{2} = [16]; parameters_idx{3} = [14]; parameters_idx{4} = [14 16 23]; 
parameters_idx{5} = [11]; parameters_idx{6} = [10]; 
%
R_train = cell(length(parameters_idx),1);
R_test = cell(length(parameters_idx),1);
for I = 1:length(parameters_idx)
    X = double( T(:,parameters_idx{I}));
    Y = double( T(:, 22 )) ;
    temp_train = NaN(K,1);
    temp_test = NaN(K,1);
    for J = 1:K
        %aa = [I J]
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

figure; hold on; grid on; set(gca , 'FontSize' , 20); 
data = [R_test{1} R_test{2} R_test{3} R_test{4} R_test{5} R_test{6}];
boxplot(data , 'notch' , 'on' , 'color' , [.2 .2 .2]);

clrs = parula(24);
for I = 1:length(parameters_idx)
    scatter(repmat(I,K,1)+rand(K,1)*0.2-0.1 , R_test{I} , 70 , clrs(4*I-2,:) , 'filled' , 'd');
end
set(gca , 'Xtick' , [1:6] , 'XtickLabel' , xticklabels);

ylabel('R^{2}');
set(gcf , 'PaperPosition' , [0 0 20 20]);
title('Fraction of strains containing a gene');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Manuscript/Fig11C');


%% legend for 11A-C
figure; hold on; grid on; set(gca , 'FontSize' , 20);
clrs = parula(16);
for I = 1:length(parameters_idx)
    scatter(repmat(I,K,1)+rand(K,1)*0.2-0.1 , R_test{I} , 70 , clrs(4*I-2,:) , 'filled' , 'd');
end
legend('1. Distance to telomeric end, log' , '2. Replication timing' , '3. Distance to telomeric end, log + Replication timing'  , ...
    '4. Degree of underreplication in metaphase' ,'location' , 'best');




%%
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');
DS.unreplicated_cdc20_bin = NaN(length(DS) , 1);
for I = 1:length(DS)
    if DS.percent_unreplicated_not_trimmed_cdc20_smooth(I) <= 20
        DS.unreplicated_cdc20_bin(I) = 1;
    elseif DS.percent_unreplicated_not_trimmed_cdc20_smooth(I) <= 40
        DS.unreplicated_cdc20_bin(I) = 2;
    else
        DS.unreplicated_cdc20_bin(I) = 3;
    end
end

%% SNP -- Trep
DS = sortrows(DS , 'Trep_spline');
unq_underrep_bin = unique(DS.unreplicated_cdc20_bin);
rt = NaN(length(DS) , length(unq_underrep_bin));
data = NaN(length(DS) , length(unq_underrep_bin));
K = 500;
for I = 1:length(unq_underrep_bin)
    idx = find(DS.unreplicated_cdc20_bin == unq_underrep_bin(I));
    D = DS(idx , :);
    temp_data = NaN(length(D)-K+1 , 1);
    temp_rt = NaN(length(D)-K+1 , 1);
    for J = 1:length(D)-K+1
        temp_data(J) = nanmean( D.freq_SNP(J:J+K-1));
        temp_rt(J) = nanmean( D.Trep_spline(J:J+K-1));
    end
    data(1:length(D)-K+1 , I) = temp_data;
    rt(1:length(D)-K+1 , I) = temp_rt;
end

clrs = parula(12);
figure; hold on; grid on; set(gca , 'FontSize' , 14);
for I = 1:length(unq_underrep_bin)
    plot(rt(:,I) , data(:, I) , 'LineWidth' , 3 , 'color' , clrs(4*I-2,:));
end
xlabel('Replication timing, min after \alpha factor release');
ylabel('SNP frequency');
set(gcf , 'PaperPosition' , [0 0 10 10]);
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Manuscript/Fig12_SNP_RT');

%% INDEL -- Trep
DS = sortrows(DS , 'Trep_spline');
unq_underrep_bin = unique(DS.unreplicated_cdc20_bin);
rt = NaN(length(DS) , length(unq_underrep_bin));
data = NaN(length(DS) , length(unq_underrep_bin));
K = 500;
for I = 1:length(unq_underrep_bin)
    idx = find(DS.unreplicated_cdc20_bin == unq_underrep_bin(I));
    D = DS(idx , :);
    temp_data = NaN(length(D)-K+1 , 1);
    temp_rt = NaN(length(D)-K+1 , 1);
    for J = 1:length(D)-K+1
        temp_data(J) = nanmean( D.freq_indel(J:J+K-1));
        temp_rt(J) = nanmean( D.Trep_spline(J:J+K-1));
    end
    data(1:length(D)-K+1 , I) = temp_data;
    rt(1:length(D)-K+1 , I) = temp_rt;
end
%
clrs = parula(12);
figure; hold on; grid on; set(gca , 'FontSize' , 14);
for I = 1:length(unq_underrep_bin)
    plot(rt(:,I) , data(:, I) , 'LineWidth' , 3 , 'color' , clrs(4*I-2,:));
end
xlabel('Replication timing, min after \alpha factor release');
ylabel('Indel frequency');
set(gcf , 'PaperPosition' , [0 0 10 10]);
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Manuscript/Fig12_indel_RT');

%% SNP -- Dist
DS = sortrows(DS , 'dist_to_the_end_kb');
unq_underrep_bin = unique(DS.unreplicated_cdc20_bin);
dist_to_the_end = NaN(length(DS) , length(unq_underrep_bin));
data = NaN(length(DS) , length(unq_underrep_bin));
K = 500;
for I = 1:length(unq_underrep_bin)
    idx = find(DS.unreplicated_cdc20_bin == unq_underrep_bin(I) & DS.dist_to_the_end_kb <= 100);
    D = DS(idx , :);
    temp_data = NaN(length(D)-K+1 , 1);
    temp_dist = NaN(length(D)-K+1 , 1);
    for J = 1:length(D)-K+1
        temp_data(J) = nanmean( D.freq_SNP(J:J+K-1));
        temp_dist(J) = nanmean( D.dist_to_the_end_kb(J:J+K-1));
    end
    data(1:length(D)-K+1 , I) = temp_data;
    dist_to_the_end(1:length(D)-K+1 , I) = temp_dist;
end
%
clrs = parula(12);
figure; hold on; grid on; set(gca , 'FontSize' , 14);
for I = 1:length(unq_underrep_bin)
    plot(dist_to_the_end(:,I) , data(:, I) , 'LineWidth' , 3 , 'color' , clrs(4*I-2,:));
end
xlabel('Distance to telomeres KB');
ylabel('SNP frequency');
set(gcf , 'PaperPosition' , [0 0 10 10]);
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Manuscript/Fig12_SNP_dist');

%% Indel - dist
DS = sortrows(DS , 'dist_to_the_end_kb');
unq_underrep_bin = unique(DS.unreplicated_cdc20_bin);
dist_to_the_end = NaN(length(DS) , length(unq_underrep_bin));
data = NaN(length(DS) , length(unq_underrep_bin));
K = 500;
for I = 1:length(unq_underrep_bin)
    idx = find(DS.unreplicated_cdc20_bin == unq_underrep_bin(I) & DS.dist_to_the_end_kb <= 100);
    D = DS(idx , :);
    temp_data = NaN(length(D)-K+1 , 1);
    temp_dist = NaN(length(D)-K+1 , 1);
    for J = 1:length(D)-K+1
        temp_data(J) = nanmean( D.freq_indel(J:J+K-1));
        temp_dist(J) = nanmean( D.dist_to_the_end_kb(J:J+K-1));
    end
    data(1:length(D)-K+1 , I) = temp_data;
    dist_to_the_end(1:length(D)-K+1 , I) = temp_dist;
end
clrs = parula(12);
figure; hold on; grid on; set(gca , 'FontSize' , 14);
for I = 1:length(unq_underrep_bin)
    plot(dist_to_the_end(:,I) , data(:, I) , 'LineWidth' , 3 , 'color' , clrs(4*I-2,:));
end
xlabel('Distance to telomeres KB');
ylabel('Indel frequency');
set(gcf , 'PaperPosition' , [0 0 10 10]);
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Manuscript/Fig12_indel_dist');

%%
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features.mat');
DS.unreplicated_cdc20_bin = NaN(length(DS) , 1);
for I = 1:length(DS)
    if DS.percent_unreplicated_not_trimmed_cdc20(I) <= 20
        DS.unreplicated_cdc20_bin(I) = 1;
    elseif DS.percent_unreplicated_not_trimmed_cdc20(I) <= 40
        DS.unreplicated_cdc20_bin(I) = 2;
    else
        DS.unreplicated_cdc20_bin(I) = 3;
    end
end
%% gene conservation -- Trep
DS = sortrows(DS , 'Trep_spline');
unq_underrep_bin = unique(DS.unreplicated_cdc20_bin);
rt = NaN(length(DS) , length(unq_underrep_bin));
data = NaN(length(DS) , length(unq_underrep_bin));
K = 50;
for I = 1:length(unq_underrep_bin)
    idx = find(DS.unreplicated_cdc20_bin == unq_underrep_bin(I) & ~isnan(DS.frac_conservation));
    D = DS(idx , :);
    temp_data = NaN(length(D)-K+1 , 1);
    temp_rt = NaN(length(D)-K+1 , 1);
    for J = 1:length(D)-K+1
        temp_data(J) = nanmean( D.frac_conservation(J:J+K-1) < 1);
        temp_rt(J) = nanmean( D.Trep_spline(J:J+K-1));
    end
    data(1:length(D)-K+1 , I) = temp_data;
    rt(1:length(D)-K+1 , I) = temp_rt;
end

clrs = parula(12);
figure; hold on; grid on; set(gca , 'FontSize' , 14);
for I = 1:length(unq_underrep_bin)
    plot(rt(:,I) , data(:, I)*100 , 'LineWidth' , 3 , 'color' , clrs(4*I-2,:));
end
xlabel('Replication timing, min after \alpha factor release');
ylabel('% of genes, that not presented in all species');
set(gcf , 'PaperPosition' , [0 0 10 10]);
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Manuscript/Fig12_fracConservation_RT');

%%
DS = sortrows(DS , 'dist_to_the_end_kb');
unq_underrep_bin = unique(DS.unreplicated_cdc20_bin);
dist_to_the_end = NaN(length(DS) , length(unq_underrep_bin));
data = NaN(length(DS) , length(unq_underrep_bin));
K = 50;
for I = 1:length(unq_underrep_bin)
    idx = find(DS.unreplicated_cdc20_bin == unq_underrep_bin(I) & DS.dist_to_the_end_kb <= 100 & ~isnan(DS.frac_conservation));
    D = DS(idx , :);
    temp_data = NaN(length(D)-K+1 , 1);
    temp_dist = NaN(length(D)-K+1 , 1);
    for J = 1:length(D)-K+1
        temp_data(J) = nanmean( D.frac_conservation(J:J+K-1) < 1);
        temp_dist(J) = nanmean( D.dist_to_the_end_kb(J:J+K-1));
    end
    data(1:length(D)-K+1 , I) = temp_data;
    dist_to_the_end(1:length(D)-K+1 , I) = temp_dist;
end
clrs = parula(12);
figure; hold on; grid on; set(gca , 'FontSize' , 14);
for I = 1:length(unq_underrep_bin)
    plot(dist_to_the_end(:,I) , data(:, I)*100 , 'LineWidth' , 3 , 'color' , clrs(4*I-2,:));
end
xlabel('Distance to telomeres KB');
ylabel('% of genes, that not presented in all species');
set(gcf , 'PaperPosition' , [0 0 10 10]);
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Manuscript/Fig12_fracConservation_dist');

%% legend
figure; hold on; grid on; set(gca , 'FontSize' , 20);
for I = 1:length(unq_underrep_bin)
    plot(dist_to_the_end(:,I) , data(:, I)*100 , 'LineWidth' , 3 , 'color' , clrs(4*I-2,:));
end
legend('<20%' , '20-40%' , '>40%');





%%
labels = {'%unrep: < 10%' , '%unrep: 10-25' , '%unrep: 25-40%' , '%unrep: > 40%'};
figure; 
unq_underrep_bin = unique(DS.unreplicated_cdc20_bin);
clrs = parula(4*length(unq_underrep_bin));
dist_grid = [0:5:100];
rep_time_grid = [10:2:70];

subplot(2,2,1); hold on; grid on; set(gca , 'FontSize' , 12);
for I = 1:length(unq_underrep_bin)
    data = NaN(length(dist_grid) - 1 , 1);
    for J = 1:length(dist_grid)-1
        idx = find(DS.unreplicated_cdc20_bin == unq_underrep_bin(I) &...
            DS.dist_to_the_end_kb >= dist_grid(J) & DS.dist_to_the_end_kb < dist_grid(J+1));
        data(J) = nanmean(DS.freq_SNP(idx));
    end
    idx = find(~isnan(data));
    plot(dist_grid(idx) , data(idx) , 'color', clrs(4*I-2,:) , 'LineWidth' , 2);
    scatter(dist_grid(idx) , data(idx) , 20 , clrs(4*I-2,:) , 'filled');
end
ylabel('SNP freq');
xlabel('Distance to the end');    

subplot(2,2,2); hold on; grid on; set(gca , 'FontSize' , 12);
for I = 1:length(unq_underrep_bin)
    data = NaN(length(rep_time_grid) - 1 , 1);
    for J = 1:length(rep_time_grid)-1
        idx = find(DS.unreplicated_cdc20_bin == unq_underrep_bin(I) &...
            DS.Trep_spline >= rep_time_grid(J) & DS.Trep_spline < rep_time_grid(J+1));
        data(J) = nanmean(DS.freq_SNP(idx));
    end
    idx = find(~isnan(data));
    
    plot(rep_time_grid(idx) , data(idx) , 'color', clrs(4*I-2,:) , 'LineWidth' , 2);
    scatter(rep_time_grid(idx) , data(idx) , 20 , clrs(4*I-2,:) , 'filled');
end
ylabel('SNP freq');
xlabel('Replication timing');    
xlim([0 65]);

subplot(2,2,3); hold on; grid on; set(gca , 'FontSize' , 12);
for I = 1:length(unq_underrep_bin)
    data = NaN(length(dist_grid) - 1 , 1);
    for J = 1:length(dist_grid)-1
        idx = find(DS.unreplicated_cdc20_bin == unq_underrep_bin(I) &...
            DS.dist_to_the_end_kb >= dist_grid(J) & DS.dist_to_the_end_kb < dist_grid(J+1));
        data(J) = nanmean(DS.freq_indel(idx));
    end
    idx = find(~isnan(data));
    plot(dist_grid(idx) , data(idx) , 'color', clrs(4*I-2,:) , 'LineWidth' , 2);
    scatter(dist_grid(idx) , data(idx) , 20 , clrs(4*I-2,:) , 'filled');
end
ylabel('Indel freq');
xlabel('Distance to the end');    

subplot(2,2,4); hold on; grid on; set(gca , 'FontSize' , 12);
for I = 1:length(unq_underrep_bin)
    data = NaN(length(rep_time_grid) - 1 , 1);
    for J = 1:length(rep_time_grid)-1
        idx = find(DS.unreplicated_cdc20_bin == unq_underrep_bin(I) &...
            DS.Trep_spline >= rep_time_grid(J) & DS.Trep_spline < rep_time_grid(J+1));
        data(J) = nanmean(DS.freq_indel(idx));
    end
    idx = find(~isnan(data));
    
    plot(rep_time_grid(idx) , data(idx) , 'color', clrs(4*I-2,:) , 'LineWidth' , 2);
    scatter(rep_time_grid(idx) , data(idx) , 20 , clrs(4*I-2,:) , 'filled');
end
ylabel('Indel freq');
xlabel('Replication timing');    
xlim([0 65]);

print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Manuscript/Fig12_for_SNP_and_INDEL');
%%
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features.mat');
DS.unreplicated_cdc20_bin = NaN(length(DS) , 1);
for I = 1:length(DS)
    if DS.percent_unreplicated_not_trimmed_cdc20(I) <= 10
        DS.unreplicated_cdc20_bin(I) = 1;
    elseif DS.percent_unreplicated_not_trimmed_cdc20(I) <= 25
        DS.unreplicated_cdc20_bin(I) = 2;
    elseif DS.percent_unreplicated_not_trimmed_cdc20(I) <= 40
        DS.unreplicated_cdc20_bin(I) = 3;
    else
        DS.unreplicated_cdc20_bin(I) = 4;
    end
end
%%
idx = find(~isnan(DS.frac_conservation) & ~isnan(DS.percent_unreplicated_not_trimmed_cdc20) );
DS = DS(idx , :);

figure; 
unq_underrep_bin = unique(DS.unreplicated_cdc20_bin);
clrs = parula(4*length(unq_underrep_bin));
dist_grid = [0:3:30];
rep_time_grid = [10:5:70];

subplot(2,2,1); hold on; grid on; set(gca , 'FontSize' , 12);
for I = 1:length(unq_underrep_bin)
    data = NaN(length(dist_grid) - 1 , 1);
    for J = 1:length(dist_grid)-1
        idx = find(DS.unreplicated_cdc20_bin == unq_underrep_bin(I) &...
            DS.dist_to_the_end_kb >= dist_grid(J) & DS.dist_to_the_end_kb < dist_grid(J+1));
        data(J) = nanmean(DS.frac_conservation(idx) < 1);
    end
    idx = find(~isnan(data));
    plot(dist_grid(idx) , data(idx) , 'color', clrs(4*I-2,:) , 'LineWidth' , 2);
    scatter(dist_grid(idx) , data(idx) , 20 , clrs(4*I-2,:) , 'filled');
end
ylabel('SNP freq');
xlabel('Distance to the end');    

subplot(2,2,2); hold on; grid on; set(gca , 'FontSize' , 12);
for I = 1:length(unq_underrep_bin)
    data = NaN(length(rep_time_grid) - 1 , 1);
    for J = 1:length(rep_time_grid)-1
        idx = find(DS.unreplicated_cdc20_bin == unq_underrep_bin(I) &...
            DS.Trep_spline >= rep_time_grid(J) & DS.Trep_spline < rep_time_grid(J+1));
        data(J) = nanmean(DS.frac_conservation(idx) < 1);
    end
    idx = find(~isnan(data));
    
    plot(rep_time_grid(idx) , data(idx) , 'color', clrs(4*I-2,:) , 'LineWidth' , 2);
    scatter(rep_time_grid(idx) , data(idx) , 20 , clrs(4*I-2,:) , 'filled');
end
ylabel('SNP freq');
xlabel('Replication timing');    
xlim([0 65]);



















%%
idx = find(~isnan(DS.frac_conservation) & ~isnan(DS.percent_unreplicated_not_trimmed_cdc20) );
D = DS(idx , :);
figure; hold on; grid on; set(gca , 'FontSize' , 12);
unq_underreplication_bin = unique(DS.unreplicated_cdc20_bin);
for I = 1:length(unq_underreplication_bin)
    idx = find(D.unreplicated_cdc20_bin == unq_underreplication_bin(I));
    scatter(unq_underreplication_bin(I) , nanmean(D.frac_conservation(idx) < 1) , 20 , [.2 .2 .2] , 'filled');
end

%%
clrs = parula(12);
idx = find(~isnan(DS.frac_conservation) & ~isnan(DS.percent_unreplicated_not_trimmed_cdc20) );
D = DS(idx , :);

figure; 
D = sortrows(D , 'percent_unreplicated_not_trimmed_cdc20');
K = 50;
data_x = [];
data_y = [];
for I = 1:length(D)-K
    T = D(I:I+K , :);
    data_x = [data_x nanmedian(T.percent_unreplicated_not_trimmed_cdc20)];
    data_y = [data_y nanmean(T.frac_conservation < 1)];
end
subplot(2,1,1); hold on; grid on; set(gca , 'FontSize' , 12);
plot(data_x , data_y , 'LineWidth' , 3 , 'color' , clrs(3,:));
xlabel('% unreplicted cells, cdc20');
ylabel('% of genes, not fully presented');

D = sortrows(D , 'dist_to_the_end_kb');
data_x = [];
data_y = [];
for I = 1:length(D)-K
    T = D(I:I+K , :);
    data_x = [data_x nanmedian(T.dist_to_the_end_kb)];
    data_y = [data_y nanmean(T.frac_conservation < 1)];
end
subplot(2,1,2); hold on; grid on; set(gca , 'FontSize' , 12);
plot(data_x , data_y , 'LineWidth' , 3 , 'color' , clrs(6,:));
xlabel('distance to the end');
ylabel('% of genes, not fully presented');

%% eCDFs
idx = find(~isnan(DS.frac_conservation) & ~isnan(DS.percent_unreplicated_not_trimmed_cdc20) );
D = DS(idx , :);
%%
K = 100;
idx = find(~isnan(DS.frac_conservation) & ~isnan(DS.percent_unreplicated_not_trimmed_cdc20) );
D = DS(idx , :);
D.DM_percentUnreplicated_over_Distance = NaN(length(D) , 1);
D.DM_fracConservation_over_Distance = NaN(length(D) , 1);
D = sortrows(D , 'dist_to_the_end_kb');
data_x = [];
data_y_fracConservation = [];
data_y_percentUnreplicated = [];
for I = 1:length(D)-K
    T = D(I:I+K , :);
    data_x = [data_x nanmean(T.dist_to_the_end_kb)];
    data_y_fracConservation = [data_y_fracConservation nanmedian(T.frac_conservation)];
    data_y_percentUnreplicated = [data_y_percentUnreplicated nanmedian(T.percent_unreplicated_not_trimmed_cdc20)];
end
for I = 1:length(D)
    D.DM_percentUnreplicated_over_Distance(I) = D.percent_unreplicated_not_trimmed_cdc20(I)-interp1(data_x , data_y_percentUnreplicated , D.dist_to_the_end_kb(I) );
    D.DM_fracConservation_over_Distance(I) = D.frac_conservation(I)-interp1(data_x , data_y_fracConservation , D.dist_to_the_end_kb(I) );
end
D.DM_Distance_over_PercentUnreplication = NaN(length(D) , 1);
D.DM_fracConservation_over_PercentUnreplication = NaN(length(D) , 1);
D = sortrows(D , 'percent_unreplicated_not_trimmed_cdc20');
data_x = [];
data_y_fracConservation = [];
data_y_dist = [];
for I = 1:length(D)-K
    T = D(I:I+K , :);
    data_x = [data_x nanmean(T.percent_unreplicated_not_trimmed_cdc20)];
    data_y_fracConservation = [data_y_fracConservation nanmedian(T.frac_conservation)];
    data_y_dist = [data_y_dist nanmedian(T.dist_to_the_end_kb)];
end
for I = 1:length(D)
    D.DM_Distance_over_PercentUnreplication(I) = D.dist_to_the_end_kb(I)-interp1(data_x , data_y_dist , D.percent_unreplicated_not_trimmed_cdc20(I) );
    D.DM_fracConservation_over_PercentUnreplication(I) = D.frac_conservation(I)-interp1(data_x , data_y_fracConservation , D.percent_unreplicated_not_trimmed_cdc20(I) );
end

%%

clrs = parula(12);
figure; 
subplot(2,2,1);hold on; grid on; set(gca , 'FontSize' , 12);
idx1 = find(~isnan(D.frac_conservation) & ...
    ~isnan(D.percent_unreplicated_not_trimmed_cdc20) & ...
    D.frac_conservation < 1);
[y,x] = ecdf(D.percent_unreplicated_not_trimmed_cdc20(idx1));
plot(x,y , 'LineWidth' , 3 , 'color' , clrs(2,:));
idx2 = find(~isnan(D.frac_conservation) & ...
    ~isnan(D.percent_unreplicated_not_trimmed_cdc20) & ...
    D.frac_conservation == 1);
[y,x] = ecdf(D.percent_unreplicated_not_trimmed_cdc20(idx2));
plot(x,y , 'LineWidth' , 3 , 'color' , clrs(6,:));
[~,p] = ttest2(D.percent_unreplicated_not_trimmed_cdc20(idx1),D.percent_unreplicated_not_trimmed_cdc20(idx2));
title(sprintf('P = %.2f.' , p));
xlim([-50 60]);
xlabel('Underreplication');
ylabel('eCDF');

subplot(2,2,2);hold on; grid on; set(gca , 'FontSize' , 12);
idx1 = find(~isnan(D.frac_conservation) & ...
    ~isnan(D.percent_unreplicated_not_trimmed_cdc20) & ...
    D.frac_conservation < 1);
[y,x] = ecdf(D.dist_to_the_end_kb(idx1));
plot(x,y , 'LineWidth' , 3 , 'color' , clrs(2,:));
idx2 = find(~isnan(D.frac_conservation) & ...
    ~isnan(D.percent_unreplicated_not_trimmed_cdc20) & ...
    D.frac_conservation == 1);
[y,x] = ecdf(D.dist_to_the_end_kb(idx2));
plot(x,y , 'LineWidth' , 3 , 'color' , clrs(6,:));
[~,p] = ttest2(D.dist_to_the_end_kb(idx1),D.dist_to_the_end_kb(idx2));
title(sprintf('P = %.2f.' , p));
xlim([0 600]);
xlabel('Distance');
ylabel('eCDF');

subplot(2,2,3);hold on; grid on; set(gca , 'FontSize' , 12);
idx1 = find(~isnan(D.frac_conservation) & ...
    ~isnan(D.percent_unreplicated_not_trimmed_cdc20) & ...
    D.frac_conservation < 1);
[y,x] = ecdf(D.DM_percentUnreplicated_over_Distance(idx1));
plot(x,y , 'LineWidth' , 3 , 'color' , clrs(2,:));
idx2 = find(~isnan(D.frac_conservation) & ...
    ~isnan(D.percent_unreplicated_not_trimmed_cdc20) & ...
    D.frac_conservation == 1);
[y,x] = ecdf(D.DM_percentUnreplicated_over_Distance(idx2));
plot(x,y , 'LineWidth' , 3 , 'color' , clrs(6,:));
plot(x,y , 'LineWidth' , 3 , 'color' , clrs(6,:));
[~,p] = ttest2(D.DM_percentUnreplicated_over_Distance(idx1),D.DM_percentUnreplicated_over_Distance(idx2));
title(sprintf('P = %.2f.' , p));
xlim([-30 30]);
xlabel('DM underreplication (over distance)');
ylabel('eCDF');

subplot(2,2,4);hold on; grid on; set(gca , 'FontSize' , 12);
idx1 = find(~isnan(D.frac_conservation) & ...
    ~isnan(D.percent_unreplicated_not_trimmed_cdc20) & ...
    D.frac_conservation < 1);
[y,x] = ecdf(D.DM_Distance_over_PercentUnreplication(idx1));
plot(x,y , 'LineWidth' , 3 , 'color' , clrs(2,:));
idx2 = find(~isnan(D.frac_conservation) & ...
    ~isnan(D.percent_unreplicated_not_trimmed_cdc20) & ...
    D.frac_conservation == 1);
[y,x] = ecdf(D.DM_Distance_over_PercentUnreplication(idx2));
plot(x,y , 'LineWidth' , 3 , 'color' , clrs(6,:));
[~,p] = ttest2(D.DM_Distance_over_PercentUnreplication(idx1),D.DM_Distance_over_PercentUnreplication(idx2));
title(sprintf('P = %.2f.' , p));
xlim([-200 500]);
ylabel('eCDF');
xlabel('DM distance (over underreplication)');

set(gcf , 'PaperPosition' , [0 0 20 20]);
print('-dpsc2' , 'eCDF__presented_genes__VS__not__diff_metrics.eps');
%% boxplots

clrs = parula(12);
figure; 

subplot(2,2,1);hold on; grid on; set(gca , 'FontSize' , 12);
idx1 = find(~isnan(D.frac_conservation) & ...
    ~isnan(D.percent_unreplicated_not_trimmed_cdc20) & ...
    D.frac_conservation < 1);
[y,x] = ecdf(D.percent_unreplicated_not_trimmed_cdc20(idx1));
plot(x,y , 'LineWidth' , 3 , 'color' , clrs(2,:));
idx2 = find(~isnan(D.frac_conservation) & ...
    ~isnan(D.percent_unreplicated_not_trimmed_cdc20) & ...
    D.frac_conservation == 1);
[y,x] = ecdf(D.percent_unreplicated_not_trimmed_cdc20(idx2));
plot(x,y , 'LineWidth' , 3 , 'color' , clrs(6,:));
[~,p] = ttest2(D.percent_unreplicated_not_trimmed_cdc20(idx1),D.percent_unreplicated_not_trimmed_cdc20(idx2));
title(sprintf('P = %.2f.' , p));

subplot(2,2,2);hold on; grid on; set(gca , 'FontSize' , 12);
idx1 = find(~isnan(D.frac_conservation) & ...
    ~isnan(D.percent_unreplicated_not_trimmed_cdc20) & ...
    D.frac_conservation < 1);
[y,x] = ecdf(D.dist_to_the_end_kb(idx1));
plot(x,y , 'LineWidth' , 3 , 'color' , clrs(2,:));
idx2 = find(~isnan(D.frac_conservation) & ...
    ~isnan(D.percent_unreplicated_not_trimmed_cdc20) & ...
    D.frac_conservation == 1);
[y,x] = ecdf(D.dist_to_the_end_kb(idx2));
plot(x,y , 'LineWidth' , 3 , 'color' , clrs(6,:));
[~,p] = ttest2(D.dist_to_the_end_kb(idx1),D.dist_to_the_end_kb(idx2));
title(sprintf('P = %.2f.' , p));

subplot(2,2,3);hold on; grid on; set(gca , 'FontSize' , 12);
idx1 = find(~isnan(D.frac_conservation) & ...
    ~isnan(D.percent_unreplicated_not_trimmed_cdc20) & ...
    D.frac_conservation < 1);
[y,x] = ecdf(D.DM_percentUnreplicated_over_Distance(idx1));
%plot(x,y , 'LineWidth' , 3 , 'color' , clrs(2,:));
idx2 = find(~isnan(D.frac_conservation) & ...
    ~isnan(D.percent_unreplicated_not_trimmed_cdc20) & ...
    D.frac_conservation == 1);
data = NaN(length(idx2) , 2);
data(:,1) = D.DM_percentUnreplicated_over_Distance(idx2);
data(1:length(idx1) , 2) = D.DM_percentUnreplicated_over_Distance(idx1);
boxplot(data)

[y,x] = ecdf(D.DM_percentUnreplicated_over_Distance(idx2));
ylim([-10 10]);
%plot(x,y , 'LineWidth' , 3 , 'color' , clrs(6,:));
[~,p] = ttest2(D.DM_percentUnreplicated_over_Distance(idx1),D.DM_percentUnreplicated_over_Distance(idx2));
title(sprintf('P = %.2f.' , p));

subplot(2,2,4);hold on; grid on; set(gca , 'FontSize' , 12);
idx1 = find(~isnan(D.frac_conservation) & ...
    ~isnan(D.percent_unreplicated_not_trimmed_cdc20) & ...
    D.frac_conservation < 1);
idx2 = find(~isnan(D.frac_conservation) & ...
    ~isnan(D.percent_unreplicated_not_trimmed_cdc20) & ...
    D.frac_conservation == 1);
data = NaN(length(idx2) , 2);
data(:,1) = D.DM_Distance_over_PercentUnreplication(idx2);
data(1:length(idx1) , 2) = D.DM_Distance_over_PercentUnreplication(idx1);
boxplot(data)

%[y,x] = ecdf(D.DM_Distance_over_PercentUnreplication(idx2));
%plot(x,y , 'LineWidth' , 3 , 'color' , clrs(6,:));
[~,p] = ttest2(D.DM_Distance_over_PercentUnreplication(idx1),D.DM_Distance_over_PercentUnreplication(idx2));
title(sprintf('P = %.2f.' , p));




%%
K = 100;
idx = find(~isnan(DS.frac_conservation) & ~isnan(DS.percent_unreplicated_not_trimmed_cdc20) );
D = DS(idx , :);
D.DM_percentUnreplicated_over_Distance = NaN(length(D) , 1);
D.DM_fracConservation_over_Distance = NaN(length(D) , 1);
D = sortrows(D , 'dist_to_the_end_kb');

data_x = [];
data_y_fracConservation = [];
data_y_percentUnreplicated = [];
for I = 1:length(D)-K
    T = D(I:I+K , :);
    data_x = [data_x nanmean(T.dist_to_the_end_kb)];
    data_y_fracConservation = [data_y_fracConservation nanmedian(T.frac_conservation)];
    data_y_percentUnreplicated = [data_y_percentUnreplicated nanmedian(T.percent_unreplicated_not_trimmed_cdc20)];
end
for I = 1:length(D)
    D.DM_percentUnreplicated_over_Distance(I) = D.percent_unreplicated_not_trimmed_cdc20(I)-interp1(data_x , data_y_percentUnreplicated , D.dist_to_the_end_kb(I) );
    D.DM_fracConservation_over_Distance(I) = D.frac_conservation(I)-interp1(data_x , data_y_fracConservation , D.dist_to_the_end_kb(I) );
end

%
D.DM_Distance_over_PercentUnreplication = NaN(length(D) , 1);
D.DM_fracConservation_over_PercentUnreplication = NaN(length(D) , 1);
D = sortrows(D , 'percent_unreplicated_not_trimmed_cdc20');
data_x = [];
data_y_fracConservation = [];
data_y_dist = [];
for I = 1:length(D)-K
    T = D(I:I+K , :);
    data_x = [data_x nanmean(T.percent_unreplicated_not_trimmed_cdc20)];
    data_y_fracConservation = [data_y_fracConservation nanmedian(T.frac_conservation)];
    data_y_dist = [data_y_dist nanmedian(T.dist_to_the_end_kb)];
end
for I = 1:length(D)
    D.DM_Distance_over_PercentUnreplication(I) = D.dist_to_the_end_kb(I)-interp1(data_x , data_y_dist , D.percent_unreplicated_not_trimmed_cdc20(I) );
    D.DM_fracConservation_over_PercentUnreplication(I) = D.frac_conservation(I)-interp1(data_x , data_y_fracConservation , D.percent_unreplicated_not_trimmed_cdc20(I) );
end

%%
clrs = parula(12);

figure; 
D = sortrows(D , 'DM_percentUnreplicated_over_Distance');
K = 100;
data_x = [];
data_y = [];
for I = 1:length(D)-K
    T = D(I:I+K , :);
    data_x = [data_x nanmedian(T.DM_percentUnreplicated_over_Distance)];
    data_y = [data_y nanmean(T.frac_conservation < 1)];
end
subplot(2,1,1); hold on; grid on; set(gca , 'FontSize' , 12);
plot(data_x , data_y , 'LineWidth' , 3 , 'color' , clrs(3,:));
xlabel('% unreplicted cells, cdc20');
ylabel('% of genes, not fully presented');

D = sortrows(D , 'DM_Distance_over_PercentUnreplication');
data_x = [];
data_y = [];
for I = 1:length(D)-K
    T = D(I:I+K , :);
    data_x = [data_x nanmedian(T.DM_Distance_over_PercentUnreplication)];
    data_y = [data_y nanmean(T.frac_conservation < 1)];
end
subplot(2,1,2); hold on; grid on; set(gca , 'FontSize' , 12);
plot(data_x , data_y , 'LineWidth' , 3 , 'color' , clrs(6,:));
xlabel('distance to the end');
ylabel('% of genes, not fully presented');

%%







idx = find(~isnan(DS.frac_conservation) & ~isnan(DS.percent_unreplicated_not_trimmed_cdc20) );
D = DS(idx , :);
D = sortrows(D , 'percent_unreplicated_not_trimmed_cdc20');
figure; hold on; grid on; set(gca , 'FontSize' , 12);
dscatter(D.percent_unreplicated_not_trimmed_cdc20 , D.frac_conservation );
K = 50;
data_x = [];
data_y = [];
for I = 1:length(D)-K
    T = D(I:I+K , :);
    data_x = [data_x nanmedian(T.percent_unreplicated_not_trimmed_cdc20)];
    data_y = [data_y nanmedian(T.frac_conservation)];
end
plot(data_x , data_y , 'LineWidth' , 3 , 'color' , 'r');
ylim([.8 1.05]);

%%


boxplot(DS.frac_conservation(idx) , DS.unreplicated_cdc20_bin(idx) , 'color' , [.2 .2 .2]);
unq_underreplication_bin = unique(DS.unreplicated_cdc20_bin);
for I = 1:length(unq_underreplication_bin)
    idx = find(DS.unreplicated_cdc20_bin == unq_underreplication_bin(I));
    scatter(repmat(I , length(idx) , 1) , DS.frac_conservation(idx) , 20 , clrs(4*I-2,:) , 'filled');
end

%%
idx = find(~isnan(DS.frac_conservation) & ~isnan(DS.percent_unreplicated_not_trimmed_cdc20) );
figure; hold on; grid on; set(gca , 'FontSize' , 12);
scatter(DS.percent_unreplicated_not_trimmed_cdc20(idx) , DS.frac_conservation(idx) , 20 , [.2 .2 .2] , 'filled');







%%
idx = find(~isnan(DS.frac_conservation) & ~isnan(DS.percent_unreplicated_not_trimmed_cdc20) &...
    DS.frac_conservation < 1);
D = DS(idx , :);

[p,stat] = anovan(D.frac_conservation , ...
    {D.percent_unreplicated_not_trimmed_cdc20 , ...
     D.dist_to_the_end_kb , D.Trep_spline } , ...
    'continuous' , [1 2 3 ]);

%%







%%
figure; hold on;
idx = find(~isnan(DS.frac_conservation) & ~isnan(DS.percent_unreplicated_not_trimmed_cdc20));
C = corrcoef(DS.frac_conservation(idx) , DS.dist_to_the_end_kb(idx))

%%
DS.subtelomeric_bool = DS.dist_to_the_end_kb <= 82;
DS.arm_bool = zeros(length(DS) , 1);
for I = 1:length(DS)
    if strcmp(DS.arm{I} , 'left')
        DS.arm_bool(I) = -1*DS.chr_num(I);
    elseif strcmp(DS.arm{I} , 'right')
        DS.arm_bool(I) = DS.chr_num(I);
    end
end 

%%
idx = find(~isnan(DS.frac_conservation) & ~isnan(DS.percent_unreplicated_not_trimmed_cdc20) &...
    ~isnan(DS.percent_unreplicated_not_trimmed_dbf2) & DS.arm_bool ~= 0);

[~, stats] = anovan(DS.frac_conservation(idx) , ...
    {DS.arm_bool(idx) , DS.subtelomeric_bool(idx) , ...
    DS.dist_to_the_end_kb(idx) , DS.Trep_spline(idx) , ...
    DS.percent_unreplicated_not_trimmed_cdc20(idx) , ...
    DS.percent_unreplicated_not_trimmed_dbf2(idx) }, ...
    'continuous' , [3:6]);



%%
figure; hold on; 
idx = find(~isnan(DS.frac_conservation) & ~isnan(DS.percent_unreplicated_not_trimmed_cdc20) &...
    ~isnan(DS.percent_unreplicated_not_trimmed_dbf2) & DS.arm_bool ~= 0 & DS.frac_conservation < 1);
[y , x] = ksdensity(DS.percent_unreplicated_not_trimmed_cdc20(idx) , [-10:2:60]);
plot(x , y , 'LineWidth' , 3);

idx = find(~isnan(DS.frac_conservation) & ~isnan(DS.percent_unreplicated_not_trimmed_cdc20) &...
    ~isnan(DS.percent_unreplicated_not_trimmed_dbf2) & DS.arm_bool ~= 0 & DS.frac_conservation < 1);
[y , x] = ksdensity(DS.dist_to_the_end_kb(idx) , [0:10:300]);
plot(x , y , 'LineWidth' , 3);



%%
figure; hold on; 
scatter( DS.percent_unreplicated_not_trimmed_cdc20 , DS.dist_to_the_end_kb)




%%
idx = find(~isnan(DS.frac_conservation) & ~isnan(DS.percent_unreplicated_not_trimmed_cdc20) );
D = DS(idx , :);
figure; hold on; grid on; set(gca , 'FontSize' , 12);
scatter(D.dist_to_the_end_kb , D.frac_conservation)

%%
idx = find(~isnan(DS.frac_conservation) & ~isnan(DS.percent_unreplicated_not_trimmed_cdc20));
D = DS(idx , :);
clrs = parula(12);
figure; hold on; grid on; set(gca , 'FontSize' , 12);

C = corrcoef(D.frac_conservation , D.percent_unreplicated_not_trimmed_cdc20);
bar(1 , C(1,2) , 'FaceColor' , clrs(2,:));

C = corrcoef(D.frac_conservation , D.percent_unreplicated_not_trimmed_dbf2);
bar(2 , C(1,2) , 'FaceColor' , clrs(4,:));

C = corrcoef(D.frac_conservation , D.dist_to_the_end_kb);
bar(3 , C(1,2) , 'FaceColor' , clrs(6,:));

C = corrcoef(D.frac_conservation , D.dist_to_ARS);
bar(4 , C(1,2) , 'FaceColor' , clrs(6,:));

C = corrcoef(D.frac_conservation , D.Trep_spline);
bar(5 , C(1,2) , 'FaceColor' , clrs(6,:));

%%
idx = find(~isnan(DS.frac_conservation) & ~isnan(DS.percent_unreplicated_not_trimmed_cdc20) &...
    DS.frac_conservation < 1);
T = DS(idx , :);
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

parameters_idx = cell(3,1);
parameters_idx{1} = [10]; parameters_idx{2} = [14]; parameters_idx{3} = [10 14]; 
%
R_train = cell(length(parameters_idx),1);
R_test = cell(length(parameters_idx),1);
for I = 1:length(parameters_idx)
    X = double( T(:,parameters_idx{I}));
    Y = double( T(:, 23)) ;
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

figure; hold on; grid on; set(gca , 'FontSize' , 14); hold on; grid on; set(gca , 'FontSize' , 14);
data = [R_test{1} R_train{1} R_test{2} R_train{2} R_test{3} R_train{3} ];
boxplot(data , 'notch' , 'on' , 'color' , [.2 .2 .2]);

clrs = parula(16);
%clrs_set = [clrs2(2,:); clrs2(10,:); clrs1(3,:) ; clrs4(2,:); clrs3(4,:); clrs5(3,:); clrs2(6,:); clrs1(3,:)];
%labels = {'dist' , 'dist+Trep' , 'dist+chr', 'dist+chr+Trep' , 'dist+ARS' , 'dist+expression+Trep+chr' , 'all'};
for I = 1:length(parameters_idx)
    scatter(repmat(2*I-1,K,1)+rand(K,1)*0.2-0.1 , R_test{I} , 70 , clrs(4*I-2,:) , 'filled' , 'd');
    scatter(repmat(2*I,K,1)+rand(K,1)*0.2-0.1 , R_train{I} , 70 , clrs(4*I-2,:) , 'filled' , 'o');
%     [~ , p] = ttest2(R_test{I} , R_train{I});
%     if p < 0.05
%         scatter(3*I-1.5 , 1 , 100 , 'k'  , 'x');
%     end
end
xticklabels = {'1.Test' , '1.Train' , '2.Test' , '2.Train' , '3.Test' , '3.Train' };
set(gca , 'Xtick' , [1:8] , 'XtickLabel' , xticklabels);

ylabel('R^{2}');
xtickangle(-45);
set(gcf , 'PaperPosition' , [0 0 20 20]);
title('Fraction of unreplicated cells, {\it dbf2}.');

%%
idx = find(~isnan(DS.frac_conservation) & ~isnan(DS.percent_unreplicated_not_trimmed_cdc20) );
D = DS(idx , :);

dist_to_the_end_bin = [0 50 100 150 200 ];
figure; hold on; grid on; set(gca , 'FontSize' , 12);
clrs = parula(4*length(dist_to_the_end_bin));
for I = 1:length(dist_to_the_end_bin) - 1
    
    idx = find(D.dist_to_the_end_kb > dist_to_the_end_bin(I) & ...
        D.dist_to_the_end_kb < dist_to_the_end_bin(I+1) &...
        D.percent_unreplicated_not_trimmed_cdc20 > 16);
    bar(3*I-2 , nanmean(D.frac_conservation(idx) < 1) , 'FaceColor' , clrs(4*I-3,:) , ...
        'Display' , 'replicated');
    
    idx = find(D.dist_to_the_end_kb > dist_to_the_end_bin(I) & ...
        D.dist_to_the_end_kb < dist_to_the_end_bin(I+1) &...
        D.percent_unreplicated_not_trimmed_cdc20 < 16);
    bar(3*I-1 , nanmean(D.frac_conservation(idx) < 1) , 'FaceColor' , clrs(4*I-1,:) , ...
        'Display' , 'unreplicated');
    
end
set(gca , 'Xtick' , [1.5 4.5 7.5 10.5] , 'XtickLabel' , {'0-50' , '50-100' , '100-150' , '150-200'});
legend('location' , 'EastOutside'); 
xlabel('Distance to the end');
ylabel('Fraction of not entirely conserved genes');
title('Conditional probability');
print('-dpng' , 'condProb__replicated_or_not__gene_conservation');
%%
dist_to_the_end_bin = [0:5:150];
figure; hold on; grid on; set(gca , 'FontSize' , 12);
clrs = parula(4*length(dist_to_the_end_bin));
clrs1 = summer(12);
clrs2 = hot(12);
for I = 1:length(dist_to_the_end_bin) - 1
    
    idx = find(D.dist_to_the_end_kb > dist_to_the_end_bin(I) & ...
        D.dist_to_the_end_kb < dist_to_the_end_bin(I+1) &...
        D.percent_unreplicated_not_trimmed_cdc20 > 16);
    scatter(I , nanmean(D.frac_conservation(idx) < 1) , 100 , clrs1(3,:) , 'filled');
    
    idx = find(D.dist_to_the_end_kb > dist_to_the_end_bin(I) & ...
        D.dist_to_the_end_kb < dist_to_the_end_bin(I+1) &...
        D.percent_unreplicated_not_trimmed_cdc20 < 16);
     scatter(I , nanmean(D.frac_conservation(idx) < 1) , 100 , clrs2(3,:) , 'filled');
    
end

%%
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
%%
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
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig3/SNP_InDel__underrep' , '-r300');

%%
DS.underrep_DM_bin = NaN(length(DS) , 1);
for I = 1:length(DS)
    if DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman(I) <= 0
        DS.underrep_DM_bin(I) = 1;
    else 
        DS.underrep_DM_bin(I) = 2;
    end
end
%%
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
%xlabel('% unreplicated cells');
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
%xlabel('% unreplicated cells');
ylabel('InDel frequency');
legend('location' , 'SouthOutside');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig3/SNP_InDel__underrepDM' , '-r300');

%%
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features.mat');
idx = find(~isnan(DS.frac_conservation)); DS = DS(idx , :);
DS.underrep_DM_bin = NaN(length(DS) , 1);
for I = 1:length(DS)
    if DS.percent_unreplicated_not_trimmed_cdc20_DM(I) <= 0
        DS.underrep_DM_bin(I) = 1;
    else 
        DS.underrep_DM_bin(I) = 2;
    end
end
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

%%
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
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig3/FracConservation__underrep' , '-r300');

%%
clrs1 = autumn(8);
figure('units','centimeters','position',[5 5 9 9]);
legend_titles = {'<0' , '>0'};
hold on; grid on; set(gca , 'FontSize' , 12);
h1 = boxplot(DS.frac_conservation*100 , DS.underrep_DM_bin , 'color' , [.2 .2 .2] , 'symbol' , '');
ylim([80 100]);
set(h1 , 'LineWidth' , 1.5);
h = findobj(gca,'Tag','Box'); set(h1 , 'LineWidth' , 1.5);
for j=1:2
    patch(get(h(2-j+1),'XData'),get(h(2-j+1),'YData'), clrs1(3*j,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
set(gca , 'Xtick' , []); 
%xlabel('% unreplicated cells');
%ylabel('% strains with preserved gene');
legend('location' , 'SouthOutside');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig3/FracConservation__underrepDM' , '-r300');

