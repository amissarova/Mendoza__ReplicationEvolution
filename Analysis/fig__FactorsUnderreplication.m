%% all figures regarding subtelomeric underreplication
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');

%% dist and rep time
DS.dist_log = log2(DS.dist_to_the_end_kb);
DS.dist_ARS_log = log2(DS.dist_to_ARS);
idx = find(~isnan(DS.percent_unreplicated_not_trimmed_cdc20_smooth) );
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

%
parameters_idx = cell(3,1);
parameters_idx{1} = [23]; parameters_idx{2} = [11]; parameters_idx{3} = [23 11]; 
R_train = cell(length(parameters_idx),1);
R_test = cell(length(parameters_idx),1);
for I = 1:length(parameters_idx)
    X = double( T(:,parameters_idx{I}));
    Y = double( T(:, 18)) ;
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


legend_titles = {'log (distance to the end)' , 'Replication time' , 'both'};
clrs = parula(12); clrs_set = [clrs(7,:) ; clrs(10,:) ; clrs(3,:)];
figure('units','centimeters','position',[5 5 8 8]); hold on; grid on; set(gca , 'FontSize' , 10 ); 
data = [R_test{1} R_test{2} R_test{3}];
h1 = boxplot(data , 'color' , [.2 .2 .2], 'symbol','');
h = findobj(gca,'Tag','Box');set(h1 , 'LineWidth' , 1.9);
for j=1:length(h)
    patch(get(h(3-j+1),'XData'),get(h(3-j+1),'YData'), clrs_set(3-j+1,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
legend('location' , 'best');
title('All gDNA');
set(gca , 'Xtick' , []);
%print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2D_all' , '-r300');

%%
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');
DS.dist_log = log2(DS.dist_to_the_end_kb);
DS.dist_ARS_log = log2(DS.dist_to_ARS);
idx = find(~isnan(DS.percent_unreplicated_not_trimmed_cdc20_smooth) );
T = DS(idx , :);
rand_vec = rand(length(T) , 1);
T.random_num = NaN(length(T) , 1);
for I = 1:length(T)
    T.random_num(I) = rand_vec(I);
end
T = sortrows(T , 'random_num');
K = 25;
N = round(length(T)/K);
all_idx = [1:length(T)];
idx = cell(K , 1);
for I = 1:K-1
    temp_idx = all_idx(1:N);
    all_idx = setdiff(all_idx , temp_idx);
    idx{I} = temp_idx;
end 
idx{K} = all_idx;

%
parameters_idx = cell(3,1);
parameters_idx{1} = [23]; parameters_idx{2} = [23 11]; 
parameters_idx{3} = [23 11 37];  
R_train = cell(length(parameters_idx),1);
R_test = cell(length(parameters_idx),1);
for I = 1:length(parameters_idx)
    X = double( T(:,parameters_idx{I}));
    Y = double( T(:, 18)) ;
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

%%
legend_titles = {'log (distance to the end)' , '+Replication time' , '+Transc. rate' , '+G4' , '+GC'};
clrs = parula(12); clrs1 = hot(12); clrs2 = summer(12); clrs3 = lines(6); 
clrs_set = [clrs(7,:) ; clrs(10,:) ; clrs(3,:) ; clrs1(3,:) ; clrs2(2,:) ; clrs3(4,:)];
figure('units','centimeters','position',[5 5 8 8]); hold on; grid on; set(gca , 'FontSize' , 10 ); 
data = [R_test{1} R_test{2} R_test{3} ];
h1 = boxplot(data , 'color' , [.2 .2 .2], 'symbol','');
h = findobj(gca,'Tag','Box');set(h1 , 'LineWidth' , 1.9);
N = 3;
for j=1:length(h)
    patch(get(h(N-j+1),'XData'),get(h(N-j+1),'YData'), clrs_set(N-j+1,:) ,'FaceAlpha' , .5);
end
legend(legend_titles,'location' , 'SouthOutside');
set(gca , 'Xtick' , []);
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2D_3features' , '-r300');





%%
DS = sortrows(DS, 'dist_to_the_end_kb');
data_x = DS.dist_to_the_end_kb;
data_y = DS.percent_unreplicated_not_trimmed_cdc20_smooth;
S = CalcDM_Newman06( data_x , data_y  , 200 , 1 );
DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman = S.DM_diff(:);
DS.median_values_DM = S.DM_ypred(:);
%% Original dist to the end -- underrep
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');
idx = find(DS.dist_to_the_end_kb <= 250);
DS = DS(idx , :);
clrs2 = hot(12); clrs3 = lines(6); clrs4 = parula(12); clrs5 = summer(12);
figure('units','centimeters','position',[5 5 8 8]); hold on; grid on;
h1 = scatter(DS.dist_to_the_end_kb , DS.percent_unreplicated_not_trimmed_cdc20_smooth , ...
    5 , clrs5(7,:) , 'filled' , 'MarkerFaceAlpha',.1 , 'Display' , 'all 200-bp windows');
h2 = plot(DS.dist_to_the_end_kb , DS.median_values_DM , 'LineWidth' , 3.5 , 'color' , clrs5(1,:) , 'Display' , 'running median');
ylim([-20 65]);
idx = find(DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman < 10 & ...
    DS.percent_unreplicated_not_trimmed_cdc20_smooth > 30 & ...
    DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman > 5, 15);
idx = idx(end);
h3 = plot([DS.dist_to_the_end_kb(idx) DS.dist_to_the_end_kb(idx)] , ...
    [DS.median_values_DM(idx) DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx)] , ...
    'LineWidth' , 1.5 , 'color' , clrs2(3,:), 'Display' , sprintf('DM = %.2f' , DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman(idx)));

idx = find(DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman < -6 & ...
    DS.percent_unreplicated_not_trimmed_cdc20_smooth > -6, 105);
idx = idx(end);
h4 = plot([DS.dist_to_the_end_kb(idx) DS.dist_to_the_end_kb(idx)] , ...
    [DS.median_values_DM(idx) DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx)] , ...
    'LineWidth' , 1.5 , 'color' , clrs3(4,:) ,'Display' , sprintf('DM = %.2f' , DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman(idx)));
%legend('location' , 'ne');

idx = find(DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman < 10 & ...
    DS.percent_unreplicated_not_trimmed_cdc20_smooth > 30 & ...
    DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman > 5, 15);
idx = idx(end);
h5 = scatter(DS.dist_to_the_end_kb(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx),...
    30 , clrs2(3,:) , 'filled');
h6 = plot([DS.dist_to_the_end_kb(idx)-6 DS.dist_to_the_end_kb(idx)+6] , ...
    [DS.median_values_DM(idx) DS.median_values_DM(idx)] , ...
    'LineWidth' , 1.5 , 'color' , clrs2(3,:));
h7 = plot([DS.dist_to_the_end_kb(idx)-6 DS.dist_to_the_end_kb(idx)+6] , ...
    [DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx)-1 DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx)-1] , ...
    'LineWidth' , 1.5 , 'color' , clrs2(3,:));

idx = find(DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman < -6 & ...
    DS.percent_unreplicated_not_trimmed_cdc20_smooth > -6, 105);
idx = idx(end);
h8 = scatter(DS.dist_to_the_end_kb(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx),...
    30 , clrs3(4,:) , 'filled');
h9 = plot([DS.dist_to_the_end_kb(idx)-6 DS.dist_to_the_end_kb(idx)+6] , ...
    [DS.median_values_DM(idx) DS.median_values_DM(idx)] , ...
    'LineWidth' , 1.5 , 'color' , clrs3(4,:));
h10 = plot([DS.dist_to_the_end_kb(idx)-6 DS.dist_to_the_end_kb(idx)+6] , ...
    [DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx)+1 DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx)+1] , ...
    'LineWidth' , 1.5 , 'color' , clrs3(4,:));
%legend([h1 h2 h3 h4] ,'location' , 'ne');
ylim([-9 60]);
set(gca , 'Ytick' , [0 20 40 60]);
xlim([0 250]);
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2E__dist' , '-r300');

%%
data_dist = log10(DS.dist_to_the_end_kb);
clrs2 = hot(12); clrs3 = lines(6); clrs4 = parula(12); clrs5 = summer(12);
figure('units','centimeters','position',[5 5 8 8]); hold on; grid on;
h1 = scatter(data_dist , DS.percent_unreplicated_not_trimmed_cdc20_smooth , ...
    20 , [.7 .7 .7] , 'filled' , 'MarkerFaceAlpha',.2 , 'Display' , 'all 200-bp windows');
h2 = plot(data_dist , DS.median_values_DM , 'LineWidth' , 3.5 , 'color' , clrs2(3,:) , 'Display' , 'running median');
ylim([-20 65]);
idx = find(DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman > 20 , 15);
idx = idx(end);
h3 = plot([data_dist(idx) data_dist(idx)] , ...
    [DS.median_values_DM(idx) DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman(idx)] , ...
    'LineWidth' , 2.5 , 'color' , clrs4(2,:), 'Display' , sprintf('DM = %.2f' , DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman(idx)));

idx = find(DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman < -10 , 1005);
idx = idx(end);
h4 = plot([data_dist(idx) data_dist(idx)] , ...
    [DS.median_values_DM(idx) DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman(idx)] , ...
    'LineWidth' , 2.5 , 'color' , clrs5(3,:) ,'Display' , sprintf('DM = %.2f' , DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman(idx)));
legend('location' , 'ne');

idx = find(DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman > 20 , 15);
idx = idx(end);
h5 = scatter(data_dist(idx) , DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman(idx),...
    30 , clrs4(2,:) , 'filled');
h6 = plot([data_dist(idx)-0.5 data_dist(idx)+0.5] , ...
    [DS.median_values_DM(idx) DS.median_values_DM(idx)] , ...
    'LineWidth' , 1.5 , 'color' , clrs4(2,:));
h7 = plot([data_dist(idx)-0.5 data_dist(idx)+0.5] , ...
    [DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman(idx)-1.6 DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman(idx)-1.6] , ...
    'LineWidth' , 1.5 , 'color' , clrs4(2,:));

idx = find(DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman < -10 , 1005);
idx = idx(end);
h8 = scatter(data_dist(idx) , DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman(idx),...
    30 , clrs5(3,:) , 'filled');
h9 = plot([data_dist(idx)-0.5 data_dist(idx)+0.5] , ...
    [DS.median_values_DM(idx) DS.median_values_DM(idx)] , ...
    'LineWidth' , 1.5 , 'color' , clrs5(3,:));
h10 = plot([data_dist(idx)-0.5 data_dist(idx)+0.5] , ...
    [DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman(idx)+1.6 DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman(idx)+1.6] , ...
    'LineWidth' , 1.5 , 'color' , clrs5(3,:));
legend([h1 h2 h3 h4] ,'location' , 'sw');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2E__log10_dist' , '-r300');


%%
figure ('units','centimeters','position',[5 5 10 10]);
subplot(1,2,2); hold on; grid on; set(gca , 'FontSize' , 10);
trep_grid = [0 40 45 50 70];
data = NaN(length(DS) , length(trep_grid)-1);
clrs = winter(4*length(trep_grid)-4);
for I = 1:length(trep_grid)-1
    idx = find(DS.Trep_spline > trep_grid(I) & DS.Trep_spline < trep_grid(I+1) );
    data(1:length(idx) , I) = DS.percent_underreplicated_cdc20_not_trimmed_DM_dist_Newman(idx);
%     scatter(repmat(I,length(idx),1)+randn(length(idx),1)*0.03-0.015 , ...
%        data(1:length(idx) , I) , 10 , clrs(4*I-2,:) , 'filled' , 'MarkerFaceAlpha' , .1);
end
legend_titles = {'<= 40 min' , '40-45 min' , '45-50 min' , '> 50 min'};
h1 = boxplot(data , 'color' , [.2 .2 .2]  , 'symbol' , '');
h = findobj(gca,'Tag','Box'); set(h1 , 'LineWidth' , 1.5);
for j=1:4
    patch(get(h(4-j+1),'XData'),get(h(4-j+1),'YData'), clrs(4*j-2,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end

ylim([-20 20]);
set(gca , 'Xtick' , []);
%set(gca , 'Xtick' , [1:length(trep_grid)-1] , 'XtickLabel' , ...
%    {'<= 40 min' , '40-45 min' , '45-50 min' , '> 50 min'} );
%xlabel('Replication timing, min after \alpha factor release');
ylabel('under-replication DM');
legend('location' , 'SouthOutside');

subplot(1,2,1); hold on; grid on; set(gca , 'FontSize' , 10);
trep_grid = [0 40 45 50 70];
data = NaN(length(DS) , length(trep_grid)-1);
clrs = winter(4*length(trep_grid)-4);
for I = 1:length(trep_grid)-1
    idx = find(DS.Trep_spline > trep_grid(I) & DS.Trep_spline < trep_grid(I+1) );
    data(1:length(idx) , I) = DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx);
    %scatter(repmat(I,length(idx),1)+randn(length(idx),1)*0.03-0.015 , ...
    %   data(1:length(idx) , I) , 10 , clrs(4*I-2,:) , 'filled' , 'MarkerFaceAlpha' , .1);
end
legend_titles = {'<= 40 min' , '40-45 min' , '45-50 min' , '> 50 min'};
h1 = boxplot(data , 'color' , [.2 .2 .2]  , 'symbol' , '');
h = findobj(gca,'Tag','Box'); set(h1 , 'LineWidth' , 1.5);
for j=1:4
    patch(get(h(4-j+1),'XData'),get(h(4-j+1),'YData'), clrs(4*j-2,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
set(h , 'LineWidth' , 1.5);
ylim([-20 70]);
set(gca , 'Xtick' , []);
%set(gca , 'Xtick' , [1:length(trep_grid)-1] , 'XtickLabel' , ...
%    {'<= 40 min' , '40-45 min' , '45-50 min' , '> 50 min'} );
%xlabel('Replication timing, min after \alpha factor release');
legend('location' , 'SouthOutside');
ylabel('% unreplicated cells');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2F' , '-r300');

%%
clrs = parula(12);
figure('units','centimeters','position',[5 5 8 8]); hold on; grid on;
scatter(DS.Trep_spline , DS.percent_unreplicated_not_trimmed_cdc20_smooth , 20 , clrs(5,:) , 'filled' , 'MarkerFaceAlpha' , .1)
xlim([10 65]);
ylim([-20 70]);
R = corrcoef(DS.Trep_spline , DS.percent_unreplicated_not_trimmed_cdc20_smooth);
title(sprintf('R = %.1f.' , R(1,2)));
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2F__scatter' , '-r300');


%%
D = DS;
unq_dist_to_the_end = [0.2 : 0.2 : 765];
median_unrep_cdc20 = NaN(length(unq_dist_to_the_end) , 1);
median_unrep_dbf2 = NaN(length(unq_dist_to_the_end) , 1);
for I = 1:length(unq_dist_to_the_end)
    idx = find(D.dist_to_the_end_kb >= unq_dist_to_the_end(I)-0.1 & D.dist_to_the_end_kb <= unq_dist_to_the_end(I)+0.1);
    median_unrep_cdc20(I) = nanmedian(D.percent_unreplicated_not_trimmed_cdc20_smooth(idx));
    median_unrep_dbf2(I) = nanmedian(D .percent_unreplicated_not_trimmed_dbf2_smooth(idx));
end
DS.percent_underreplicated_cdc20_not_trimmed_DM_dist = NaN(length(DS) , 1);
DS.percent_underreplicated_dbf2_not_trimmed_DM_dist = NaN(length(DS) , 1);
for I = 1:length(DS)
    vq = interp1(unq_dist_to_the_end , median_unrep_cdc20 , DS.dist_to_the_end_kb(I) , 'spline');
    DS.percent_underreplicated_cdc20_not_trimmed_DM_dist(I) = DS.percent_unreplicated_not_trimmed_cdc20(I) - vq;
    vq = interp1(unq_dist_to_the_end, median_unrep_dbf2 , DS.dist_to_the_end_kb(I) ,  'spline');
    DS.percent_underreplicated_dbf2_not_trimmed_DM_dist(I) = DS.percent_unreplicated_not_trimmed_dbf2(I) - vq;
end
save('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat','DS');
%%
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features.mat' , 'DS');
D = DS;
idx = find(strcmp(D.TYPE , 'ARS'));
D = D(idx , :);
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat' , 'DS');
DS.dist_to_ARS = NaN(length(DS) , 1);
for I = 1:length(DS)
    idx = find(D.chr_num == DS.chr_num(I));
    data = D.middle_point_kb(idx);
    DS.dist_to_ARS(I) = nanmin(abs(DS.middle_point_kb(I) - data));
end
save('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat','DS');

%%
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');

figure; 
idx = find(~isnan(DS.percent_unreplicated_not_trimmed_cdc20_smooth) & ...
    DS.percent_unreplicated_not_trimmed_cdc20_smooth < 10000 & ...
	DS.percent_unreplicated_not_trimmed_cdc20_smooth > -10000 & ...
    ~isnan(DS.percent_underreplicated_cdc20_not_trimmed_DM_dist) & ...
    DS.percent_underreplicated_cdc20_not_trimmed_DM_dist < 10000 & ...
	DS.percent_underreplicated_cdc20_not_trimmed_DM_dist > -10000  );
%DS = DS(idx , :);
%
subplot(3,2,1); hold on; grid on; set(gca , 'FontSize' , 12);
dscatter(DS.dist_to_the_end_kb(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx) );
[C,p] = corrcoef(DS.dist_to_the_end_kb(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx));
title(sprintf('R = %.2f, P = %.2f.' , C(1,2) , p(1,2)));
xlabel('Distance to the end');
ylim([-50 100]);
%ylabel('Degree Underreplication');

subplot(3,2,2); hold on; grid on; set(gca , 'FontSize' , 12);
dscatter(DS.dist_to_the_end_kb(idx) , DS.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx) );
[C,p] = corrcoef(DS.dist_to_the_end_kb(idx) , DS.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx));
title(sprintf('R = %.2f, P = %.2f.' , C(1,2) , p(1,2)));
xlabel('Distance to the end');
ylim([-50 50]);
%ylabel('Degree Underreplication DM');

subplot(3,2,3); hold on; grid on; set(gca , 'FontSize' , 12);
dscatter(DS.Trep_spline(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx) );
[C,p] = corrcoef(DS.Trep_spline(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx));
title(sprintf('R = %.2f, P = %.2f.' , C(1,2) , p(1,2)));
xlabel('Replication timing, minutes after release from alpha-factor arrest');
ylim([-50 100]);
%ylabel('Degree Underreplication');
xlim([0 65]);

subplot(3,2,4); hold on; grid on; set(gca , 'FontSize' , 12);
dscatter(DS.Trep_spline(idx) , DS.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx) );
[C,p] = corrcoef(DS.Trep_spline(idx) , DS.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx));
title(sprintf('R = %.2f, P = %.2f.' , C(1,2) , p(1,2)));
xlabel('Replication timing');
%ylabel('Degree Underreplication DM');
xlim([0 65]);
ylim([-100 50]);

subplot(3,2,5); hold on; grid on; set(gca , 'FontSize' , 12);
dscatter(DS.dist_to_ARS(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx) );
[C,p] = corrcoef(DS.dist_to_ARS(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx));
title(sprintf('R = %.2f, P = %.2f.' , C(1,2) , p(1,2)));
ylim([-50 100]);
xlabel('Distance to ARS');
%ylabel('Degree Underreplication');
xlim([0 45]);

subplot(3,2,6); hold on; grid on; set(gca , 'FontSize' , 12);
dscatter(DS.dist_to_ARS(idx) , DS.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx));
[C,p] = corrcoef(DS.dist_to_ARS(idx) , DS.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx));
title(sprintf('R = %.2f, P = %.2f.' , C(1,2) , p(1,2)));
xlabel('Distance to ARS');
%ylabel('Degree Underreplication DM');
xlim([0 45]);
ylim([-100 50]);
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Manuscript/Fig6A');

%%
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');

figure; 
idx = find(~isnan(DS.percent_unreplicated_not_trimmed_cdc20_smooth) & ...
    DS.percent_unreplicated_not_trimmed_cdc20_smooth < 10000 & ...
	DS.percent_unreplicated_not_trimmed_cdc20_smooth > -10000 & ...
    ~isnan(DS.percent_underreplicated_cdc20_not_trimmed_DM_dist) & ...
    DS.percent_underreplicated_cdc20_not_trimmed_DM_dist < 10000 & ...
	DS.percent_underreplicated_cdc20_not_trimmed_DM_dist > -10000  );
%DS = DS(idx , :);
%
subplot(3,1,1); hold on; grid on; set(gca , 'FontSize' , 12);
dscatter(DS.dist_to_the_end_kb(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx) );
[C,p] = corrcoef(DS.dist_to_the_end_kb(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx));
title(sprintf('R = %.2f, P = %.2f.' , C(1,2) , p(1,2)));
xlabel('KB from the telomere');
ylim([-50 100]);
colorbar;
%ylabel('Degree Underreplication');

subplot(3,1,2); hold on; grid on; set(gca , 'FontSize' , 12);
dscatter(DS.Trep_spline(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx) );
[C,p] = corrcoef(DS.Trep_spline(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx));
title(sprintf('R = %.2f, P = %.2f.' , C(1,2) , p(1,2)));
xlabel('Replication timing, minutes after release from alpha-factor arrest');
ylim([-50 100]);
%ylabel('Degree Underreplication');
xlim([0 65]);
colorbar;

subplot(3,1,3); hold on; grid on; set(gca , 'FontSize' , 12);
dscatter(DS.dist_to_ARS(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx) );
[C,p] = corrcoef(DS.dist_to_ARS(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx));
title(sprintf('R = %.2f, P = %.2f.' , C(1,2) , p(1,2)));
ylim([-50 100]);
xlabel('Distance to the closest ARS, KB');
%ylabel('Degree Underreplication');
xlim([0 45]);
colorbar;

print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Manuscript/Fig6A');


%% RT
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');
idx = find(~isnan(DS.percent_unreplicated_not_trimmed_cdc20_smooth) & ...
    DS.percent_unreplicated_not_trimmed_cdc20_smooth < 10000 & ...
	DS.percent_unreplicated_not_trimmed_cdc20_smooth > -10000 & ...
    ~isnan(DS.percent_underreplicated_cdc20_not_trimmed_DM_dist) & ...
    DS.percent_underreplicated_cdc20_not_trimmed_DM_dist < 10000 & ...
	DS.percent_underreplicated_cdc20_not_trimmed_DM_dist > -10000  );
DS = DS(idx , :);
figure; hold on; grid on; set(gca , 'FontSize' , 12);
data = DS.Trep_spline;
% trep_grid = [0 quantile(data,0.9) ...
%     quantile(data,0.95) quantile(data,0.99) quantile(data,0.995) quantile(data,1)];
trep_grid = [0 40 45 50 70];
data = NaN(length(DS) , length(trep_grid)-1);
clrs = winter(4*length(trep_grid)-4);
for I = 1:length(trep_grid)-1
    idx = find(DS.Trep_spline > trep_grid(I) & DS.Trep_spline < trep_grid(I+1) );
    data(1:length(idx) , I) = DS.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx);
    scatter(repmat(I,length(idx),1)+randn(length(idx),1)*0.03-0.015 , ...
       data(1:length(idx) , I) , 10 , clrs(4*I-2,:) , 'filled');
end
h = boxplot(data , 'color' , [.2 .2 .2] , 'notch' , 'on' , 'symbol' , '');
set(h , 'LineWidth' , 1.5);
ylim([-35 40]);
set(gca , 'Xtick' , []);
set(gca , 'Xtick' , [1:length(trep_grid)-1] , 'XtickLabel' , ...
    {'<= 40 min' , '40-45 min' , '45-50 min' , '> 50 min'} );
xlabel('Replication timing, min after \alpha factor release');
ylabel('underreplication DM');
set(gcf, 'PaperPosition' , [0 0 12 12]);
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Manuscript/Fig6B');


%% legend for 6B
figure; hold on; grid on; set(gca , 'FontSize' , 20);
clrs = winter(4*length(trep_grid)-4);
for I = 1:length(trep_grid)-1
    idx = find(DS.Trep_spline > trep_grid(I) & DS.Trep_spline < trep_grid(I+1) );
    data(1:length(idx) , I) = DS.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx);
    scatter(repmat(I,length(idx),1)+randn(length(idx),1)*0.03-0.015 , ...
       data(1:length(idx) , I) , 10 , clrs(4*I-2,:) , 'filled');
end
legend('<= 0.9' , '0.9-0.95' , '0.95-0.99' , '0.99-0.995' , '> 0.995',...
    'location' , 'best');

%% dist to ARS
idx = find(~isnan(DS.percent_unreplicated_not_trimmed_cdc20_smooth) & ...
    DS.percent_unreplicated_not_trimmed_cdc20_smooth < 10000 & ...
	DS.percent_unreplicated_not_trimmed_cdc20_smooth > -10000 & ...
    ~isnan(DS.percent_underreplicated_cdc20_not_trimmed_DM_dist) & ...
    DS.percent_underreplicated_cdc20_not_trimmed_DM_dist < 10000 & ...
	DS.percent_underreplicated_cdc20_not_trimmed_DM_dist > -10000  );
DS = DS(idx , :);
figure; hold on; grid on; set(gca , 'FontSize' , 16);
data = DS.dist_to_ARS;
ARS_grid = [0 ...
    quantile(data,0.999) quantile(data,1)];
data = NaN(length(DS) , length(ARS_grid)-1);
clrs = spring(4*length(ARS_grid)-4);
for I = 1:length(ARS_grid)-1
    idx = find(DS.dist_to_ARS > ARS_grid(I) & DS.dist_to_ARS < ARS_grid(I+1) );
    data(1:length(idx) , I) = DS.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx);
    scatter(repmat(I,length(idx),1)+randn(length(idx),1)*0.03-0.015 , ...
       data(1:length(idx) , I) , 10 , clrs(4*I-2,:) , 'filled');
end
boxplot(data , 'color' , [.2 .2 .2] , 'notch' , 'on' , 'symbol' , '');
ylim([-30 30]);
set(gca , 'Xtick' , [1:length(ARS_grid)-1] , 'XtickLabel' , ...
    {'<99.9%' , '>99.9%'});
xlabel('Quantile of distance to ARS');
ylabel('Degree underreplication DM');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Manuscript/Fig6C');



%%
idx = find(~isnan(DS.percent_unreplicated_not_trimmed_cdc20_smooth) & ...
    DS.percent_unreplicated_not_trimmed_cdc20_smooth < 10000 & ...
	DS.percent_unreplicated_not_trimmed_cdc20_smooth > -10000 & ...
    ~isnan(DS.percent_underreplicated_cdc20_not_trimmed_DM_dist) & ...
    DS.percent_underreplicated_cdc20_not_trimmed_DM_dist < 10000 & ...
	DS.percent_underreplicated_cdc20_not_trimmed_DM_dist > -10000  );
D = DS(idx , :);
figure; 
clrs1 = parula(12);
clrs2 = bone(12);

subplot(2,1,1);hold on; grid on; set(gca , 'FontSize' , 12);
idx = find(D.Trep_spline >= quantile(D.Trep_spline , 0.99));
[y,x] = ksdensity(D.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx) , [-30:2:40]);
plot(x,y/sum(y) , 'LineWidth' , 3 , 'Display' , 'Top 1% by replication time' , 'color' , clrs1(1,:));
idx = find(D.Trep_spline < quantile(D.Trep_spline , 0.99));
[y,x] =  ksdensity(D.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx) , [-30:2:40] );
plot(x,y/sum(y) , 'LineWidth' , 3 , 'Display' , 'Bottom 99% by replication time','color' , clrs1(4,:));
xlabel('Degree underreplication DM');
legend('location' , 'best');
ylabel('PDF');
[~,p] = ttest2(D.percent_underreplicated_cdc20_not_trimmed_DM_dist(D.Trep_spline >= quantile(D.Trep_spline , 0.99)),...
D.percent_underreplicated_cdc20_not_trimmed_DM_dist(D.Trep_spline < quantile(D.Trep_spline , 0.99)) );
title(sprintf('P = %.2f.' , p));


subplot(2,1,2);hold on; grid on; set(gca , 'FontSize' , 12);
idx = find(D.Trep_spline >= quantile(D.dist_to_ARS , 0.99));
[y,x] = ksdensity(D.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx) , [-30:2:40]);
plot(x,y/sum(y) , 'LineWidth' , 3 , 'Display' , 'Top 1% distance to ARS','color' , clrs2(2,:));
idx = find(D.Trep_spline < quantile(D.dist_to_ARS , 0.99));
[y,x] =  ksdensity(D.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx) , [-30:2:40] );
plot(x,y/sum(y) , 'LineWidth' , 3 , 'Display' , 'Bottom 99% distance to ARS','color' , clrs2(7,:));
xlabel('Degree underreplication DM');
legend('location' , 'best');
ylabel('PDF');
[~,p] = ttest2(D.percent_underreplicated_cdc20_not_trimmed_DM_dist(D.dist_to_ARS >= quantile(D.dist_to_ARS , 0.99)),...
D.percent_underreplicated_cdc20_not_trimmed_DM_dist(D.dist_to_ARS < quantile(D.dist_to_ARS , 0.99)) );
title(sprintf('P = %.2f.' , p));








