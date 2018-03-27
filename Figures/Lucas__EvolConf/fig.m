%% barplots showing R^2 for model to predict evolution data from distance-to-end or cdc20 arrest
% AM March 2018

%% set paths
cd('~/Develop/Mendoza__ReplicationEvolution/Figures/Lucas__EvolConf');
addpath(genpath('~/Develop/matlab'));

%% for SNP & InDel frequency 
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');

DS.dist_log = log2(DS.dist_to_the_end_kb);
DS = DS(: , {'dist_log' , 'Trep_spline' , 'percent_unreplicated_not_trimmed_cdc20_smooth' , 'freq_SNP' , 'freq_indel'});
idx = find(~isnan(DS.percent_unreplicated_not_trimmed_cdc20_smooth));
T = dataset2table( DS(idx , :) ) ;

% shuffle dataset
T = T( randperm(nrows(T)) ,:);

K = 25;
N = round(size(T,1)/K);
all_idx = [1:size(T,1)];
idx = cell(K , 1);
for I = 1:K-1
    temp_idx = all_idx(1:N);
    all_idx = setdiff(all_idx , temp_idx);
    idx{I} = temp_idx;
end 
idx{K} = all_idx;

%% fit & plot data
idx_output = [4 5];
title_names = {'SNP frequency' , 'InDel frequency'};
save_names = {'SNP__R2' , 'Indel__R2'};
parameters_idx = cell(4,1);
parameters_idx{1} = [1]; parameters_idx{2} = [2]; parameters_idx{3} = [1 2]; parameters_idx{4} = [3]; 
for Z = 1:2
    R_train = cell(length(parameters_idx),1);
    R_test = cell(length(parameters_idx),1);
    for I = 1:length(parameters_idx)
        X = table2array( T(:,parameters_idx{I}));
        Y = table2array( T(:, idx_output(Z))) ;
        temp_train = NaN(K,1);
        temp_test = NaN(K,1);
        for J = 1:K
            
            idx_train = setdiff([1:size(T,1)] , idx{J});
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

    % plot all four bars
    figure; hold on; grid on; set(gca , 'FontSize' , 8); 
    data = [R_test{1} R_test{2} R_test{3} R_test{4} ];
    mean_data = mean(data);
    std_data = std(data);
	clrs = get(gca,'ColorOrder');
    clrs(3,:) = clrs(4,:) ;
    clrs(4,:) = [0 0 0] ;
    for I = 1:numel(mean_data)
        bar(I , mean_data(I) , 'FaceColor' , clrs(I,:));
        errorbar(I , mean_data(I) , std_data(I) , 'LineWidth' , 1, 'color' , [.5 .5 .5]);
    end
    [~,p] = ttest2(data(:,3) , data(:,4));
    ylabel('R^{2}');
    set(gcf , 'PaperPosition' , [0 0 3 6]);
    title( sprintf( '%s 3v4 p=%0.06f' , title_names{Z} , p )); 
    xlim([.5 4.5]);
    set(gca , 'Xtick' , [ ]);
    print('-dpng' , save_names{Z} , '-r600');
    close ;
    
    % plot only two bars
    figure; hold on; grid on; set(gca , 'FontSize' , 8); 
    data = [R_test{3} R_test{4} ];
    mean_data = mean(data);
    std_data = std(data);
	clrs = get(gca,'ColorOrder');
    for I = 1:numel(mean_data)
        bar(I , mean_data(I) , 'FaceColor' , clrs(I,:));
        errorbar(I , mean_data(I) , std_data(I) , 'LineWidth' , 1, 'color' , [.5 .5 .5]);
    end
    [~,p] = ttest2(data(:,1) , data(:,2));
    ylabel('R^{2}');
    set(gcf , 'PaperPosition' , [0 0 2 6]);
    title( sprintf( '%s 3v4 p=%0.06f' , title_names{Z} , p )); 
    xlim([.5 2.5]);
    set(gca , 'Xtick' , [ ]);
    print('-dpng' , [ save_names{Z} '_2' ] , '-r600');
    close ;    
end



