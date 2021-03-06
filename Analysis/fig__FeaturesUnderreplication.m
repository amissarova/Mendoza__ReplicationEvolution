%% Figures 6

%%
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features.mat');

%
D = dataset('file' , 'mean_PROseq_minus.bed'); D = unique(D);
DS = join(DS , D , 'Type', 'left' , 'Keys' , {'chr' , 'start_point' , 'end_point'} , 'MergeKeys' , true);
D = dataset('file' , 'mean_PROseq_plus.bed'); D = unique(D);
DS = join(DS , D , 'Type', 'left' , 'Keys' , {'chr' , 'start_point' , 'end_point'} , 'MergeKeys' , true);
DS.mean_PROseq = NaN(length(DS) , 1);
for I = 1:length(DS)
    if strcmp(DS.mean_minus_PROseq{I} , '.') & ~strcmp(DS.mean_plus_PROseq{I} , '.')
        DS.mean_PROseq(I) = str2num(DS.mean_plus_PROseq{I});
    elseif ~strcmp(DS.mean_minus_PROseq{I} , '.') & strcmp(DS.mean_plus_PROseq{I} , '.')
        DS.mean_PROseq(I) = abs(str2num(DS.mean_minus_PROseq{I}));
    elseif ~strcmp(DS.mean_minus_PROseq{I} , '.') & ~strcmp(DS.mean_plus_PROseq{I} , '.')
        DS.mean_PROseq(I) = nanmean([abs(str2num(DS.mean_minus_PROseq{I})) str2num(DS.mean_plus_PROseq{I})]);
    end
end
        
%%
idx = find(~isnan(DS.median_PROseq) & strcmp(DS.TYPE , 'ORF') );
DS = DS(idx , :);
DS.PROseq_bin = NaN(length(DS) , 1);
for I = 1:length(DS)
    if DS.max_PROseq(I) < 76.2
        DS.PROseq_bin(I) = 1;
    else
        DS.PROseq_bin(I) = 2;
    end
end

%% Fig 2H: transcription rate VS under-rep and under-rep DM
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
%print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2H' , '-r300');
%%
figure; hold on;
scatter(DS.percent_unreplicated_not_trimmed_cdc20_smooth , DS.median_PROseq);

%%
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features.mat');
%%
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
% for I = 2:length(type_of_interest)
%     temp_data_ORF = data(:,1);
%     temp_data = data(:,I);
%     [~,p] = ttest2(temp_data_ORF(~isnan(temp_data_ORF)) , temp_data(~isnan(temp_data)));
%     if p < .005 & nanmean(temp_data) > nanmean(temp_data_ORF)
%         scatter(I-0.1,45, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
%         scatter(I+0.1,45, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
%     elseif p < .05 & nanmean(temp_data) > nanmean(temp_data_ORF)
%         scatter(I,45, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
%     end
% end
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
% for I = 2:length(type_of_interest)
%     temp_data_ORF = data_DM(:,1);
%     temp_data = data_DM(:,I);
%     [~,p] = ttest2(temp_data_ORF(~isnan(temp_data_ORF)) , temp_data(~isnan(temp_data)));
%     if p < .005 & nanmean(temp_data) > nanmean(temp_data_ORF)
%         scatter(I-0.1,45, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
%         scatter(I+0.1,45, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
%     elseif p < .05 & nanmean(temp_data) > nanmean(temp_data_ORF)
%         scatter(I,45, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
%     end
% end
ylim([-20 40]);
set(gca , 'Xtick' , '');
ylabel('under-replication DM');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2I' , '-r300');
%%
idx = find(strcmp(DS.TYPE , 'ORF') & DS.chr_num == 1 & (DS.mRNA_minus < -1 | DS.mRNA_plus > 1));
D = DS(idx , :);
D.middle_point = (D.start_point + D.end_point)/2;

%%
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features.mat');
idx = find(strcmp(DS.TYPE , 'ORF') & ~isnan(DS.median_PROseq));
DS = DS(idx , :);
%
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
%%
parameters_idx = cell(3,1);
parameters_idx{1} = [11]; parameters_idx{2} = [11 13]; parameters_idx{3} = [11 13 41]; 
R_train = cell(length(parameters_idx),1);
R_test = cell(length(parameters_idx),1);
for I = 1:length(parameters_idx)
    X = double( T(:,parameters_idx{I}));
    Y = double( T(:, 24)) ;
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
legend_titles = {'log (distance to the end)' , '+Replication time' , '+TR'};
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















%%
type_of_interest = {'ORF' , 'ARS' ,'centromere' , 'tRNA' , 'transposable_element_gene' , 'LOH' , ...
    'interstitial_del_dup' , 'terminal_del_dup'};
data = NaN(length(DS) , length(type_of_interest));
data_DM = NaN(length(DS) , length(type_of_interest));
data_dbf = NaN(length(DS) , length(type_of_interest));
data_dbf_DM = NaN(length(DS) , length(type_of_interest));
for I = 1:length(type_of_interest)
    idx = find(strcmp(DS.TYPE , type_of_interest{I}));
    data(1:length(idx) , I) = DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx);
    data_DM(1:length(idx) , I) = DS.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx);
    data_dbf(1:length(idx) , I) = DS.percent_unreplicated_not_trimmed_dbf2_smooth(idx);
    data_dbf_DM(1:length(idx) , I) = DS.percent_underreplicated_dbf2_not_trimmed_DM_dist(idx);
end
clrs1 = summer(12); clrs2 = hot(12); clrs3 = parula(12); clrs4 = lines(6); clrs5 = spring(12);
clrs_set = [.65 .65 .65 ; clrs1(3,:) ; clrs2(3,:) ; clrs3(5,:) ; ...
    clrs4(4,:) ; clrs3(10,:) ; clrs3(3,:) ; clrs5(3,:)];

figure; 
subplot(2,2,1); hold on; grid on; set(gca , 'FontSize' , 12);
for I = 1:length(type_of_interest)
    temp_data = data(:,I);
    temp_data = temp_data(~isnan(temp_data));

end
for I = 1:length(type_of_interest)
    distributionPlot(data(:,I) , 'color' , clrs_set(I,:) , 'xvalues' , I);
end
set(h , 'LineWidth' , 1.5);
for I = 2:length(type_of_interest)
    temp_data_ORF = data(:,1);
    temp_data = data(:,I);
    [~,p] = ttest2(temp_data_ORF(~isnan(temp_data_ORF)) , temp_data(~isnan(temp_data)));
    if p < .005
        scatter(I-0.1,65, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
        scatter(I+0.1,65, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
    elseif p < .05
        scatter(I,65, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
    end
end
ylim([-20 70]);
set(gca , 'Xtick' , '');
%legend('ORF' , 'ARS','centromeres' , 'tRNA' , 'transposons' , ...
%    'LOH-loci' , 'Interstitial dup/del' , 'Terminal dup/del' , 'location' , 'EastOutside');
title('Metaphase');
ylabel('Degree of underreplication');

%%

subplot(2,2,2); hold on; grid on; set(gca , 'FontSize' , 12);
for I = 1:length(type_of_interest)
    temp_data = data_dbf(:,I);
    temp_data = temp_data(~isnan(temp_data));
    scatter(repmat(I,length(temp_data) , 1)+randn(length(temp_data) , 1)*0.03-0.015 , temp_data ,...
        20 , clrs_set(I,:) , 'filled');
end
h = boxplot(data_dbf , 'color' , [.2 .2 .2] , 'symbol' , '');
set(h , 'LineWidth' , 1.5);
for I = 2:length(type_of_interest)
    temp_data_ORF = data_dbf(:,1);
    temp_data = data_dbf(:,I);
    [~,p] = ttest2(temp_data_ORF(~isnan(temp_data_ORF)) , temp_data(~isnan(temp_data)));
    if p < .005
        scatter(I-0.1,55, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
        scatter(I+0.1,55, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
    elseif p < .05
        scatter(I,55, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
    end
end
ylim([-40 60]);
set(gca , 'Xtick' , '');
%legend('ORF' , 'ARS','centromeres' , 'tRNA' , 'transposons' , ...
%    'LOH-loci' , 'Interstitial dup/del' , 'Terminal dup/del' , 'location' , 'EastOutside');
title('Late anaphase');
ylabel('Degree of underreplication');


subplot(2,2,3); hold on; grid on; set(gca , 'FontSize' , 12);
for I = 1:length(type_of_interest)
    temp_data = data_DM(:,I);
    temp_data = temp_data(~isnan(temp_data));
    scatter(repmat(I,length(temp_data) , 1)+randn(length(temp_data) , 1)*0.03-0.015 , temp_data ,...
        20 , clrs_set(I,:) , 'filled');
end
h = boxplot(data_DM , 'color' , [.2 .2 .2] , 'symbol' , '');
set(h , 'LineWidth' , 1.5);
for I = 2:length(type_of_interest)
    temp_data_ORF = data_DM(:,1);
    temp_data = data_DM(:,I);
    [~,p] = ttest2(temp_data_ORF(~isnan(temp_data_ORF)) , temp_data(~isnan(temp_data)));
    if p < .005
        scatter(I-0.1,35, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
        scatter(I+0.1,35, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
    elseif p < .05
        scatter(I,35, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
    end
end
ylim([-20 40]);
%legend('ORF' , 'ARS','centromeres' , 'tRNA' , 'transposons' , ...
%    'LOH-loci' , 'Interstitial dup/del' , 'Terminal dup/del' , 'location' , 'EastOutside');
ylabel('Degree of underreplication DM');
title('Metaphase');
set(gca , 'Xtick' , '');

subplot(2,2,4); hold on; grid on; set(gca , 'FontSize' , 12);
for I = 1:length(type_of_interest)
    temp_data = data_dbf_DM(:,I);
    temp_data = temp_data(~isnan(temp_data));
    scatter(repmat(I,length(temp_data) , 1)+randn(length(temp_data) , 1)*0.03-0.015 , temp_data ,...
        20 , clrs_set(I,:) , 'filled');
end
h = boxplot(data_dbf_DM , 'color' , [.2 .2 .2] , 'symbol' , '');
set(h , 'LineWidth' , 1.5);
for I = 2:length(type_of_interest)
    temp_data_ORF = data_dbf_DM(:,1);
    temp_data = data_dbf_DM(:,I);
    [~,p] = ttest2(temp_data_ORF(~isnan(temp_data_ORF)) , temp_data(~isnan(temp_data)));
    if p < .005
        scatter(I-0.1,35, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
        scatter(I+0.1,35, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
    elseif p < .05
        scatter(I,35, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
    end
end
ylim([-30 40]);
set(gca , 'Xtick' , '');
%legend('ORF' , 'ARS','centromeres' , 'tRNA' , 'transposons' , ...
%    'LOH-loci' , 'Interstitial dup/del' , 'Terminal dup/del' , 'location' , 'EastOutside');
title('Late anaphase');
ylabel('Degree of underreplication DM');

set(gcf , 'PaperPosition' , [0 0 20 20]);
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Manuscript/Fig5A__opt2');

%%
type_of_interest = {'ORF' , 'ARS' ,'centromere' , 'tRNA' , 'transposable_element_gene' , 'LOH' , ...
    'interstitial_del_dup' , 'terminal_del_dup'};
data = NaN(length(DS) , length(type_of_interest));
data_DM = NaN(length(DS) , length(type_of_interest));
data_dbf = NaN(length(DS) , length(type_of_interest));
data_dbf_DM = NaN(length(DS) , length(type_of_interest));
for I = 1:length(type_of_interest)
    idx = find(strcmp(DS.TYPE , type_of_interest{I}));
    data(1:length(idx) , I) = DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx);
    data_DM(1:length(idx) , I) = DS.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx);
    data_dbf(1:length(idx) , I) = DS.percent_unreplicated_not_trimmed_dbf2_smooth(idx);
    data_dbf_DM(1:length(idx) , I) = DS.percent_underreplicated_dbf2_not_trimmed_DM_dist(idx);
end

for I = 6:8
    figure; 
    subplot(2,2,1); hold on; grid on; set(gca , 'FontSize' , 12);
    temp_data_ORF = data(:,1); temp_data_ORF = temp_data_ORF(~isnan(temp_data_ORF));
    temp_data = data(:,I); temp_data = temp_data(~isnan(temp_data));
    N_ORF = length(temp_data_ORF);
    N = length(temp_data);
    bulk_means = NaN(1000 , 1);
    for J = 1:10000
        idx = randsample([1:N_ORF] , N);
        bulk_means(J) = mean(temp_data_ORF(idx));
    end
    [y,x] = hist(bulk_means , [-5:0.5:10]);    
    plot(x , y/sum(y) , 'LineWidth' , 2 , 'color' , [.2 .2 .2]);
    plot([mean(temp_data) mean(temp_data)] , [0 0.5] , 'LineWidth' , 3) ;
    title(strcat( type_of_interest{I} , ':' , sprintf('%.5f' , nanmean( mean(temp_data) > bulk_means) )));
    xlabel('meta, degree');
    
    subplot(2,2,2); hold on; grid on; set(gca , 'FontSize' , 12);
    temp_data_ORF = data_dbf(:,1); temp_data_ORF = temp_data_ORF(~isnan(temp_data_ORF));
    temp_data = data_dbf(:,I); temp_data = temp_data(~isnan(temp_data));
    N_ORF = length(temp_data_ORF);
    N = length(temp_data);
    bulk_means = NaN(1000 , 1);
    for J = 1:10000
        idx = randsample([1:N_ORF] , N);
        bulk_means(J) = mean(temp_data_ORF(idx));
    end
    [y,x] = hist(bulk_means , [-5:0.5:10]);    
    plot(x , y/sum(y) , 'LineWidth' , 2 , 'color' , [.2 .2 .2]);
    plot([mean(temp_data) mean(temp_data)] , [0 0.5] , 'LineWidth' , 3) ;
    title(strcat( type_of_interest{I} , ':' , sprintf('%.5f' , nanmean( mean(temp_data) > bulk_means) )));
    xlabel('ana, degree');
    
    subplot(2,2,3); hold on; grid on; set(gca , 'FontSize' , 12);
    temp_data_ORF = data_DM(:,1); temp_data_ORF = temp_data_ORF(~isnan(temp_data_ORF));
    temp_data = data_DM(:,I); temp_data = temp_data(~isnan(temp_data));
    N_ORF = length(temp_data_ORF);
    N = length(temp_data);
    bulk_means = NaN(1000 , 1);
    for J = 1:10000
        idx = randsample([1:N_ORF] , N);
        bulk_means(J) = mean(temp_data_ORF(idx));
    end
    [y,x] = hist(bulk_means , [-5:0.5:10]);    
    plot(x , y/sum(y) , 'LineWidth' , 2 , 'color' , [.2 .2 .2]);
    plot([mean(temp_data) mean(temp_data)] , [0 0.5] , 'LineWidth' , 3) ;
    title(strcat( type_of_interest{I} , ':' , sprintf('%.5f' , nanmean( mean(temp_data) > bulk_means) )));
    xlabel('meta, degree DM');

    subplot(2,2,4); hold on; grid on; set(gca , 'FontSize' , 12);
    temp_data_ORF = data_dbf_DM(:,1); temp_data_ORF = temp_data_ORF(~isnan(temp_data_ORF));
    temp_data = data_dbf_DM(:,I); temp_data = temp_data(~isnan(temp_data));
    N_ORF = length(temp_data_ORF);
    N = length(temp_data);
    bulk_means = NaN(1000 , 1);
    for J = 1:10000
        idx = randsample([1:N_ORF] , N);
        bulk_means(J) = mean(temp_data_ORF(idx));
    end
    [y,x] = hist(bulk_means , [-5:0.5:10]);    
    plot(x , y/sum(y) , 'LineWidth' , 2 , 'color' , [.2 .2 .2]);
    plot([mean(temp_data) mean(temp_data)] , [0 0.5] , 'LineWidth' , 3) ;
    title(strcat( type_of_interest{I} , ':' , sprintf('%.5f' , nanmean( mean(temp_data) > bulk_means) )));
    xlabel('ana, degree DM');
    
    
    save_name = strcat(type_of_interest{I}  , '--distr');
    print('-dpng' , save_name);
end



%% get which ORFs are on top/bottom by both metrics
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features.mat');
idx = find(~isnan(DS.percent_unreplicated_not_trimmed_cdc20) & strcmp(DS.TYPE , 'ORF'));
DS = DS(idx , :);

DS = sortrows(DS , 'percent_unreplicated_not_trimmed_cdc20' , 'ascend');
D = DS(: , {'ORF'});
export(D , 'file' , 'ascend_degreeUnderreplication.txt');

DS = sortrows(DS , 'percent_unreplicated_not_trimmed_cdc20' , 'descend');
D = DS(: , {'ORF'});
export(D , 'file' , 'descend_degreeUnderreplication.txt');

DS = sortrows(DS , 'percent_underreplicated_cdc20_not_trimmed_DM_dist' , 'ascend');
D = DS(: , {'ORF'});
export(D , 'file' , 'ascend_degreeUnderreplication_DM.txt');

DS = sortrows(DS , 'percent_underreplicated_cdc20_not_trimmed_DM_dist' , 'descend');
D = DS(: , {'ORF'});
export(D , 'file' , 'descend_degreeUnderreplication_DM.txt');

%%
data_mRNA = DS.max_mRNA;
data_PROseq = DS.max_PROseq;
figure; 
subplot(2,1,1); hold on; grid on; 
hist(data_mRNA , [nanmin(data_mRNA) : 1 : nanmax(data_mRNA)]);
xlim([0 10]);

subplot(2,1,2); hold on; grid on; 
hist(data_PROseq , [nanmin(data_PROseq) : 0.5 : nanmax(data_PROseq)]);
xlim([0 5]);



%%
DS.mRNA_bin = NaN(length(DS) , 1);
DS.PROseq_bin = NaN(length(DS) , 1);
for I = 1:length(DS)
    if DS.median_mRNA(I) <= quantile(DS.median_mRNA , .05)
        DS.mRNA_bin(I) = 1;
    elseif DS.median_mRNA(I) <= quantile(DS.median_mRNA , .95)
        DS.mRNA_bin(I) = 2;
    else
        DS.mRNA_bin(I) = 3;
    end
    
    if DS.median_PROseq(I) <= quantile(DS.median_PROseq , .05)
        DS.PROseq_bin(I) = 1;
    elseif DS.median_PROseq(I) <= quantile(DS.median_PROseq , .95)
        DS.PROseq_bin(I) = 2;
    else
        DS.PROseq_bin(I) = 3;
    end
end
%%
idx_orf = find(strcmp(DS.TYPE , 'ORF')); 
figure; 
clrs1 = pink(56);
clrs2 = summer(56);
subplot(2,2,1); hold on; grid on; set(gca , 'FontSize' , 12);
for I = 1:3
    idx = find(DS.mRNA_bin == I & strcmp(DS.TYPE , 'ORF'));
    data = DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx); data = data(~isnan(data));
    scatter(repmat(I , 1 , length(data))+rand(1 , length(data))*0.3-0.15 , data , 20 , clrs1(6*I+4,:) , ...
        'filled');
end
h = boxplot(DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx_orf) , DS.mRNA_bin(idx_orf) , 'color' , [.2 .2 .2] , 'symbol' , '');
set(h , 'LineWidth' , 1.9);
set(gca , 'Xtick' , [1:3] , 'XtickLabel' , {'<= 5% Q' , '5-95%' , '> 95% Q'});
xlabel('mRNAseq'); ylabel('% unreplicated');
ylim([-15 35]);

subplot(2,2,2); hold on; grid on; set(gca , 'FontSize' , 12);
for I = 1:3
    idx = find(DS.mRNA_bin == I & strcmp(DS.TYPE , 'ORF'));
    data = DS.percent_unreplicated_not_trimmed_cdc20_DM(idx); data = data(~isnan(data));
    scatter(repmat(I , 1 , length(data))+rand(1 , length(data))*0.3-0.15 , data , 20 , clrs1(6*I+4,:) , ...
        'filled');
end
h = boxplot(DS.percent_unreplicated_not_trimmed_cdc20_DM(idx_orf) , DS.mRNA_bin(idx_orf) , 'color' , [.2 .2 .2] , 'symbol' , '');
set(h , 'LineWidth' , 1.9);
set(gca , 'Xtick' , [1:3] , 'XtickLabel' , {'<= 5% Q' , '5-95%' , '> 95% Q'});
xlabel('mRNAseq'); ylabel('% unreplicated DM');
ylim([-20 20]);

subplot(2,2,3); hold on; grid on; set(gca , 'FontSize' , 12);
for I = 1:3
    idx = find(DS.PROseq_bin == I & strcmp(DS.TYPE , 'ORF'));
    data = DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx); data = data(~isnan(data));
    scatter(repmat(I , 1 , length(data))+rand(1 , length(data))*0.3-0.15 , data , 20 , clrs2(8*I+4,:) , ...
        'filled');
end
h = boxplot(DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx_orf) , DS.PROseq_bin(idx_orf) , 'color' , [.2 .2 .2] , 'symbol' , '');
set(h , 'LineWidth' , 1.9);
set(gca , 'Xtick' , [1:3] , 'XtickLabel' , {'<= 5% Q' , '5-95%' , '> 95% Q'});
xlabel('PROseq'); ylabel('% unreplicated');
ylim([-20 70]);

subplot(2,2,4); hold on; grid on; set(gca , 'FontSize' , 12);
for I = 1:3
    idx = find(DS.PROseq_bin == I & strcmp(DS.TYPE , 'ORF'));
    data = DS.percent_unreplicated_not_trimmed_cdc20_DM(idx); data = data(~isnan(data));
    scatter(repmat(I , 1 , length(data))+rand(1 , length(data))*0.3-0.15 , data , 20 , clrs2(8*I+4,:) , ...
        'filled');
end
h = boxplot(DS.percent_unreplicated_not_trimmed_cdc20_DM(idx_orf) , DS.PROseq_bin(idx_orf) , 'color' , [.2 .2 .2] , 'symbol' , '');
set(h , 'LineWidth' , 1.9);
set(gca , 'Xtick' , [1:3] , 'XtickLabel' , {'<= 5% Q' , '5-95%' , '> 95% Q'});
xlabel('PROseq'); ylabel('% unreplicated DM');
ylim([-20 20]);

print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Manuscript/Fig__Expression_median__VS_underreplication');





