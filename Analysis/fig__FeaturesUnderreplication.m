%% Figures 6

%% all figures regarding subtelomeric underreplication
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features_new.mat');
idx = find(~isnan(DS.percent_unreplicated_not_trimmed_cdc20) & DS.percent_unreplicated_not_trimmed_cdc20 < 10000 & ...
    DS.percent_unreplicated_not_trimmed_cdc20 > -10000);
DS = DS(idx , :);

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
    scatter(repmat(I,length(temp_data) , 1)+randn(length(temp_data) , 1)*0.036-0.018 , temp_data ,...
        10 , clrs_set(I,:) , 'filled');
end
h = boxplot(data , 'color' , [.2 .2 .2] , 'symbol' , '');
set(h , 'LineWidth' , 1.5);
for I = 2:length(type_of_interest)
    temp_data_ORF = data(:,1);
    temp_data = data(:,I);
    [~,p] = ttest2(temp_data_ORF(~isnan(temp_data_ORF)) , temp_data(~isnan(temp_data)));
    if p < .005 & nanmean(temp_data) > nanmean(temp_data_ORF)
        scatter(I-0.1,65, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
        scatter(I+0.1,65, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
    elseif p < .05 & nanmean(temp_data) > nanmean(temp_data_ORF)
        scatter(I,65, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
    end
end
ylim([-20 70]);
set(gca , 'Xtick' , '');
%legend('ORF' , 'ARS','centromeres' , 'tRNA' , 'transposons' , ...
%    'LOH-loci' , 'Interstitial dup/del' , 'Terminal dup/del' , 'location' , 'EastOutside');
title('Metaphase');
ylabel('Degree of underreplication');



subplot(2,2,2); hold on; grid on; set(gca , 'FontSize' , 12);
for I = 1:length(type_of_interest)
    temp_data = data_dbf(:,I);
    temp_data = temp_data(~isnan(temp_data));
    scatter(repmat(I,length(temp_data) , 1)+randn(length(temp_data) , 1)*0.036-0.018 , temp_data ,...
        10 , clrs_set(I,:) , 'filled');
end
h = boxplot(data_dbf , 'color' , [.2 .2 .2] , 'symbol' , '');
set(h , 'LineWidth' , 1.5);
for I = 2:length(type_of_interest)
    temp_data_ORF = data_dbf(:,1);
    temp_data = data_dbf(:,I);
    [~,p] = ttest2(temp_data_ORF(~isnan(temp_data_ORF)) , temp_data(~isnan(temp_data)));
    if p < .005 & nanmean(temp_data) > nanmean(temp_data_ORF)
        scatter(I-0.1,55, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
        scatter(I+0.1,55, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
    elseif p < .05 & nanmean(temp_data) > nanmean(temp_data_ORF)
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
    scatter(repmat(I,length(temp_data) , 1)+randn(length(temp_data) , 1)*0.036-0.018 , temp_data ,...
        10 , clrs_set(I,:) , 'filled');
end
h = boxplot(data_DM , 'color' , [.2 .2 .2] , 'symbol' , '');
set(h , 'LineWidth' , 1.5);
for I = 2:length(type_of_interest)
    temp_data_ORF = data_DM(:,1);
    temp_data = data_DM(:,I);
    [~,p] = ttest2(temp_data_ORF(~isnan(temp_data_ORF)) , temp_data(~isnan(temp_data)));
    if p < .005 & nanmean(temp_data) > nanmean(temp_data_ORF)
        scatter(I-0.1,35, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
        scatter(I+0.1,35, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
    elseif p < .05 & nanmean(temp_data) > nanmean(temp_data_ORF)
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
    scatter(repmat(I,length(temp_data) , 1)+randn(length(temp_data) , 1)*0.036-0.018 , temp_data ,...
        10 , clrs_set(I,:) , 'filled');
end
h = boxplot(data_dbf_DM , 'color' , [.2 .2 .2] , 'symbol' , '');
set(h , 'LineWidth' , 1.5);
for I = 2:length(type_of_interest)
    temp_data_ORF = data_dbf_DM(:,1);
    temp_data = data_dbf_DM(:,I);
    [~,p] = ttest2(temp_data_ORF(~isnan(temp_data_ORF)) , temp_data(~isnan(temp_data)));
    if p < .005 & nanmean(temp_data) > nanmean(temp_data_ORF)
        scatter(I-0.1,35, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
        scatter(I+0.1,35, 10 , [.2 .2 .2] , '*' , 'LineWidth' , 1);
    elseif p < .05 & nanmean(temp_data) > nanmean(temp_data_ORF)
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
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Manuscript/Fig5A');

%%
idx = find(strcmp(DS.TYPE , 'ORF') & DS.chr_num == 1 & (DS.mRNA_minus < -1 | DS.mRNA_plus > 1));
D = DS(idx , :);
D.middle_point = (D.start_point + D.end_point)/2;
%%
figure; hold on; grid on;
plot(D.middle_point , D.mRNA_minus);
plot(D.middle_point , D.mRNA_plus);





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
idx = find(~isnan(DS.log2_mRNA_Nagalakshimi));
D = DS(idx , :);
D.bin_expression = NaN(length(D) , 1);
for I = 1:length(D)
    if D.log2_mRNA_Nagalakshimi(I) <= 2
        D.bin_expression(I) = 0;
    elseif D.log2_mRNA_Nagalakshimi(I) <= 4
        D.bin_expression(I) = 1;
	elseif D.log2_mRNA_Nagalakshimi(I) <= 6
        D.bin_expression(I) = 2;
	elseif D.log2_mRNA_Nagalakshimi(I) <= 8
        D.bin_expression(I) = 3; 
	elseif D.log2_mRNA_Nagalakshimi(I) <= 10
        D.bin_expression(I) = 4;    
    else
        D.bin_expression(I) = 5;
    end
end
%%
unq_bin = unique(D.bin_expression);
figure; 

subplot(2,2,1); hold on; grid on; set(gca , 'FontSize' , 12);
clrs = parula(4*length(unq_bin));
for I = 1:length(unq_bin)
    idx = find(D.bin_expression == unq_bin(I));
    scatter(repmat(I , length(idx) , 1)+randn(length(idx) , 1)*0.03-0.015 , D.percent_unreplicated_not_trimmed_cdc20_smooth(idx) , 30 , clrs(4*I-2,:) , 'filled' )
    boxplot(D.percent_unreplicated_not_trimmed_cdc20_smooth , D.bin_expression , ...
        'color' , [.2 .2 .2] , 'symbol' , '');
    data = DS.percent_unreplicated_not_trimmed_cdc20_smooth(strcmp(DS.TYPE , 'ORF'));
    [~,p] = ttest2(data , D.percent_unreplicated_not_trimmed_cdc20_smooth(idx));
    if p < .005
        scatter(I-0.1 , 25 , 30,[.2 .2 .2] , '*' , 'LineWidth' , 3);
        scatter(I+0.1 , 25 , 30,[.2 .2 .2] , '*' , 'LineWidth' , 3);
    elseif p < .05
        scatter(I , 25 , 30,[.2 .2 .2] , '*' , 'LineWidth' , 3);
    end
end
xlabel('Expression');
ylabel('Degree underrep, meta');
ylim([-20 30]);

subplot(2,2,2); hold on; grid on; set(gca , 'FontSize' , 12);
clrs = parula(4*length(unq_bin));
for I = 1:length(unq_bin)
    idx = find(D.bin_expression == unq_bin(I));
    scatter(repmat(I , length(idx) , 1)+randn(length(idx) , 1)*0.03-0.015 , D.percent_unreplicated_not_trimmed_cdc20_smooth(idx) , 30 , clrs(4*I-2,:) , 'filled' )
    boxplot(D.percent_unreplicated_not_trimmed_cdc20_smooth , D.bin_expression , ...
        'color' , [.2 .2 .2] , 'symbol' , '');
    data = DS.percent_unreplicated_not_trimmed_dbf2_smooth(strcmp(DS.TYPE , 'ORF'));
    [~,p] = ttest2(data , D.percent_unreplicated_not_trimmed_dbf2_smooth(idx));
    if p < .005
        scatter(I-0.1 , 25 , 30,[.2 .2 .2] , '*' , 'LineWidth' , 3);
        scatter(I+0.1 , 25 , 30,[.2 .2 .2] , '*' , 'LineWidth' , 3);
    elseif p < .05
        scatter(I , 25 , 30, [.2 .2 .2] , '*' , 'LineWidth' , 3);
    end
end
xlabel('Expression');
ylabel('Degree underrep, ana');
ylim([-25 30]);

subplot(2,2,3); hold on; grid on; set(gca , 'FontSize' , 12);
clrs = parula(4*length(unq_bin));
for I = 1:length(unq_bin)
    idx = find(D.bin_expression == unq_bin(I));
    scatter(repmat(I , length(idx) , 1)+randn(length(idx) , 1)*0.03-0.015 , D.percent_unreplicated_not_trimmed_cdc20_smooth(idx) , 30 , clrs(4*I-2,:) , 'filled' )
    boxplot(D.percent_underreplicated_cdc20_not_trimmed_DM_dist , D.bin_expression , ...
        'color' , [.2 .2 .2] , 'symbol' , '');
    data = DS.percent_underreplicated_cdc20_not_trimmed_DM_dist(strcmp(DS.TYPE , 'ORF'));
    [~,p] = ttest2(data , D.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx));
    if p < .005
        scatter(I-0.1 , 25 , 30, [.2 .2 .2] , '*' , 'LineWidth' , 3);
        scatter(I+0.1 , 25 , 30, [.2 .2 .2] , '*' , 'LineWidth' , 3);
    elseif p < .05
        scatter(I , 25 , 30, [.2 .2 .2] , '*' , 'LineWidth' , 3);
    end
end
xlabel('Expression');
ylabel('Degree underrep DM, meta');
ylim([-20 40]);

subplot(2,2,4); hold on; grid on; set(gca , 'FontSize' , 12);
clrs = parula(4*length(unq_bin));
for I = 1:length(unq_bin)
    idx = find(D.bin_expression == unq_bin(I));
    scatter(repmat(I , length(idx) , 1)+randn(length(idx) , 1)*0.03-0.015 , D.percent_unreplicated_not_trimmed_cdc20_smooth(idx) , 30 , clrs(4*I-2,:) , 'filled' )
    boxplot(D.percent_underreplicated_dbf2_not_trimmed_DM_dist , D.bin_expression , ...
        'color' , [.2 .2 .2] , 'symbol' , '');
	data = DS.percent_underreplicated_dbf2_not_trimmed_DM_dist(strcmp(DS.TYPE , 'ORF'));
    [~,p] = ttest2(data , D.percent_underreplicated_dbf2_not_trimmed_DM_dist(idx));
    if p < .005
        scatter(I-0.1 , 25 , 30, [.2 .2 .2] , '*' , 'LineWidth' , 3);
        scatter(I+0.1 , 25 , 30, [.2 .2 .2] , '*' , 'LineWidth' , 3);
    elseif p < .05
        scatter(I , 25 ,30,  [.2 .2 .2] , '*' , 'LineWidth' , 3);
    end
end
xlabel('Expression');
ylabel('Degree underrep DM, ana');
ylim([-25 30]);

%% how expression correlates with underreplication
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features.mat');
D = readtable('~/Develop/Mendoza__ReplicationEvolution/Data/ExternalData/Yeast_mRNA_abundance_Lipson_Causey.xls');
D = table2dataset(D);
DS = join(DS , D , 'Keys' , 'ORF' , 'Type' , 'Left' , 'MergeKeys' , true);
%
idx = find(strcmp(DS.TYPE , 'ORF') & ~isnan(DS.Avg) & ~isnan(DS.percent_unreplicated_not_trimmed_cdc20));
D = DS(idx , :);

%%
%quantile_expression_bins = [quantile(data , 0) quantile(data , 0.2) quantile(data , 0.4) quantile(data , 0.6) quantile(data , 0.8)  quantile(data , 1)];
figure; 
K = 5;
data = D.Avg;
quantile_expression_bins = [quantile(data , 0) ...
    quantile(data , 0.1)  quantile(data , 0.2)  ...
    quantile(data , 0.95) quantile(data , 0.99)  ...
    quantile(data , 1)];

clrs = parula(4*K);
subplot(2,1,1);hold on; grid on; set(gca , 'FontSize' , 12);
data = NaN(length(D) , K);
for I = 1:K
    idx = find(D.Avg >= quantile_expression_bins(I) ...
        & D.Avg < quantile_expression_bins(I+1));
    data(1:length(idx) , I) = D.percent_unreplicated_not_trimmed_cdc20(idx);
    scatter(repmat(I , length(idx) , 1)+randn(length(idx) , 1)*0.03 - 0.015 , ...
       D.percent_unreplicated_not_trimmed_cdc20(idx), 10 , clrs(4*I-2,:) , 'filled');
end
boxplot(data , 'color' , [.2 .2 .2] , 'symbol' , '' , 'notch' , 'on');
%plot([1:K] , nanmedian(data) , 'LineWidth' , 2 , 'color' , [.2 .2 .2]);
ylim([-20 50]);
set(gca , 'Xtick' , [1:5] , 'XtickLabel' , {'0-20%' , '20-40%' , '40-60%' , '60-80%' , '80-100%'});
xlabel('Expression, quantile');
ylabel('Degree underreplication');

subplot(2,1,2);hold on; grid on; set(gca , 'FontSize' , 12);
data = D.Avg;
quantile_expression_bins = [quantile(data , 0) ...
    quantile(data , 0.1)  quantile(data , 0.2)  ...
    quantile(data , 0.95) quantile(data , 0.99)  ...
    quantile(data , 1)];
data = NaN(length(D) , K);
for I = 1:K
    idx = find(D.Avg >= quantile_expression_bins(I) ...
        & D.Avg < quantile_expression_bins(I+1));
    data(1:length(idx) , I) = D.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx);
    scatter(repmat(I , length(idx) , 1)+randn(length(idx) , 1)*0.03 - 0.015 , ...
       D.percent_unreplicated_not_trimmed_cdc20(idx), 10 , clrs(4*I-2,:) , 'filled');
end
boxplot(data , 'color' , [.2 .2 .2] , 'symbol' , '' , 'notch' , 'on');
ylim([-20 60]);
set(gca , 'Xtick' , [1:5] , 'XtickLabel' , {'0-20%' , '20-40%' , '40-60%' , '60-80%' , '80-100%'});
xlabel('Expression, quantile');
ylabel('Degree underreplication DM');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Manuscript/Fig5C');
%%
K = 5;
data = D.Avg;
quantile_expression_bins = [quantile(data , 0) ...
    quantile(data , 0.2)  quantile(data , 0.4)  ...
    quantile(data , 0.6) quantile(data , 0.8)  ...
    quantile(data , 1)];


data = NaN(length(D) , K);
for I = 1:K
    idx = find(D.Avg >= quantile_expression_bins(I) ...
        & D.Avg < quantile_expression_bins(I+1));
    data(1:length(idx) , I) = D.percent_unreplicated_not_trimmed_cdc20(idx);

end
p_vals = NaN(K , K);
for I = 1:K-1
    for J = I+1:K
        [~,p] = ttest2(data(:,I) , data(:,J));
        p_vals(I,J) = p;
    end
end














