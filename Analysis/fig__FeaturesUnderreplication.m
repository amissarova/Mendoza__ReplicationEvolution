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





