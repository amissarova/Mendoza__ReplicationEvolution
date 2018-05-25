%% Pipeline to generate panels for fig2 and potential Supplements for fig 2
addpath(genpath('~/Develop/matlab'));
cd ~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/

%% Fig 2A: example for one chromosome, under-replication VS distance to the end (pMet3-cdc20 and dbf2-2)
load('~/Develop/Mendoza__ReplicationEvolution/Data/MutantChr200_Michi.mat');
DS = G;
chr_num = 5; K = 75; unq_mutant = {'cdc20' , 'dbf2'};
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
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2A' , '-r300');

%% Fig 2B: histograms of under-replication, binned by replication timing

load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp_new.mat');
clrs1 = summer(12); clrs2 = pink(12); clrs3 = lines(6); clrs4 = winter(12);
clrs_set = [clrs1(6,:) ; clrs4(7,:) ; clrs3(4,:) ; clrs4(7,:)];
figure('units','centimeters','position',[5 5 10 10]); hold on; grid on; set(gca , 'FontSize' , 15);
data = DS.Trep_spline; 
%trep_thresh = [nanmin(data) quantile(data , .25) quantile(data , .5) quantile(data , .75) nanmax(data)];
trep_thresh = [0 30 45 65];
legend_names = {'<30 minutes' , '30-45 minutes' , '>45 minutes'};

for I = 1:length(trep_thresh)-1
    idx = find(DS.Trep_spline >= trep_thresh(I) & DS.Trep_spline < trep_thresh(I+1) );
    if ~isempty(idx)
        [y,x] = ksdensity(DS.percent_underreplicated_cdc20(idx)*100 , [-20:1:70]);
        plot(x,y/sum(y) , 'LineWidth' , 2.5 , 'color' , clrs_set(I,:));
    end
end
xlim([-20 70]);
legend(legend_names , 'location' , 'East');
% ylabel('Probability density'); xlabel('Under-replication, %');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2B' , '-r300');

%% fig 2C: for all 200bp-windows, G4-quadruplexes, absolute under-rep and DM for poor and rich regions

load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp_new.mat');
legend_titles = {'G4-poor' , 'G4-rich'};
clrs1 = lines(6); clrs2 = parula(12); clrs_set = [clrs1(4,:) ; clrs2(6,:)];
fh = figure('units','centimeters','position',[5 5 4 8 ]); hold on; grid on; set(gca , 'FontSize' , 10);
thresh_G4 = 80;
idx = DS.G4 > thresh_G4;
h1 = boxplot(DS.percent_underreplicated_cdc20*100 , idx ,'symbol','','color' , [.2 .2 .2]);
h = findobj(gca,'Tag','Box');set(h1 , 'LineWidth' , 1.9);
for j=1:length(h)
    patch(get(h(2-j+1),'XData'),get(h(2-j+1),'YData'), clrs_set(2-j+1,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
[~,p] = ttest2(DS.percent_underreplicated_cdc20(DS.G4 <= thresh_G4) , ...
    DS.percent_underreplicated_cdc20(DS.G4 > thresh_G4) );
%title(sprintf('P = %.4f.' , p));
ylim([-15 85]);
ylabel('Under-replication, %');
set(gca,'xtick',[]);
legend('location', 'SouthOutside');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2C' , '-r300');

%% fig 2D: High and low transcription rate have higher fraction of regions that are underreplicated
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp_new.mat');
figure('units','centimeters','position',[5 5 9 12 ]); 
subtel = [1 0]; 
%underrep = [10 15 20];
underrep = 15;
K = 0;
xticks = {'<= 4.48' , '4.48 - 5.97' , '5.97 - 10.44' , '10.44 - 14.92' , '>= 14.92'};
clrs1 = summer(12); clrs2 = lines(6); clrs3 = winter(12);
clrs_set = [clrs1(4,:) ; clrs2(4,:) ; clrs3(7,:)];
for I2 = 1:length(underrep)
	for I1 = 1:2
        K = K + 1;
        subplot(2 , length(underrep) , K); hold on; grid on;
        idx = find(DS.subtelomeric_bool == subtel(I1) & ~isnan(DS.max_PROseq));
        D = DS(idx , :);
        unq_bin = unique(D.PROseq_bin);
        data = NaN(length(unq_bin) , 1);
        for J = 1:length(unq_bin)
            idx = find(D.PROseq_bin == unq_bin(J));
            data(J) = nanmean(D.percent_underreplicated_cdc20(idx)*100 > underrep(I2))*100;
        end
        plot(unq_bin , data , 'LineWidth' , 1.5 , 'color' , clrs_set(I2,:));
        scatter(unq_bin , data , 30 , clrs_set(I2,:) , 'filled' );
        set(gca , 'Xtick' , [1:5] );
%         if K == 2 
%             title('Subtelomeric under-replicated regions');
%         elseif K == 5
%             title('Not subtelomeric under-replicated regions');
%         end
    end
end
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2D' , '-r300');

%% fig 2D LEGEND:
figure('units','centimeters','position',[5 5 16 8 ]); hold on;
set(gca , 'FontSize' , 20);
clrs1 = summer(12); clrs2 = lines(6); clrs3 = winter(12);
clrs_set = [clrs1(4,:) ; clrs2(4,:) ; clrs3(7,:)];

for I2 = 1:length(underrep)
    plot([0 0 ] , [1 1] , 'LineWidth' , 3 , 'color' , clrs_set(I2,:) , 'Display' , sprintf('%d%%' , underrep(I2)) );
end
legend('location' , 'best');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2D_LEGEND' , '-r300');
 

%% Don't run !!! for each GO ontology - save name, avg&Std -- dist to the end, avg&Std -- under-rep, and sample size

load('~/Develop/Data/YeastGO.mat');
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features_new.mat');
idx = find(strcmp(DS.TYPE , 'ORF')); DS = DS(idx , :);
T = struct2ds(YeastGO);
T.SampleSize = NaN(length(T) , 1);
T.mean_dist_to_the_end = NaN(length(T) , 1);
T.cv_dist_to_the_end = NaN(length(T) , 1);
T.mean_underreplication = NaN(length(T) , 1);
T.cv_underreplication = NaN(length(T) , 1);

for I = 1:length(T)
    temp_orf = T.ORF{I};
    idx = NaN(length(temp_orf) , 1);
    for J = 1:length(temp_orf)
        temp_idx = find(strcmp(DS.ORF , temp_orf{J}));
        if ~isempty(temp_idx)
            idx(J) = temp_idx;
        end
    end
    idx = idx(~isnan(idx));
    T.SampleSize(I) = length(idx);
    T.mean_dist_to_the_end(I) = nanmean(DS.dist_to_the_end_kb(idx));
    T.cv_dist_to_the_end(I) = nanstd(DS.dist_to_the_end_kb(idx))/nanmean(DS.dist_to_the_end_kb(idx));
    
    T.mean_underreplication(I) = nanmean(DS.percent_underreplicated_cdc20(idx)*100);
    T.cv_underreplication(I) = nanstd(DS.percent_underreplicated_cdc20(idx)*100)/nanmean(DS.percent_underreplicated_cdc20(idx)*100);
end
save('~/Develop/Mendoza__ReplicationEvolution/Data/YeastGOwUnderrep.mat' , 'T');

%% fig 2E: for each GO: dist to the end VS under-replication
load('~/Develop/Mendoza__ReplicationEvolution/Data/YeastGOwUnderrep.mat');
figure('units','centimeters','position',[5 5 13 13]); hold on; grid on; set(gca , 'FontSize' , 16)
idx = find(T.SampleSize > 5);
scatter(T.mean_dist_to_the_end(idx) , T.mean_underreplication(idx) , 30 , [.5 .5 .5] , 'filled' , ...
    'MarkerFaceAlpha' , .5);
idx = find(T.SampleSize > 5 & T.mean_underreplication > 20);
clrs1 = parula(12); clrs2 = summer(12); clrs3 = lines(6); clrs4 = spring(12); clrs5 = pink(12); clrs6 = hot(12);
clrs_set = [clrs1(3,:) ; clrs3(4,:) ; clrs1(10,:) ; clrs2(4,:) ; clrs4(5,:) ; ...
    clrs1(7,:) ; clrs5(5,:) ; clrs6(3,:) ; clrs2(1,:)];
for I = 1:length(idx)
    %errorbar(T.mean_dist_to_the_end(idx(I)) , T.mean_underreplication(idx(I)) , ...
    %    T.cv_underreplication(idx(I)) , 'LineWidth' , 1 , 'color' , [.3 .3 .3]);
    %herrorbar(T.mean_dist_to_the_end(idx(I)) , T.mean_underreplication(idx(I)) , ...
    %    T.cv_dist_to_the_end(idx(I)) );
    scatter(T.mean_dist_to_the_end(idx(I)) , T.mean_underreplication(idx(I)) , 100 , clrs_set(I,:) , 'filled','d' ,...
        'Display' , T.FamilyName{idx(I)} );
end
%legend('location' , 'EastOutside');
ylim([-5 40]); xlim([0 450]);
xlabel('Distance to telomeric start, kb');
ylabel('Under-replication, %');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2E' , '-r300');

%% fig 2E LEGEND: all highlighted GOs
load('~/Develop/Mendoza__ReplicationEvolution/Data/YeastGOwUnderrep.mat');
figure('units','centimeters','position',[5 5 13 13]); hold on; grid on; set(gca , 'FontSize' , 16)
idx = find(T.SampleSize > 5 & T.mean_underreplication > 20);
clrs1 = parula(12); clrs2 = summer(12); clrs3 = lines(6); clrs4 = spring(12); clrs5 = pink(12); clrs6 = hot(12);
clrs_set = [clrs1(3,:) ; clrs3(4,:) ; clrs1(10,:) ; clrs2(4,:) ; clrs4(5,:) ; ...
    clrs1(7,:) ; clrs5(5,:) ; clrs6(3,:) ; clrs2(1,:)];
for I = 1:length(idx)

    scatter(T.mean_dist_to_the_end(idx(I)) , T.mean_underreplication(idx(I)) , 100 , clrs_set(I,:) , 'filled','d' ,...
        'Display' , T.FamilyName{idx(I)} );
end
legend('location' , 'best');
ylim([-5 40]); xlim([0 450]);
xlabel('Distance to telomeric start, kb');
ylabel('Under-replication, %');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2E_LEGEND' , '-r300');


%% Fig 2F: under-replication (absolute and DM): ORF VS tRNA, transposons and fragile sites

load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features_new.mat');
DS.type_general = cell(length(DS) , 1);
for I = 1:length(DS)
    temp_type = DS.TYPE{I};
    if ~strcmp(temp_type , 'interstitial_del_dup') & ~strcmp(temp_type , 'terminal_del_dup')
        DS.type_general{I} = temp_type;
    else
        DS.type_general{I} = 'del_dup';
    end
end
type_of_interest = {'ORF' , 'transposable_element_gene', 'del_dup'};
data = NaN(length(DS) , length(type_of_interest));
for I = 1:length(type_of_interest)
    idx = find(strcmp(DS.type_general , type_of_interest{I}));
    data(1:length(idx) , I) = DS.percent_underreplicated_cdc20(idx)*100;
end
clrs1 = summer(12); clrs2 = hot(12); clrs3 = parula(12); clrs4 = lines(6); clrs5 = spring(12);
clrs_set = [.65 .65 .65 ; clrs1(3,:) ; clrs2(3,:) ; clrs3(5,:) ; ...
    clrs4(4,:) ; clrs3(10,:) ; clrs3(3,:) ; clrs5(3,:)];
legend_titles = {'ORFs'  , 'transposable elements' , ...
    'fragile sites'};
figure('units','centimeters','position',[5 5 10 10]); hold on; grid on; set(gca , 'FontSize' , 10);
h1 = boxplot(data , 'color' , [.2 .2 .2] , 'symbol' , '');
set(h1 , 'LineWidth' , 1.5);
h = findobj(gca,'Tag','Box'); 
N = length(h);
for j=1:N
    patch(get(h(N-j+1),'XData'),get(h(N-j+1),'YData'), clrs_set(j,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
legend('location' , 'nw');
ylim([-10 60]);
set(gca , 'Xtick' , '');
%ylabel('Under-replication, %');

print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2F' , '-r300');

%% Fig 2G: essential genes VS dist to the end and under-replication
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features_new.mat');
idx = find(strcmp(DS.TYPE , 'ORF')); D = DS(idx , :);
for I = 1:length(D)
    if isempty(D.Phenotype_Observable{I})
        D.Phenotype_Observable{I} = 'Viable';
    end
end
D = D(: , {'ORF' , 'percent_underreplicated_cdc20' , 'Phenotype_Observable'});
%D = unique(D);
figure('units','centimeters','position',[5 5 5 8]); hold on; grid on;
h1 = boxplot(D.percent_underreplicated_cdc20*100 , D.Phenotype_Observable , 'color' , [.2 .2 .2] , 'symbol' , '');
h = findobj(gca,'Tag','Box');set(h1 , 'LineWidth' , 1.5);
clrs1 = spring(12); clrs2 = winter(12); clrs_set = [clrs1(6,:) ; clrs2(6,:)];
N = 2;
legend_names = {'Essential' , 'Not essential'};
for j=1:length(h)
    patch(get(h(N-j+1),'XData'),get(h(N-j+1),'YData'), clrs_set(N-j+1,:) ,'FaceAlpha' , .5 , 'Display' , ...
        legend_names{j});
end
legend('location' , 'SouthOutside');
set(gca , 'Xtick' , []);
ylim([-10 20]);
ylabel('under-replication , %');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2G' , '-r300');
%% Fig 2H: under-replication against relative fitness from Baryshnikova
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features_new.mat');
idx = find(strcmp(DS.TYPE , 'ORF')); D = DS(idx , :);
D = D(: , {'ORF' , 'dist_to_the_end_kb' , 'percent_underreplicated_cdc20'});
T = dataset('file' , '~/Develop/ExternalData/Baryshnikova10.tab');
T = T(: , {'ORF' , 'fitness'}); T = unique(T);
D = join(D , T , 'Type' , 'left' , 'Keys' , 'ORF' , 'MergeKeys' , true);
idx = find(~isnan(D.fitness)); D = D(idx , :);
clrs = winter(12);
figure('units','centimeters','position',[5 5 5 8]); hold on; grid on;
scatter(D.percent_underreplicated_cdc20*100 , D.fitness , 20 , clrs(7,:) , 'filled' , 'MarkerFaceAlpha' , .2);
xlim([-15 60]);
xlabel('Under-replication , %');
ylabel('Relative to WT fitness');
D = sortrows(D , 'percent_underreplicated_cdc20');
%S = CalcDM_Newman06( D.percent_underreplicated_cdc20*100 , D.fitness  , 100 , 0 );
%plot(D.percent_underreplicated_cdc20*100 , S.DM_ypred , 'LineWidth' , 3 , 'color' , clrs(5,:));    
R = corrcoef(D.percent_underreplicated_cdc20*100 , D.fitness); 
title(sprintf('R = %.1f.' , R(1,2)));
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2H' , '-r300');

%% 2H - boxplot option: grid how in 2I
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features_new.mat');
idx = find(strcmp(DS.TYPE , 'ORF')); D = DS(idx , :);
D = D(: , {'ORF' , 'dist_to_the_end_kb' , 'percent_underreplicated_cdc20'});
T = dataset('file' , '~/Develop/ExternalData/Baryshnikova10.tab');
T = T(: , {'ORF' , 'fitness'}); T = unique(T);
D = join(D , T , 'Type' , 'left' , 'Keys' , 'ORF' , 'MergeKeys' , true);
D.underrep_bin = NaN(length(DS) , 1);
for I = 1:length(D)
    if D.percent_underreplicated_cdc20(I)*100 <= 20
        D.underrep_bin(I) = 1;
    elseif D.percent_underreplicated_cdc20(I)*100 <= 30
        D.underrep_bin(I) = 2;
    elseif D.percent_underreplicated_cdc20(I)*100 <= 40
        D.underrep_bin(I) = 3;
    else
        D.underrep_bin(I) = 4;
    end
end


clrs1 = winter(8);
figure('units','centimeters','position',[5 5 9 9]);
legend_titles = {'<20' , '20-30' , '30-40' , '>40'};
hold on; grid on; set(gca , 'FontSize' , 12);
h1 = boxplot(D.fitness , D.underrep_bin , 'color' , [.2 .2 .2] , 'symbol' , '');
%ylim([0 .175]);
set(h1 , 'LineWidth' , 1.5);
h = findobj(gca,'Tag','Box'); set(h1 , 'LineWidth' , 1.5);
for j=1:4
    patch(get(h(4-j+1),'XData'),get(h(4-j+1),'YData'), clrs1(2*j,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
set(gca , 'Xtick' , []); 
%xlabel('% unreplicated cells');
ylabel('Relative to WT fitness');
legend('location' , 'SouthOutside');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2H__boxplot' , '-r300');

%% Fig 2I: boxplot for SNPs and InDel frequency (for 200bp-windows) for different under-replication bins

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
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2I' , '-r300');

%% Fig 2J: boxplot for gene preservation (for ORFs) for different under-replication bins

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
legend_titles = {'<20' , '20-30' , '30-40' , '>40'};
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
legend('location' , 'SouthOutside');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2J' , '-r300');
 
%% AlsuSup1: all chromosomes, all repeatable replicates of cdc20


load('~/Develop/Mendoza__ReplicationEvolution/Data/MutantChr200_Michi.mat');
DS = G;
unq_mutant = {'cdc20'}; chr_nums = [1:16];
clrs1 = hot(12); clrs2 = lines(6); clrs3 = parula(12); clrs4 = summer(12); clrs5 = spring(12);
clrs_set = [clrs1(3,:) ; clrs2(4,:) ; clrs3(2,:) ; clrs3(10,:) ; clrs4(2,:) ; clrs3(6,:)];
for I = 1:length(unq_mutant)
    figure('units','centimeters','position',[5 5 15 15]);
    for J = 1:length(chr_nums)
        subplot(4,4,J); hold on; grid on; set(gca , 'FontSize' , 6);
        idx = find( strcmp(DS.MutantID , unq_mutant{I}) & DS.chr_num == chr_nums(J) );
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

print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/SupAlsu2' , '-r300');


%% get GOs for highly under-replicated
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features_new.mat');
idx = find(strcmp(DS.TYPE , 'ORF')); D = DS(idx , :);
D = sortrows(D , 'percent_underreplicated_cdc20' , 'descend');
D = D(: , 'ORF');

export(D , 'file' , 'ORFs_underrep_descent.tab')

%% additional for Michi: heatshock genes against others

orfs = {'YAL005C' , 'YLL024C' , 'YBL075C' , 'YER103W' , 'YDL229W' , 'YNL209W' , ...
    'YLL026W' , 'YLL026W' , 'YDR258C' , 'YLR259C'};

load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features_new.mat');
idx = find(strcmp(DS.TYPE , 'ORF')); DS = DS(idx , :);
DS.heatshock_bool = zeros(length(DS) , 1);
for I = 1:length(DS)
    if ~isempty(find(strcmp(DS.ORF{I} , orfs)))
        DS.heatshock_bool(I) = 1;
    end
end
figure; hold on; grid on;
h1 = boxplot(DS.percent_underreplicated_cdc20 , DS.heatshock_bool , 'color' , [.2 .2 .2] , 'symbol' , '');
set(h1 , 'LineWidth' , 1.5);
h = findobj(gca,'Tag','Box'); set(h1 , 'LineWidth' , 1.5);
N = length(h);
for j=1:N
    patch(get(h(N-j+1),'XData'),get(h(N-j+1),'YData'), clrs1(2*j,:) ,'FaceAlpha',.5 );
end
xlabel('Heatshock bool');
ylabel('Under-replication');
legend('location' , 'SouthOutside');
    







