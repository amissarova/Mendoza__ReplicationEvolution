%% Different Versions of Figure regarding CDK inhibition --> finishing DNA replication

%% Plot: distance to the end VS under-rep, example for one chromosome (5, left arm)
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_Tsveti.mat');
DS = G;
chr_num = 5; K = 75; unq_mutant = {'M' , 'S'};
clrs1 = parula(12); clrs_set = [.6 .6 .6 ; clrs1(5,:)];
figure('units','centimeters','position',[5 5 8 8]);
hold on; grid on; set(gca , 'FontSize' , 12);
for I = 1:length(unq_mutant)
	idx = find( strcmp(DS.MutantID , unq_mutant{I}) & DS.chr_num == chr_num );
	data = NaN( length(idx) , K);
    for J = 1:length(idx)
        y = DS.unreplicated_cells_not_trimmed{idx(J)};
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
    set(h,'facealpha',.2);
end
xlim([0 K]);
ylim([-9 40]);
set(gca , 'Xtick' , [0 25 50 75]);
set(gca , 'Ytick' , [0 20 40 60]);
title(strcat ( DS.chr{idx(1)} , ', left'));
xlabel('Distance to the end, kbp');
ylabel('% unreplicated cells');
%print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig4__CDKinhibit_Tsveti/4_CDK_inhibit_1chr' , '-r300');

%% Plot: distance to the end VS under-rep, all chromosomal ends together
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_Tsveti.mat');
DS = G;
K = 75 ; unq_mutant = {'M' , 'S'};
clrs1 = parula(12); clrs_set = [.6 .6 .6 ; clrs1(5,:)];
figure('units','centimeters','position',[5 5 8 8]);
hold on; grid on; set(gca , 'FontSize' , 12);
legend_names = {'Metaphase arrest' , 'Metaphase arrest + CDK inhibition'};
for I = 1:length(unq_mutant)
    idx = find(strcmp(DS.MutantID , unq_mutant{I}));
    data = NaN( 2*length(idx) , K);	
    for J = 1:length(idx)
        y = DS.unreplicated_cells_not_trimmed{idx(J)};
        start_point_kb = DS.start_point_kb{idx(J)};
        for Z = 1:K
            idx_current_kb = find(start_point_kb == Z);
            data(2*J-1 , Z) = nanmedian(y(idx_current_kb));
            idx_current_kb = find(start_point_kb == nanmax(start_point_kb) - Z + 1);
            data(2*J , Z) = nanmedian(y(idx_current_kb));
        end  
    end
    

	x = [1:K];
	mean_y = nanmedian(data); std_y = nanstd(data);
	plot(x , mean_y , 'LineWidth' , 3 , 'color' , clrs_set(I,:) ,...
            'Display' , strcat( 'mean,' , unq_mutant{I} ));
    h1 = fill([x';flipud(x')],[mean_y'-std_y';flipud(mean_y'+std_y')], clrs_set(I,:) ,...
            'linestyle','none' ,  'Display' , strcat( '2xstd,' , unq_mutant{I} ) );
    set(h1,'facealpha',.2);
end
xlim([0 K]);
ylim([-9 40]);
set(gca , 'Xtick' , [0 25 50 75]);
set(gca , 'Ytick' , [0 20 40 60]);
title('All subtelomeric regions');
xlabel('Distance to the end, kbp');
ylabel('% unreplicated cells');
%print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig4__CDKinhibit_Tsveti/4_CDK_inhibit_all_chr' , '-r300');

%% Scatter, each dot -- one 200bp window, x-axis: underrep in M, y-axis: underrep in M + CDK inhibtion
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_Tsveti.mat');
DS = sortrows(DS , 'percent_underrep_Tsveti_M');
clrs = parula(12); 
idx = find(DS.percent_underrep_Tsveti_M >= 0);
D = DS(idx , :);
figure('units','centimeters','position',[5 5 8 8]); hold on; grid on; set(gca , 'FontSize' , 10);
scatter(D.percent_underrep_Tsveti_M , D.percent_underrep_Tsveti_S , 30 , clrs(5,:) , 'filled' , ...
    'MarkerFaceAlpha' , .2);
xlabel('Metaphase');
ylabel('Metaphase + CDK inhibition');
R = corrcoef(D.percent_underrep_Tsveti_M , D.percent_underrep_Tsveti_S );
[~,p] = ttest2(D.percent_underrep_Tsveti_M , D.percent_underrep_Tsveti_S );
title(sprintf('R = %.2f. p < 0.001.' , R(1,2)));
plot([0 50 ] , [0 50] , 'LineWidth' , 2 , 'LineStyle' , '--' , 'color' , [.2 .2 .2]);
xlim([0 50]);
%print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig4__CDKinhibit_Tsveti/4_CDK_inhibit_scatter' , '-r300');











