%% Different Versions of Figure regarding CDK inhibition --> finishing DNA replication

%% Plot: distance to the end VS under-rep, example for one chromosome (5, left arm)
K = 75; unq_mutant = {'M' , 'S'};
clrs1 = parula(12); clrs2 = lines(6) ; clrs_set = [clrs1(6,:) ; clrs1(10,:)]; 
str = cell(2,1);
load( '~/Develop/Mendoza__ReplicationEvolution/Data/MutantChr200_Tsveti1.mat') ;
D = G;
load( '~/Develop/Mendoza__ReplicationEvolution/Data/MutantChr200_Tsveti2.mat' );
D = [D; G];
figure('units','centimeters','position',[5 5 15 10]);
hold on; grid on; set(gca , 'FontSize' , 16);
count = 0;
%for I1 = 1:length(str)
    for I2 = 1:length(unq_mutant)
        idx = find( strcmp(D.MutantID , unq_mutant{I2}));
        if ~isempty(idx)
            count = count + 1;
            data = NaN( 2*length(idx) , K);
            for J = 1:length(idx)
                y = D.percent_underreplicated{idx(J)}*100;
                start_point_kb = D.start_point_kb{idx(J)};
                start_point_kb = round(start_point_kb);
                for Z = 1:K
                    idx_current_kb = find(start_point_kb == Z);
                    data(2*J-1 , Z) = nanmedian(y(idx_current_kb));
                    idx_current_kb = find(start_point_kb == nanmax(start_point_kb) - Z + 1);
                    data(2*J , Z) = nanmedian(y(idx_current_kb));
                end
            end
            x = [1:K];
            mean_y = nanmedian(data); std_y = nanstd(data);
            plot(x , mean_y , 'LineWidth' , 3 , 'color' , clrs_set(count,:));
            h = fill([x';flipud(x')],[mean_y'-std_y';flipud(mean_y'+std_y')], clrs_set(count,:) ,...
                'linestyle','none');
            set(h,'facealpha',.15);
        end
	end
%end
xlim([0 K]);
ylim([-9 40]);
set(gca , 'Xtick' , [0 25 50 75]);
set(gca , 'Ytick' , [0 20 40 60]);
title('All subtelomeric regions');
%xlabel('Distance to the end, kbp');
%ylabel('Under-replication, %');
%print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig4__CDKinhibit_Tsveti/4_CDK_inhibit_all_chr' , '-r600');
