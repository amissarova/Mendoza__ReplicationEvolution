
%%
cd ~/Develop/Mendoza__ReplicationEvolution/
addpath(genpath('~/Develop/matlab'));
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_CLB2.mat');
DS = sortrows(DS , {'ReplicateNum' , 'Strain' , 'Condition' , 'TimePoint'});

%% Figure for the main body and Supp: for each replicate of y632: generate vector of pairs (GalRaf and Raf) of distributions (ksdensity) along time

% violet for GalRaf, Green for Raf
clrs2 = lines(6); clrs1 = summer(12); clrs_set = [clrs2(4,:) ; clrs1(4,:)];
for Z = 1:4
    figure('units','centimeters','position',[4 4 5 20]);
    % induce different final timepoints coz replicate 4 was a bit slower 
    if Z == 4
        idx = find(DS.TimeH <= 1.7 & DS.Strain == 632 & DS.ReplicateNum == Z & DS.TimeH >= 0);
    else
        idx = find(DS.TimeH <= 1.5 & DS.Strain == 632 & DS.ReplicateNum == Z & DS.TimeH >= 0);
    end
    unq_tp = unique(DS.TimePoint(idx));
    for I = 1:length(unq_tp)
        subplot(length(unq_tp) , 1 , I); hold on; grid on; set(gca , 'FontSize' , 10);
        idx = find(DS.TimePoint == unq_tp(I) & DS.Strain == 632 & DS.ReplicateNum == Z);
        for J = 1:length(idx)
            data = DS.trimmed_data{idx(J)};
            [y,x] = ksdensity(data , [2^13-1000:200:2^15]);
            plot(x , y/sum(y) , 'LineWidth' , 1.5 , 'color' , clrs_set(J,:) ); 
        end
        %title(sprintf('%d minutes after HU release' , DS.TimeMin(idx(1))));
    end
    save_name = sprintf('stat_ksdensity__rep%d' , Z);
    print('-dpng' , save_name);
end

%% Figure for the Supp: for y49 = generate vector of pairs (GalRaf and Raf) of distributions (ksdensity) along time
% violet for GalRaf, Green for Raf
clrs2 = lines(6); clrs1 = summer(12); clrs_set = [clrs2(4,:) ; clrs1(4,:)];
figure('units','centimeters','position',[4 4 5 20]);
idx = find(DS.TimeH <= 1.5 & DS.Strain == 49 & DS.TimeH >= 0);
unq_tp = unique(DS.TimePoint(idx));
for I = 1:length(unq_tp)
	subplot(length(unq_tp) , 1 , I); hold on; grid on; set(gca , 'FontSize' , 10);
	idx = find(DS.TimePoint == unq_tp(I) & DS.Strain == 49);
	for J = 1:length(idx)
        data = DS.trimmed_data{idx(J)};
        [y,x] = ksdensity(data , [2^13-1000:200:2^15]);
        plot(x , y/sum(y) , 'LineWidth' , 1.5 , 'color' , clrs_set(J,:) ); 
    end
	%title(sprintf('%d minutes after HU release' , DS.TimeMin(idx(1))));
end
save_name = 'stat_ksdensity__y49';
print('-dpng' , save_name);

%% additional for Lucas: scatter of G2-peaks, x-axis -- GalRaf; y-axis -- Raf
% can change timepoints if you wish
clrs = parula(16);
figure; hold on; grid on; 
for Z = 1:4
    G2_data = NaN(2,1);
    for J = 1:2
        if Z == 1
            idx = find(DS.TimeMin == 45 & DS.ReplicateNum == Z & DS.Strain == 632 & strcmp(DS.Condition , unq_condition{J}) );
        elseif Z < 4
            idx = find(DS.TimeMin == 60 & DS.ReplicateNum == Z & DS.Strain == 632 & strcmp(DS.Condition , unq_condition{J}) );
        else
            idx = find(DS.TimeMin == 80 & DS.ReplicateNum == Z & DS.Strain == 632 & strcmp(DS.Condition , unq_condition{J}) );
        end
        D = DS(idx , :);
        data = D.trimmed_data{1};
        data = data(data > 2^14);
        G2_data(J) = modefit(data , 0 , [2^14:200:2^15]);
    end
    scatter(G2_data(1) , G2_data(2) , 100 , clrs(4*Z-2,:) , 'filled' , 'Display' , sprintf('Rep #%d, TimeMin %d' , Z , D.TimeMin(1)));
end
G2_data = NaN(2,1);
for J = 1:2
    idx = find(DS.TimeMin == 60 & DS.Strain == 49 & strcmp(DS.Condition , unq_condition{J}) );
	D = DS(idx , :);
	data = D.trimmed_data{1};
	data = data(data > 2^14);
	G2_data(J) = modefit(data , 0 , [2^14:200:2^15]);
end
scatter(G2_data(1) , G2_data(2) , 100 , 'r' , 'filled' , 'Display' ,sprintf( 'y49, TimeMin %d' , D.TimeMin(1)));
plot([16000 26000] , [16000 26000] , 'LineWidth' , 2 , 'color' ,[.2 .2 .2] , 'LineStyle' , '--' , ...
    'Display' , 'X==Y');
xlabel('GalRaf'); ylabel('Raf'); legend('location' , 'nw'); title('G2 peak');
print('-dpng' , 'G2_across_replicates_early_timepoints');






