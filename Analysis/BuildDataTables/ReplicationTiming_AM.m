%% first run ReplicationTiming to get dataset

%% alsu's violin style
pctile_thresholds = [20 25 30 35 40 45 ];
legend_text = [0 pctile_thresholds] ; 
Classes = zeros(height(T),1);
for I = 1:numel(pctile_thresholds)
    Classes(T.percent_underreplicated_cdc20*100 >= pctile_thresholds(I)) = I ; 
end
unq_class = unique(Classes); K = length(unq_class); clrs_set = parula(K);

fh = figure('units','centimeters','position',[5 5 20 20]); 

% Raghuraman01
subplot(2,2,1); hold on; grid on;
data = T.Raghuraman01;
data_violin = NaN(length(data) , K);
for I = 1:K
    idx = find(Classes == unq_class(I));
    data_violin(1:length(idx) , I) = data(idx);
end

violin(data_violin , 'facecolor',clrs_set,'edgecolor',[.25 .25 .25],'bw',2,'mc',[] , 'medc' , 'k');
title({'Raghuraman et. al. 2001' 'RM14-3a, \alpha-factor & cdc7-1) ; 6 reps'})
ylabel('Replication timing (min after G1)')
set(gca,'Xtick', [1:K] , 'XtickLabel' , legend_text);
xlabel('Under-replication (METpr-CDC20)')

data =  T.Raghuraman01 - abs(min(T.Raghuraman01))   ;
data = 100 * (data ./ max(data(~isinf(data))) );
data_violin = NaN(length(data) , K);
for I = 1:K
    idx = find(Classes == unq_class(I));
    data_violin(1:length(idx) , I) = data(idx);
end
subplot(2,2,2); hold on; grid on; 
violin(data_violin , 'facecolor',clrs_set,'edgecolor',[.25 .25 .25],'bw',2,'mc',[] , 'medc' , 'k');
title({'Raghuraman et. al. 2001' 'RM14-3a, \alpha-factor & cdc7-1) ; 6 reps'})
ylabel('Replication timing (% of latest)')
set(gca,'Xtick', [1:K] , 'XtickLabel' , legend_text);
xlabel('Under-replication (METpr-CDC20)')
set(gca,'ytick',0:10:100)


data = -1*T.Koren10 ;
data = data + abs(min(data)) ; 
data = 100 * (data ./ max(data) );
data_violin = NaN(length(data) , K);
for I = 1:K
    idx = find(Classes == unq_class(I));
    data_violin(1:length(idx) , I) = data(idx);
end
subplot(2,2,3); hold on; grid on; 
violin(data_violin , 'facecolor',clrs_set,'edgecolor',[.25 .25 .25],'bw',2,'mc',[] , 'medc' , 'k'); 
title({'Koren et. al. 2010' '(BY4741, sort-array S & G1) ; 4 reps'})
ylabel('Replication timing (G1/S)  (% of latest)')
set(gca,'Xtick', [1:K] , 'XtickLabel' , legend_text);
xlabel('Under-replication (METpr-CDC20)')
set(gca,'ytick',0:10:100)

%data =  T.Muller14_G2 ./ T.Muller14_S  ;
%data = 100 * (data ./ max(data(~isinf(data))) );
data =  -1*(T.Muller14-100) ; 
data_violin = NaN(length(data) , K);
for I = 1:K
    idx = find(Classes == unq_class(I));
    data_violin(1:length(idx) , I) = data(idx);
end
subplot(2,2,4); hold on; grid on; 
violin(data_violin , 'facecolor',clrs_set,'edgecolor',[.25 .25 .25],'bw',2,'mc',[] , 'medc' , 'k'); 
title({'Muller et. al. 2014' '(w303, sort-seq S & G2) ; 1 rep'} ) ; 
ylabel('Replication timing (G2/S) (% of latest)')
set(gca,'Xtick', [1:K] , 'XtickLabel' , legend_text);
xlabel('Under-replication (METpr-CDC20)')
ylim([30 105])
set(gca,'ytick',0:10:100)

print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/S7A__allRT__bin_by_underrep' , '-r600');

%% all RTs, binning by RT
fh = figure('units','centimeters','position',[5 5 20 20]); 

% Raghuraman01
subplot(2,2,1); hold on; grid on;
pctile_thresholds = [30 45 60];
legend_text = [0 pctile_thresholds] ; 
Classes = zeros(height(T),1);
data_RT = T.Raghuraman01;
for I = 1:numel(pctile_thresholds)
    Classes(data_RT >= pctile_thresholds(I)) = I ; 
end
unq_class = unique(Classes); K = length(unq_class); clrs_set = parula(K);
data = T.percent_underreplicated_cdc20*100;
data_violin = NaN(length(data) , K);
for I = 1:K
    idx = find(Classes == unq_class(I));
    data_violin(1:length(idx) , I) = data(idx);
end

violin(data_violin , 'facecolor',clrs_set,'edgecolor',[.25 .25 .25],'bw',2,'mc',[] , 'medc' , 'k');
title({'Raghuraman et. al. 2001' 'RM14-3a, \alpha-factor & cdc7-1) ; 6 reps'})
xlabel('Replication timing (min after G1)')
set(gca,'Xtick', [1:K] , 'XtickLabel' , legend_text);
ylabel('Under-replication (METpr-CDC20)');
ylim([-20 65]);

subplot(2,2,2); hold on; grid on;
pctile_thresholds = [50 60 70 90 99];
legend_text = [0 pctile_thresholds] ; 
Classes = zeros(height(T),1);
data_RT =  T.Raghuraman01 - abs(min(T.Raghuraman01))   ;
data_RT = 100 * (data_RT ./ max(data_RT(~isinf(data_RT))) );
for I = 1:numel(pctile_thresholds)
    Classes(data_RT >= pctile_thresholds(I)) = I ; 
end
unq_class = unique(Classes); K = length(unq_class); clrs_set = parula(K);
data = T.percent_underreplicated_cdc20*100;
data_violin = NaN(length(data) , K);
for I = 1:K
    idx = find(Classes == unq_class(I));
    data_violin(1:length(idx) , I) = data(idx);
end

violin(data_violin , 'facecolor',clrs_set,'edgecolor',[.25 .25 .25],'bw',2,'mc',[] , 'medc' , 'k');
title({'Raghuraman et. al. 2001' 'RM14-3a, \alpha-factor & cdc7-1) ; 6 reps'})
xlabel('Replication timing (% of latest)');
set(gca,'Xtick', [1:K] , 'XtickLabel' , legend_text);
ylabel('Under-replication (METpr-CDC20)');
ylim([-20 65]);

subplot(2,2,3); hold on; grid on;
pctile_thresholds = [50 60 70 90 99];
legend_text = [0 pctile_thresholds] ; 
Classes = zeros(height(T),1);
data_RT = -1*T.Koren10 ;
data_RT = data_RT + abs(min(data_RT)) ; 
data_RT = 100 * (data_RT ./ max(data_RT) );
for I = 1:numel(pctile_thresholds)
    Classes(data_RT >= pctile_thresholds(I)) = I ; 
end
unq_class = unique(Classes); K = length(unq_class); clrs_set = parula(K);
data = T.percent_underreplicated_cdc20*100;
data_violin = NaN(length(data) , K);
for I = 1:K
    idx = find(Classes == unq_class(I));
    data_violin(1:length(idx) , I) = data(idx);
end
violin(data_violin , 'facecolor',clrs_set,'edgecolor',[.25 .25 .25],'bw',2,'mc',[] , 'medc' , 'k');
title({'Koren et. al. 2010' '(BY4741, sort-array S & G1) ; 4 reps'})
xlabel('Replication timing (G1/S)  (% of latest)');
set(gca,'Xtick', [1:K] , 'XtickLabel' , legend_text);
ylabel('Under-replication (METpr-CDC20)');
ylim([-20 65]);


subplot(2,2,4); hold on; grid on;
pctile_thresholds = [50 60 70 90 99];
legend_text = [0 pctile_thresholds] ; 
Classes = zeros(height(T),1);
data_RT =  -1*(T.Muller14-100) ; 
for I = 1:numel(pctile_thresholds)
    Classes(data_RT >= pctile_thresholds(I)) = I ; 
end
unq_class = unique(Classes); K = length(unq_class); clrs_set = parula(K);
data = T.percent_underreplicated_cdc20*100;
data_violin = NaN(length(data) , K);
for I = 1:K
    idx = find(Classes == unq_class(I));
    data_violin(1:length(idx) , I) = data(idx);
end

violin(data_violin , 'facecolor',clrs_set,'edgecolor',[.25 .25 .25],'bw',2,'mc',[] , 'medc' , 'k');
title({'Muller et. al. 2014' '(w303, sort-seq S & G2) ; 1 rep'})
xlabel('Replication timing (G2/S) (% of latest)');
set(gca,'Xtick', [1:K] , 'XtickLabel' , legend_text);
ylabel('Under-replication (METpr-CDC20)');
ylim([-20 65]);
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/S7A__allRT__bin_by_RT' , '-r600');

%%



