%% builds a single dataset w/rep timing from everyone's data
DATADIR = '~/Develop/Mendoza__ReplicationEvolution/Data/ExternalData/RepTiming/' ;
load( [DATADIR '../../DS_stat__200bp_new.mat']);
ChrLengths = readtable( [DATADIR 'genome.chrlengths'],'Delimiter','\t','FileType','text');
ChrLengths.Properties.VariableNames = {'chr_num' 'len'};

% regions to delete
%rDNA:
%chrXII, 450.6 - 491 kb

%deletion:
%chrII, 258.8 - 259.4 kb
%chrII, 265.6 - 266 kb
%chrIII, 148.8 - 151.2 kb
%chrIII, 316.6 kb
%chrIV, 2.2 kb
%chrIV, 461.8 kb
%chrVII, 289.6 kb
%chrXV, 721.8 - 722.4 kb
c = {'chrXII' 'chrII' 'chrII' 'chrIII' 'chrIII' 'chrIV' 'chrIV' 'chrVII' 'chrXV'};
s = [450.6 258.8  265.6 148.8  316.6  2.2 461.8  289.6 721.8  ] ; 
e = [ 491  259.4  266   151.2  316.6  2.2 461.8  289.6 722.4 ];
R = table(c',s',e');

for I = 1:height(R)
    to_delete = ( strcmp(DS.chr,R.Var1{I}) & ( DS.start_point>=R.Var2(I)*1000 &  DS.end_point<=R.Var3(I)*1000 ));
    DS = DS( ~to_delete ,:);
    R.RemovedN(I) = sum(to_delete);
end

R = readtable([ DATADIR 'Raghuraman01_200.txt'] );
R.Properties.VariableNames = {'chr_num' 'start_point' 'Raghuraman01'} ; 



T = innerjoin( R , dataset2table(DS(:,{'chr_num' 'start_point' 'dist_to_the_end_kb' 'Trep_spline' 'G4' 'GC'...
    'percent_underreplicated_cdc20' 'percent_underreplicated_dbf2' })) ...
    );

%% optionally, limit to 50kb from telomere

T = T(T.dist_to_the_end_kb > 50 ,:);
%% alsu's violin style
pctile_thresholds = [20 25 30 35 40 45 ];
legend_text = [0 pctile_thresholds] ; 
Classes = zeros(height(T),1);
for I = 1:numel(pctile_thresholds)
    Classes(T.percent_underreplicated_cdc20*100 >= pctile_thresholds(I)) = I ; 
end
unq_class = unique(Classes); K = length(unq_class); clrs_set = autumn(K);

%fh = figure('units','centimeters','position',[5 5 15 30]); 
%
vec_fig = [5 5 9 5];
% Raghuraman01
%subplot(5,1,1); hold on; grid on;
fh = figure('units','centimeters','position',vec_fig); hold on; grid on;
data = T.Raghuraman01;
data_violin = NaN(length(data) , K);
for I = 1:K
    idx = find(Classes == unq_class(I));
    data_violin(1:length(idx) , I) = data(idx);
end
violin(data_violin , 'facecolor',clrs_set,'edgecolor',[.25 .25 .25],'bw',2,'mc',[] , 'medc' , 'k');
%title({'Raghuraman et. al. 2001' 'RM14-3a, \alpha-factor & cdc7-1) ; 6 reps'})
%ylabel('Replication timing (min after G1)')
set(gca,'Xtick', [1:K] , 'XtickLabel' , legend_text);
xlabel('Under-replication (METpr-CDC20)')


data =  T.Raghuraman01 - abs(min(T.Raghuraman01))   ;
data = 100 * (data ./ max(data(~isinf(data))) );
data_violin = NaN(length(data) , K);
for I = 1:K
    idx = find(Classes == unq_class(I));
    data_violin(1:length(idx) , I) = data(idx);
end
%subplot(5,1,2); hold on; grid on; 
fh = figure('units','centimeters','position',vec_fig); hold on; grid on;
violin(data_violin , 'facecolor',clrs_set,'edgecolor',[.25 .25 .25],'bw',2,'mc',[] , 'medc' , 'k');
%title({'Raghuraman et. al. 2001' 'RM14-3a, \alpha-factor & cdc7-1) ; 6 reps'})
%ylabel('Replication timing (% of latest)')
set(gca,'Xtick', [1:K] , 'XtickLabel' , legend_text);
xlabel('Under-replication (METpr-CDC20)')
set(gca,'ytick',0:25:100)


vec_fig = [5 5 9 5];
% Raghuraman01
%subplot(5,1,1); hold on; grid on;
fh = figure('units','centimeters','position',vec_fig); hold on; grid on; 
pctile_thresholds = [30 45 60];
legend_text = [0 pctile_thresholds] ; 
Classes = zeros(height(T),1);
data_RT = T.Raghuraman01;
for I = 1:numel(pctile_thresholds)
    Classes(data_RT >= pctile_thresholds(I)) = I ; 
end
unq_class = unique(Classes); K = length(unq_class); clrs_set = cool(K);
data = T.percent_underreplicated_cdc20*100;
data_violin = NaN(length(data) , K);
for I = 1:K
    idx = find(Classes == unq_class(I));
    data_violin(1:length(idx) , I) = data(idx);
end

violin(data_violin , 'facecolor',clrs_set,'edgecolor',[.25 .25 .25],'bw',2,'mc',[] , 'medc' , 'k');
%title({'Raghuraman et. al. 2001' 'RM14-3a, \alpha-factor & cdc7-1) ; 6 reps'})
xlabel('Replication timing (min after G1)')
set(gca,'Xtick', [1:K] , 'XtickLabel' , legend_text);
%ylabel('Under-replication (METpr-CDC20)');
ylim([-20 65]);

%subplot(5,1,2); hold on; grid on;
fh = figure('units','centimeters','position',vec_fig); hold on; grid on; 
pctile_thresholds = [50 60 70 90 99];
legend_text = [0 pctile_thresholds] ; 
Classes = zeros(height(T),1);
data_RT =  T.Raghuraman01 - abs(min(T.Raghuraman01))   ;
data_RT = 100 * (data_RT ./ max(data_RT(~isinf(data_RT))) );
for I = 1:numel(pctile_thresholds)
    Classes(data_RT >= pctile_thresholds(I)) = I ; 
end
unq_class = unique(Classes); K = length(unq_class); clrs_set = winter(K);
data = T.percent_underreplicated_cdc20*100;
data_violin = NaN(length(data) , K);
for I = 1:K
    idx = find(Classes == unq_class(I));
    data_violin(1:length(idx) , I) = data(idx);
end

violin(data_violin , 'facecolor',clrs_set,'edgecolor',[.25 .25 .25],'bw',2,'mc',[] , 'medc' , 'k');
%title({'Raghuraman et. al. 2001' 'RM14-3a, \alpha-factor & cdc7-1) ; 6 reps'})
xlabel('Replication timing (% of latest)');
set(gca,'Xtick', [1:K] , 'XtickLabel' , legend_text);
%ylabel('Under-replication (METpr-CDC20)');
ylim([-20 65]);
