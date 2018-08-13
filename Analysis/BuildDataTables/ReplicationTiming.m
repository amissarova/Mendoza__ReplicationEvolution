%% builds a single dataset w/rep timing from everyone's data
DATADIR = '~/Develop/Mendoza__ReplicationEvolution/Data/ExternalData/RepTiming/' ;

%% load our data
load( [DATADIR '../../DS_stat__200bp_new.mat']);

%% chr lengths
ChrLengths = readtable( [DATADIR 'genome.chrlengths'],'Delimiter','\t','FileType','text');
ChrLengths.Properties.VariableNames = {'chr_num' 'len'};

%% load .bed files made w/bedools map mean
K = readtable([ DATADIR 'Koren10_WT_200.txt'] );
K.Properties.VariableNames = {'chr_num' 'start_point' 'Koren10'} ; 
R = readtable([ DATADIR 'Raghuraman01_200.txt'] );
R.Properties.VariableNames = {'chr_num' 'start_point' 'Raghuraman01'} ; 
A = readtable([ DATADIR 'Alvino07_200.txt'] );
A.Properties.VariableNames = {'chr_num' 'start_point' 'Alvino07'} ; 
M = readtable([ DATADIR 'Muller14.txt'] );
M.Properties.VariableNames = {'chr' 'start_point' 'stop_point' 'Muller14_exponential' 'Muller14_stationary' 'Muller14_G2' 'Muller14_S'} ; 
M.chr_num = cellfun( @(X) roman2num(regexprep(X,'chr','')) , M.chr); 
M = M( ~isnan(M.chr_num),:);

%% Muller14 normalization
G = grpstats( M ,'chr_num',{@modefit @sum @mean},'DataVars',{'Muller14_exponential' 'Muller14_stationary' 'Muller14_G2' 'Muller14_S'});
M.Muller14_exponential_norm = NaN(height(M),1);
M.Muller14_stationary_norm =  NaN(height(M),1);
M.Muller14_G2_norm =  NaN(height(M),1);
M.Muller14_S_norm =  NaN(height(M),1);
for I = 1:height(M)
    idx =  G.chr_num==M.chr_num(I) ;
    M.Muller14_exponential_norm(I) = M.Muller14_exponential(I) / G.modefit_Muller14_exponential(idx) ; 
    M.Muller14_stationary_norm(I) = M.Muller14_stationary(I) / G.modefit_Muller14_stationary(idx) ; 
    M.Muller14_G2_norm(I) = M.Muller14_G2(I) / G.modefit_Muller14_G2(idx) ; 
    M.Muller14_S_norm(I) = M.Muller14_S(I) / G.modefit_Muller14_S(idx) ; 
end
%%

T = outerjoin(K,R,'Key',{'chr_num' 'start_point'},'MergeKeys',true);
T = outerjoin(T,A,'Key',{'chr_num' 'start_point'},'MergeKeys',true);
T = outerjoin(T,M,'Key',{'chr_num' 'start_point'},'MergeKeys',true);

T = innerjoin( T , dataset2table(DS(:,{'chr_num' 'start_point' 'dist_to_the_end_kb' 'Trep_spline' 'G4' 'GC'...
    'percent_underreplicated_cdc20' 'percent_underreplicated_dbf2' })) ...
    );

%%% Question to ask: what fraction of our underrep regions are
%%% late-replicating in various papers? 

% % % %% Under-rep THRESHOLDS by % of the genome
% % % pctile_thresholds = [80 90 95 99];
% % % legend_text = {'<80%' '>20%' '>10%' '>5%' '>1%'} ; 
% % % Classes = zeros(height(T),1);
% % % for I = 1:numel(pctile_thresholds)
% % %     THRESH = prctile(T.percent_underreplicated_cdc20 , pctile_thresholds(I)); 
% % %     Classes(T.percent_underreplicated_cdc20 >= THRESH) = I ; 
% % %     fprintf('%d\t%0.0f%% = %0.02f\n' , I , pctile_thresholds(I) , 100*THRESH);
% % % end
% % % 

%% Under-rep THRESHOLDS by value
pctile_thresholds = [20 25 30 35 40 45 ];
legend_text = [0 pctile_thresholds] ; 
Classes = zeros(height(T),1);
for I = 1:numel(pctile_thresholds)
    Classes(T.percent_underreplicated_cdc20*100 >= pctile_thresholds(I)) = I ; 
end
%% boxplots
fh = figure('units','centimeters','position',[5 5 20 20]); 

% Raghuraman01
subplot(2,2,1)
bh = boxplot( T.Raghuraman01 , Classes,'Symbol','');
title({'Raghuraman et. al. 2001' 'RM14-3a, \alpha-factor & cdc7-1) ; 6 reps'})
ylabel('Replication timing (min after G1)')
set(gca,'xticklabel',legend_text)
xlabel('Under-replication (METpr-CDC20)')
yval = get(bh(5,end),'YData') ;
line( xlim , [yval(1) yval(1)] ,'LineStyle','--','Color',[.5 .5 .5])


data =  T.Raghuraman01 - abs(min(T.Raghuraman01))   ;
data = 100 * (data ./ max(data(~isinf(data))) );
subplot(2,2,2)
bh = boxplot( data , Classes,'Symbol','');
title({'Raghuraman et. al. 2001' 'RM14-3a, \alpha-factor & cdc7-1) ; 6 reps'})
ylabel('Replication timing (% of latest)')
set(gca,'xticklabel',legend_text)
xlabel('Under-replication (METpr-CDC20)')
yval = get(bh(5,end),'YData') ;
line( xlim , [yval(1) yval(1)] ,'LineStyle','--','Color',[.5 .5 .5])
set(gca,'ytick',0:10:100)


data = -1*T.Koren10 ;
data = data + abs(min(data)) ; 
data = 100 * (data ./ max(data) );
subplot(2,2,3)
bh = boxplot( data , Classes,'Symbol','') ; 
title({'Koren et. al. 2010' '(BY4741, sort-array S & G1) ; 4 reps'})
ylabel('Replication timing (G1/S)  (% of latest)')
set(gca,'xticklabel',legend_text)
xlabel('Under-replication (METpr-CDC20)')
yval = get(bh(5,end),'YData') ;
line( xlim , [yval(1) yval(1)] ,'LineStyle','--','Color',[.5 .5 .5])
set(gca,'ytick',0:10:100)

data =  T.Muller14_G2 ./ T.Muller14_S  ;
data = 100 * (data ./ max(data(~isinf(data))) );
subplot(2,2,4)
bh = boxplot(data, Classes,'Symbol','') ;
title({'Muller et. al. 2014' '(w303, sort-seq S & G2) ; 1 rep'} ) ; 
ylabel('Replication timing (G2/S) (% of latest)')
set(gca,'xticklabel',legend_text)
xlabel('Under-replication (METpr-CDC20)')
ylim([40 90])
yval = get(bh(5,end),'YData') ;
line( xlim , [yval(1) yval(1)] ,'LineStyle','--','Color',[.5 .5 .5])
set(gca,'ytick',0:10:100)