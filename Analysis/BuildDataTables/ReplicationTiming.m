%% builds a single dataset w/rep timing from everyone's data
DATADIR = '~/Develop/Mendoza__ReplicationEvolution/Data/ExternalData/RepTiming/' ;

%% load our data
load( [DATADIR '../../DS_stat__200bp_new.mat']);

%% chr lengths
ChrLengths = readtable( [DATADIR 'genome.chrlengths'],'Delimiter','\t','FileType','text');
ChrLengths.Properties.VariableNames = {'chr_num' 'len'};


%% regions to delete
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



%% Muller14 .wig files
M14tc = table();
for TimePoint = [25 30 35 40 45 50 90]
    R = readtable( [ DATADIR '/Muller14__GSE48212/' sprintf('t%d.bed',TimePoint) ] ,'FileType','text');
    R.chr_num = str2double(regexprep(R.Var1,'chr',''));
    R.M14tc_Ratio = R.Var5;
    R.start_point = R.Var2-500 ;  % move from 500 to 000
    R.TimePoint = repmat( TimePoint ,  height(R) , 1);
    R = sortrows(R,{'chr_num' 'start_point'},'ascend');
    M14tc = vertcat( M14tc , R( : , { 'TimePoint' 'chr_num' 'start_point' 'M14tc_Ratio'}));
    R.start_point = R.start_point + 200 ; % copy from 500 to 200
    M14tc = vertcat( M14tc , R( : , { 'TimePoint' 'chr_num' 'start_point' 'M14tc_Ratio'}));
    R.start_point = R.start_point + 200 ; % copy from 500 to 400
    M14tc = vertcat( M14tc , R( : , { 'TimePoint' 'chr_num' 'start_point' 'M14tc_Ratio'}));
    R.start_point = R.start_point + 200 ; % copy from 500 to 600
    M14tc = vertcat( M14tc , R( : , { 'TimePoint' 'chr_num' 'start_point' 'M14tc_Ratio'}));
    R.start_point = R.start_point + 200 ; % copy from 500 to 800
    M14tc = vertcat( M14tc , R( : , { 'TimePoint' 'chr_num' 'start_point' 'M14tc_Ratio'}));
end
M14tc = M14tc( ~isnan(M14tc.chr_num) ,:); %remove mito

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
% G = grpstats( M ,'chr_num',{@modefit @sum @mean},'DataVars',{'Muller14_exponential' 'Muller14_stationary' 'Muller14_G2' 'Muller14_S'});
% M.Muller14_exponential_norm = NaN(height(M),1);
% M.Muller14_stationary_norm =  NaN(height(M),1);
% M.Muller14_G2_norm =  NaN(height(M),1);
% M.Muller14_S_norm =  NaN(height(M),1);
% for I = 1:height(M)
%     idx =  G.chr_num==M.chr_num(I) ;
%     M.Muller14_exponential_norm(I) = M.Muller14_exponential(I) / G.modefit_Muller14_exponential(idx) ; 
%     M.Muller14_stationary_norm(I) = M.Muller14_stationary(I) / G.modefit_Muller14_stationary(idx) ; 
%     M.Muller14_G2_norm(I) = M.Muller14_G2(I) / G.modefit_Muller14_G2(idx) ; 
%     M.Muller14_S_norm(I) = M.Muller14_S(I) / G.modefit_Muller14_S(idx) ; 
% end
M.Muller14 = NaN(height(M),1);
for chrI = 1:16
    idx = M.chr_num == chrI ; 
    X = M.start_point(idx)./1000;
    Y = M.Muller14_S(idx) ./ M.Muller14_G2(idx) ; 
    Yn = M.Muller14_S(idx) ;
    Yd = M.Muller14_G2(idx);
    Y = Yn ./ Yd ; 
    [ Ys  , P ] = csaps( X , Y  , .1 , X) ; 
    Ys = (Ys - min(abs(Ys)));
    Ys = ( Ys ./ max(Ys) ) * 100 ; 
    M.Muller14(idx) = Ys ; 
end


%%

T = outerjoin(K,R,'Key',{'chr_num' 'start_point'},'MergeKeys',true);
T = outerjoin(T,A,'Key',{'chr_num' 'start_point'},'MergeKeys',true);
T = outerjoin(T,M,'Key',{'chr_num' 'start_point'},'MergeKeys',true);
T = outerjoin(T,M14tc(M14tc.TimePoint==50,:),'Key',{'chr_num' 'start_point'},'MergeKeys',true);

T = innerjoin( T , dataset2table(DS(:,{'chr_num' 'start_point' 'dist_to_the_end_kb' 'Trep_spline' 'G4' 'GC'...
    'percent_underreplicated_cdc20' 'percent_underreplicated_dbf2' })) ...
    );


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
%line( xlim , [yval(1) yval(1)] ,'LineStyle','--','Color',[.5 .5 .5])


data =  T.Raghuraman01 - abs(min(T.Raghuraman01))   ;
data = 100 * (data ./ max(data(~isinf(data))) );
subplot(2,2,2)
bh = boxplot( data , Classes,'Symbol','');
title({'Raghuraman et. al. 2001' 'RM14-3a, \alpha-factor & cdc7-1) ; 6 reps'})
ylabel('Replication timing (% of latest)')
set(gca,'xticklabel',legend_text)
xlabel('Under-replication (METpr-CDC20)')
yval = get(bh(5,end),'YData') ;
%line( xlim , [yval(1) yval(1)] ,'LineStyle','--','Color',[.5 .5 .5])
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

%data =  T.Muller14_G2 ./ T.Muller14_S  ;
%data = 100 * (data ./ max(data(~isinf(data))) );
data =  -1*(T.Muller14-100) ; 
subplot(2,2,4)
bh = boxplot(data, Classes,'Symbol','') ;
title({'Muller et. al. 2014' '(w303, sort-seq S & G2) ; 1 rep'} ) ; 
ylabel('Replication timing (G2/S) (% of latest)')
set(gca,'xticklabel',legend_text)
xlabel('Under-replication (METpr-CDC20)')
ylim([40 100])
yval = get(bh(5,end),'YData') ;
line( xlim , [yval(1) yval(1)] ,'LineStyle','--','Color',[.5 .5 .5])
set(gca,'ytick',0:10:100)

%%  Mulller14 timecourse
figure;
data = 100*(2+(-1*T.M14tc_Ratio));
bh = boxplot(data, Classes,'Symbol','') ;
title({'Muller et. al. 2014' '(w303, 50'') ; 1 rep'} ) ; 
ylabel('%Under-replication (G1/50'')')
set(gca,'xticklabel',legend_text)
xlabel('Under-replication (METpr-CDC20)')
ylim([0 100])
yval = get(bh(5,end),'YData') ;
line( xlim , [yval(1) yval(1)] ,'LineStyle','--','Color',[.5 .5 .5])
set(gca,'ytick',0:10:100)