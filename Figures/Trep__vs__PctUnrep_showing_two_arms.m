%% show that not all late replicating regions are underreplicated in metaphase
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');
GC = readtable('~/Google Drive/CareyLab/ExternalData/Yeast/genome_200bpwindows_GC.bed','FileType','text','ReadVariableNames',false);
GC.Properties.VariableNames = {'chr' 'start_point' 'end_point' 'GC'};
T = join(dataset2table(DS) , GC , 'Key' , {'chr' 'start_point' 'end_point'});
%% remove NaNs & Infs
X  = T.Trep_spline ; 
Y1 = T.percent_unreplicated_not_trimmed_cdc20 ;
Y2 = T.percent_unreplicated_not_trimmed_dbf2 ;

idx = ~isnan(X) & ~isinf(X) & ~isnan(Y1) & ~isinf(Y1) & ~isnan(Y2) & ~isinf(Y2) ;

T = T( idx , :);

%% bin to 1kb windows to reduce noise

T.end_point_kb_rounded = round(T.end_point ./ 1000); 
G = grpstats( T , {'end_point_kb_rounded' 'chr'}, {'mean' } ...
    , 'DataVars',{'dist_to_the_end_kb' 'percent_unreplicated_not_trimmed_cdc20' 'Trep_spline' 'GC' 'percent_unreplicated_not_trimmed_dbf2' 'dist_to_ARS'});

T.end_point_kb_rounded = [] ; 

%% ARS cherrypicking to reduce noise
G2 = G( G.mean_dist_to_the_end_kb>50,:);
figure; hold on ;
dscatter( G2.mean_dist_to_ARS , G2.mean_percent_unreplicated_not_trimmed_cdc20); grid on ;
[~,p] = ttest2( G2.mean_percent_unreplicated_not_trimmed_cdc20(G2.mean_dist_to_ARS>60)...
    , G2.mean_percent_unreplicated_not_trimmed_cdc20(G2.mean_dist_to_ARS<20)) 

%% linear model
mdl = fitglm( G( G.mean_dist_to_the_end_kb>75,:) , 'ResponseVar', ...
    'mean_percent_unreplicated_not_trimmed_dbf2' ,'PredictorVars',{'mean_GC' 'mean_dist_to_the_end_kb' 'mean_Trep_spline' 'mean_dist_to_ARS'})
%% anova
%X = [ T.percent_unreplicated_not_trimmed_dbf2 T.dist_to_the_end_kb T.dist_to_ARS T.GC];
%X = zscore(X( T.dist_to_the_end_kb > 75 ,:));

X = [ G.mean_percent_unreplicated_not_trimmed_dbf2 G.mean_dist_to_the_end_kb G.mean_dist_to_ARS G.mean_GC];
X = zscore(X( G.mean_dist_to_the_end_kb > 75 ,:));

[~,tbl,stats] = anovan( X(:,1) , X(:,2:end) ,'Continuous',1:3,'varnames',{'toEnd' 'toARS' 'GC'},'model','linear');
(cell2mat(tbl(2:end-1,2))) ./ sum(cell2mat(tbl(2:end-1,2)))
%%
X  = G.mean_Trep_spline ; 
Y = G.mean_percent_unreplicated_not_trimmed_cdc20 ;

fh = figure('units','centimeters','position',[5 5 12 9 ]);
hold on ;
clrs = get(gca,'ColorOrder');

thresh_idx = Y > 10 ;
sh = scatter(  X(thresh_idx) , Y(thresh_idx) , 10 ,...
    clrs(1,:) ,'MarkerFaceColor',clrs(1,:),'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5 );

sh = scatter(  X , Y , 10 ,...
    'k' ,'MarkerFaceColor','k','MarkerEdgeAlpha',0.2,'MarkerFaceAlpha',0.2 );


xlim([10 65])
ylim([-30 70])
xlabel('Replication timing (minutes after {\alpha}factor release)')
ylabel('% of cells in which DNA is unreplicated')
% %%
% fh = figure('units','centimeters','position',[5 5 12 9 ]);
% 
% 
% dscatter(  X(idx) , Y(idx) , 'MSIZE' , 5 , 'FILLED' , true , 'BINS' , [200 200] ,'SMOOTHING' , 5 );
% 
% xlim([10 65])
% ylim([-30 70])
% xlabel('Replication timing (minutes after {\alpha}factor release)')
% ylabel('% of cells in which DNA is unreplicated')
%
print('-dpng2','~/Downloads/not all late replicating regions are underreplicated in metaphase_v2.png','-r300');