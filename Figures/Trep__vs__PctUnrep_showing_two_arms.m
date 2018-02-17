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
G = sortrows( G , 'mean_percent_unreplicated_not_trimmed_cdc20' , 'descend');
%% ARS cherrypicking to reduce noise
G2 = G( G.mean_dist_to_the_end_kb>75,:);
%figure; hold on ;
%dscatter( G2.mean_dist_to_ARS , G2.mean_percent_unreplicated_not_trimmed_cdc20); grid on ;
%[~,p] = ttest2( G2.mean_percent_unreplicated_not_trimmed_cdc20(G2.mean_dist_to_ARS>35)...
%    , G2.mean_percent_unreplicated_not_trimmed_cdc20(G2.mean_dist_to_ARS<20)) 
fh = figure('units','centimeters','position',[5 5 8 8 ]);
hold on ;
GT =  G2.mean_percent_unreplicated_not_trimmed_cdc20( G2.mean_dist_to_ARS>35) ; 
LT =  G2.mean_percent_unreplicated_not_trimmed_cdc20( G2.mean_dist_to_ARS<35) ; 
[f,x]=ksdensity(GT,-25:0.1:30);plot(x,100*f,'-','LineWidth',3,'DisplayName','> 35kb');
[f,x]=ksdensity(LT,-25:0.1:30);plot(x,100*f,'-','LineWidth',3,'DisplayName','< 35kb');
[~,p]=ttest2( GT,LT);
title( sprintf('1kb win, >75kb from end, p=%0.03f' , p))
xlim([-25 30])
legend('location','nw')
set(gca,'ytick',[0 max(ylim)])
xlabel('% of cells not replicated')
ylabel('% of 1kb loci')

THRESH = 20  ;
data = 100*[mean(GT>THRESH) mean(LT>THRESH)];
fh = figure('units','centimeters','position',[5 5 2 8 ]);
bar( diag(data),'stacked');
%ylim([0 11])
set(gca,'ytick',[0 max(ylim)])
ylabel('% of loci > 20')
set(gca,'xtick',[])
xlim([0.5 2.5])

%%
fh = figure('units','centimeters','position',[5 5 15 8 ]);
dscatter( G2.mean_dist_to_ARS ,G2.mean_percent_unreplicated_not_trimmed_cdc20);
grid on ;
ylim([-20 40])
ylabel('% of cells not replicated')
xlabel('dist to ARS')
title('>75kb from end, 1kb windows')
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

%% tree model feature importance
% [trainedModel, validationRMSE , partitionedModel , r2 ] = trainRegressionModel(trainingData , predictorNames )
r2 = NaN( 3 , 4);
T.end_point_kb_rounded = round(T.end_point ./ 1000); 
G = grpstats( T , {'end_point_kb_rounded' 'chr'}, {'mean' } ...
    , 'DataVars',{ 'dist_to_the_end_kb' 'percent_unreplicated_not_trimmed_cdc20' 'Trep_spline' 'GC' 'percent_unreplicated_not_trimmed_dbf2' 'dist_to_ARS'});

T.end_point_kb_rounded = [] ; 

idx = G.mean_dist_to_the_end_kb >= 0 ; 
predictorNames = {'mean_dist_to_the_end_kb', 'mean_Trep_spline' , 'mean_dist_to_ARS' , 'mean_GC'};
for I = 1:nrows(r2)
[~,~,~,r2(I,1)] = trainRegressionModel(G(idx,:) , predictorNames) ;
[~,~,~,r2(I,2)] = trainRegressionModel(G(idx,:) , predictorNames([2 3 4])) ; 
[~,~,~,r2(I,3)] = trainRegressionModel(G(idx,:) , predictorNames([1 3 4])) ;
[~,~,~,r2(I,4)] = trainRegressionModel(G(idx,:) , predictorNames([1 2 4])) ;
[~,~,~,r2(I,5)] = trainRegressionModel(G(idx,:) , predictorNames([1 2 3])) ;
[~,~,~,r2(I,6)] = trainRegressionModel(G(idx,:) , predictorNames([1])) ;
[~,~,~,r2(I,7)] = trainRegressionModel(G(idx,:) , predictorNames([2])) ;
[~,~,~,r2(I,8)] = trainRegressionModel(G(idx,:) , predictorNames([3])) ;
[~,~,~,r2(I,9)] = trainRegressionModel(G(idx,:) , predictorNames([4])) ;
[~,~,~,r2(I,9)] = trainRegressionModel(G(idx,:) , predictorNames([1 2])) ;

end

%%
data = mean(r2(:,1)) - mean(r2(:,2:end))  ;
fh = figure('units','centimeters','position',[5 5 25 5 ]);
bar(  data )
set(gca,'xticklabel',regexprep( predictorNames , '_' ,' '))
set(gca,'yscale','log')