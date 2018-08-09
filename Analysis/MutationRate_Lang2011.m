%% load data
PROJECTDIR = '~/Develop/Mendoza__ReplicationEvolution/' ; 
T = readtable([ PROJECTDIR 'Data/ExternalData/Lang11.xlsx'] );
load([ PROJECTDIR '/Data/DS_stat__features_new.mat'] );
T = innerjoin( T , dataset2table(DS) , 'LeftKey','Locus','RightKey','ORF') ; 
clear 'DS' ;


%% build predictive model on all the data
Q = table();
Q.URA3_10_8 = T.URA3_10_8 ;
Q.x_RepTime = zscore(T.x_RepTime);
Q.x_GC = zscore(T.x_GC);
Q.percent_underreplicated_cdc20 = zscore(T.percent_underreplicated_cdc20);
Q.percent_underreplicated_cdc20_DM = zscore(T.percent_underreplicated_cdc20_DM);
Q.dist_to_ARS = zscore(T.dist_to_ARS);
Q.x_RecRate = zscore(T.x_RecRate);
%Q.dist_to_the_end_kb = zscore(T.dist_to_the_end_kb);
Q.dist_to_the_end_kb_log = zscore(log10(T.dist_to_the_end_kb));

GLM = stepwiseglm( Q , 'linear', 'ResponseVar' ,'URA3_10_8' ) 
GLM_rep  = fitglm( Q , 'linear', 'ResponseVar' ,'URA3_10_8', 'PredictorVars', {'x_RepTime'} ) ;

%% plot
figure;
hold on ;
plot(Q.URA3_10_8,Q.x_RepTime,'o','DisplayName','RepTiming')
plot(Q.URA3_10_8,Q.percent_underreplicated_cdc20,'o','DisplayName','CDC20')
legend('location','best')
xlabel('URA3 mutation rate')
ylabel('Rep timing or CDC20 under-rep')

%% cross-val
mat = NaN(200,4);
for I = 1:nrows(mat)
    idx = false(height(Q),1);
    idx(randsample( height(Q) , round(height(Q)*0.8) )) = true ; % cross val indexes
    mdl = fitglm( Q(idx,:) , 'linear', 'ResponseVar' ,'URA3_10_8' , 'PredictorVars', GLM.PredictorNames ); %
    mat(I,1) = corr(Q.URA3_10_8(idx) , predict(mdl,Q(idx,2:end)))^2 ; 
    mat(I,2) = corr(Q.URA3_10_8(~idx) , predict(mdl,Q(~idx,2:end)))^2 ; 
    
    mdl = fitglm( Q(idx,:) , 'linear', 'ResponseVar' ,'URA3_10_8' , 'PredictorVars', {'x_RepTime'} );
    mat(I,3) = corr(Q.URA3_10_8(idx) , predict(mdl,Q(idx,2:end)))^2 ;     
    mat(I,4) = corr(Q.URA3_10_8(~idx) , predict(mdl,Q(~idx,2:end)))^2 ; 

%    mdl = fitglm( Q(idx,:) , 'linear', 'ResponseVar' ,'URA3_10_8' , 'PredictorVars', {'percent_underreplicated_cdc20'} );
%    mat(I,5) = corr(Q.URA3_10_8(idx) , predict(mdl,Q(idx,2:end)))^2 ;     
%    mat(I,6) = corr(Q.URA3_10_8(~idx) , predict(mdl,Q(~idx,2:end)))^2 ; 
end
[~,p] = ttest2( mat(:,2) ,mat(:,4));
%%
fh = figure('units','centimeters','position',[5 5 8 9]);
boxplot(mat,'notch','on' ,'Positions',[0.9 1.1 1.4 1.6],'Symbol','')
set(gca,'xtick',[1 1.5])
set(gca,'xticklabel',{'all params' 'rep timing only' 'under-rep only' })
grid on ;
%title( sprintf('p = %0.03f' , p))
ylabel('R^2 predict URA3 mutation rate (Lang et al. 2011)')
clrs = get(gca,'ColorOrder');
h = findobj(gca,'Tag','Box');
for j=1:2:length(h)
patch( get(h(j),'XData') , get(h(j),'YData') , clrs(2,:),'FaceAlpha',0.75);
end
for j=2:2:length(h)
patch( get(h(j),'XData') , get(h(j),'YData') , clrs(1,:),'FaceAlpha',0.75);
end