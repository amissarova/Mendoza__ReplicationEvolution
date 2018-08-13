%% Script to load data and make figures from the 10min EdU pulse 
%   from cells grown 30min prior in alpha-factor
DATADIR = '~/Develop/Mendoza__ReplicationEvolution/Data/EdU_10minPulse/';
% load data
T = readtable([ DATADIR 'stat_EdU_aF_1.txt' ] );
T = T( ~isnan(T.Area) , :);
% propagage mother/daughter, stage, Cell#, and anything else that is only
% recorded for the first row of each cell
for I = 1:height(T)
    if ~isnan(T.CellNum(I))
        CellNum = T.CellNum(I);
        Stage = T.Stage{I};
        MD = T.Mother_daughter{I};
        mtDNARepObs = T.mtDNARepObserved(I);
    else
        T.CellNum(I) = CellNum ; 
        T.Stage{I} = Stage ; 
        T.Mother_daughter{I} = MD ; 
        T.mtDNARepObserved(I) = mtDNARepObs ; 
    end
end
T.Stage = categorical(T.Stage); 
T.CellPart = categorical(T.CellPart); 
T.Fluor = categorical(T.Fluor); 

T = T( T.Fluor=='EdU',:);

%% boxplots -- build dataset & set arbitrary parameters

% only take rows of interest
T = T( (T.Fluor=='EdU' & T.CellPart ~= 'cell' ),:);
for I = 1:height(T)
    if regexpcmp(T.Label{I},'MAX_Control')
        T.Stage(I) = 'no edu' ;
    end
end


idx_nuc = T.CellPart == 'nuclei'   ;
idx_bck = T.CellPart =='cytoplasm' ;

R = table();
R.EdU_Nuc  =  T.Mean(idx_nuc) ;
R.EdU_Background =  T.Mean(idx_bck) ;
R.EdU_NucMax = T.Max(idx_nuc);
R.EdU_BckMax = T.Max(idx_bck);
R.EdU_NucMedian = T.Median(idx_nuc);
R.EdU_BckMedian = T.Median(idx_bck);

R.EdUsig_diff  =  R.EdU_Nuc - R.EdU_Background ; 
R.EdUsig_ratio = R.EdU_Nuc ./ R.EdU_Background ; 
R.EdUsig = R.EdUsig_diff   ; 

R.Stage  =  T.Stage(idx_nuc);
R.MD      =  T.Mother_daughter(idx_nuc);
R.CellNum =  T.CellNum(idx_nuc);

R.EDU_FLAG=  regexpcmp(T.Label(idx_nuc),'MAX_Control');
R.Stage( R.EDU_FLAG )  = {'no edu'} ; 

MAXCUTOFF = prctile(R.EdUsig , 98) ; 
MINCUTOFF = 1 ; 

THRESH = prctile( R.EdUsig(R.EDU_FLAG) , 100 ); % for drawing line and for FE test

unique_groups = {'eG1' 'G1' 'S' 'M' 'A' 'T'  'no edu'};

%% make plot
fh = figure('units','centimeters','position',[5 5 25 10]);
hold on ;

clrs = get(gca,'ColorOrder');
clrs = parula(numel(unique_groups));

bh = boxplot( R.EdUsig , R.Stage  , 'ColorGroup',clrs  ,'GroupOrder',unique_groups ,'Symbol','' ,'notch','on'); 

ylabel('EdU (nuclear / bckgrnd)')


for I = 1:numel(unique_groups)
    data = R.EdUsig( R.Stage == unique_groups{I} );
    dataMAXCUTOFF = data; dataMAXCUTOFF(data>MAXCUTOFF) = MAXCUTOFF ; 
    dataMAXCUTOFF(data<MINCUTOFF) = MINCUTOFF ; 

    x = get(bh(1,I),'XData'); x = x(1);    
    sh = scatter( random('uniform',x-0.2,x+0.2,size(data)) , dataMAXCUTOFF ,10,'k','MarkerFaceColor', get(bh(5,I),'Color') );

    txt = sprintf('%0.0f%% (%d/%d)' , ( sum(data>THRESH) / numel(data))*100 , sum(data>THRESH) , numel(data));
    text( x-0.25 , MAXCUTOFF+0.2 , txt);
 
     [~,p] =  fishertest( [ sum(data>THRESH)  sum(data<=THRESH) ; sum(R.EdUsig(R.EDU_FLAG)>THRESH) sum(R.EdUsig(R.EDU_FLAG)<=THRESH) ] ) ;
     p = (numel(unique_groups)-1) * p ; 
%     if p<0.1
%         text( x-0.25 , MAXCUTOFF*0.8  , sprintf('F.E. p=%0.02f' , p) );
%     end
%     
     [~,pT] = ttest2( data , R.EdUsig(R.EDU_FLAG) ) ;
     pT = (numel(unique_groups)-1) * pT ; 
%     if pT<0.1
%         text( x-0.25 , MAXCUTOFF*0.7 , sprintf('T_{test}   p=%0.02f' , pT ) );
%     end
    fprintf('%s\tFE=%0.03f\tTt=%0.03f\n' , unique_groups{I} , p , pT);
end
fprintf('\n');

line( xlim , [THRESH THRESH], 'LineStyle','--','Color',[.7 .7 .7])
ylim([MINCUTOFF  MAXCUTOFF+0.3])

%%
fh = figure('units','centimeters','position',[5 5 10 10]);
gh = gscatter( R.EdU_Background,R.EdU_Nuc , R.Stage , parula(numel(unique_groups)) , '...+..s')
xlabel('EdU Background')
ylabel('EdU Nuclear')
lh = line([10 80],[10 80],'LineStyle','-','Color','k');
set(lh,'HandleVisibility','off');
lh = line([10 80],[30 100],'LineStyle','--','Color','k');
set(lh,'HandleVisibility','off');
%set(gca,'ytick',0:5:40)

%%
grps = {'G1' 'S' 'M' 'A' 'T' 'eG1'};

Y = NaN(numel(grps),2);
for I = 1:numel(grps)
    data = R.EdUsig( R.Stage == grps{I});
    medians = bootstrp( 1e4 , @median , data);
    means = bootstrp( 1e4 , @mean , data);
    Y( I , 1 ) = median(data);
    Y( I , 2 ) = std(medians) ;
    
    Y( I , 1 ) = mean(data);
    Y( I , 2 ) = prctile(means,95)-mean(data) ;    
end
fh = figure('units','centimeters','position',[5 5 10 15]);
clrs = get(gca,'ColorOrder');
hold on ;

errorbar( 0:numel(grps)+1 ,  repmat( mean(control_data),1,8) , repmat(std(means),1,8) , 'LineStyle','-','DisplayName','no EdU control','LineWidth',2)

errorbar( 0:numel(grps)+1 , [Y(end,1) Y(:,1)' Y(end,1)], [ Y(end,2) Y(:,2)' Y(end,2)]  ,'o-k','MarkerFaceColor','k','LineWidth',3,'DisplayName','data')
plot( 0:numel(grps)+1 ,  [ Y(4,1) Y(4,1) Y(2,1) Y(4,1) Y(4,1) Y(4,1) Y(4,1)  Y(4,1)] , 's-','Color',[.7 .7 .7],'DisplayName','expectation (cartoon)','LineWidth',2)
set(gca,'xtick',1:numel(grps))
set(gca,'xticklabel',grps)
ylim([6 60])
ylabel('EdU signal (nuclear - background)')
xlabel('cell cycle stage')
legend('location','ne')
xlim([0.5 6.5])

control_data = R.EdUsig( R.Stage=='no edu') ; 
means = bootstrp( 1e4 , @mean , control_data);

%%
% figure; 
% [p,t,stats] = anova1(  R.EdUsig , R.Stage , 'off');
% [c,m,h,nms] = multcompare(stats,'Display','on');

% 
% %% taking mean of nuclei for cells that haven't yet finished cytokin
% idx = strcmp(R.Stages,'A')  | strcmp(R.Stages,'T') ;
% G = grpstats( R(idx,:) , {'CellNum' 'Stages' 'EDU_FLAG'} , {'mean'} , 'DataVars' ,'EdUsig');
% G.GroupCount = [] ;
% G.Properties.VariableNames = {'CellNum' 'Stages' 'EDU_FLAG' 'EdUsig'}; 
% G.Properties.RowNames = {} ; 
% G = vertcat(G,R(~idx,G.Properties.VariableNames)); 
% 
% fh = figure('units','centimeters','position',[5 5 25 10]);
% hold on ;
% bh = boxplot( G.EdUsig , G.Stages , 'ColorGroup',clrs  ,'GroupOrder',unique_groups ,'Symbol','' ); 
% 
% ylabel('EdU (nuclear - bckgrnd)')
% 
% for I = 1:numel(unique_groups)
%     data = G.EdUsig( strcmp(G.Stages,unique_groups{I}));
%     dataMAXCUTOFF = data; dataMAXCUTOFF(data>MAXCUTOFF) = MAXCUTOFF ; 
%     x = get(bh(1,I),'XData'); x = x(1);
%     txt = sprintf('%0.0f%% (%d/%d)' , ( sum(data>THRESH) / numel(data))*100 , sum(data>THRESH) , numel(data));
%     text( x-0.25 , 44 , txt);
%     [~,p] = fishertest( [ sum(data>THRESH)  sum(data<=THRESH) ; sum(R.EdUsig(R.EDU_FLAG)>THRESH) sum(R.EdUsig(R.EDU_FLAG)<=THRESH) ] ) ;
%     p = (numel(unique_groups)-1) * p ; 
%     sh = scatter( random('normal',x,0.05,size(data)) , dataMAXCUTOFF ,10,'k','MarkerFaceColor', get(bh(5,I),'Color') );
% 
%     if p<0.1
%         text( x-0.25 , 42 , sprintf('F.E. p=%0.02f' ,p));
%     end
%     [~,p] = ttest2( data , R.EdUsig(R.EDU_FLAG) ) ;
%     p = (numel(unique_groups)-1) * p ; 
%     if p<0.1
%         text( x-0.25 , 39.5 , sprintf('T_{test}   p=%0.02f' , p) );
%     end
%     
% end
% 
% line( xlim , [THRESH THRESH], 'LineStyle','--','Color',[.7 .7 .7])
% ylim([0 49.5])
% set(gca,'ytick',0:5:40)
% title('join nuclei by mean')
% % figure; 
% % [p,t,stats] = anova1(  G.EdUsig , G.Stages , 'off');
% % [c,m,h,nms] = multcompare(stats,'Display','on');
