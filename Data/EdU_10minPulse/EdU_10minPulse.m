%% Script to load data and make figures from the 10min EdU pulse 
%   from cells grown 30min prior in alpha-factor
DATADIR = '~/Develop/Mendoza__ReplicationEvolution/Data/EdU_10minPulse/';
%% load data
T = readtable([ DATADIR 'stat_EdU_aF_20180802_slide2.txt' ] );

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

%% boxplots

% only take rows of interest
Q = T( strcmp(T.Fluor,'EdU') & ~strcmp(T.CellPart,'cell'),:);
%Q = Q( ~strcmp(Q.Mother_daughter,'1') , :);
for I = 1:height(Q)
    if regexpcmp(Q.Label{I},'MAX_Control')
        Q.Stage{I} = 'no EdU' ;
    end
end

% alternating rows have signal & background (not a very robust way to do this analysis. 

% confirm alternating nuc/bck by area
%     some HUGE nuclei & backgrounds. WTF!? 
figure; hold on ;ecdf( Q.Area(1:2:end));ecdf( Q.Area(2:2:end));xlim([0 15]);legend({'Nuc' 'bckgrnd'},'location','best'); title(' look into this!');xlabel('ROI area (?units?)')


EdUsig  =  Q.Mean(1:2:end)  - Q.Mean(2:2:end) ;
Stages = Q.Stage(1:2:end);
MD = Q.Mother_daughter(1:2:end);
EDU = arrayfun(@(X){'+EdU'},1:numel(MD))' ; 
EDU(regexpcmp(Q.Label(1:2:end),'MAX_Control')) = {'no edu'};
mtREP = arrayfun(@(X)num2str(X),Q.mtDNARepObserved(1:2:end),'UniformOutput',false) ; 

%EdUsig(EdUsig<=0) = 0.1 ;
EdUsig(EdUsig>40) = 40 ;
%groups_cat = strcat(Stages,',',MD) ; 
groups_cat = strcat(Stages) ; 
%groups_cat = strcat(Stages  , ',' , EDU ) ; 
%groups_cat = strcat(EDU ,' ' , Stages   ) ; 

unique_groups = unique(groups_cat) ; 
fh = figure('units','centimeters','position',[5 5 20 10]);
hold on ;

clrs = get(gca,'ColorOrder');
%clrs = repmat( clrs(1:2,:) , 5 , 1);
clrs = parula(numel(unique_groups));

%bh = boxplot( EdUsig , groups_cat ,'GroupOrder',unique_groups ,'Symbol','' ); 
bh = boxplot( EdUsig , groups_cat, 'ColorGroup',clrs  ,'GroupOrder',unique_groups ,'Symbol','' ); 

[~,p] = ttest2( EdUsig( strcmp(Stages,'A') ) , EdUsig(  strcmp(Stages,'G1'))) ; 
%[~,p] = fishertest( [ sum(EdUsig( strcmp(Stages,'G1') )>15) sum(EdUsig( strcmp(Stages,'G1'))< 15) ; ...
%    sum(EdUsig( strcmp(Stages,'T') )>15) sum(EdUsig( strcmp(Stages,'T') )<15) ])

set(gca,'yscale','lin')
ylabel('EdU (nuclear - bckgrnd)')
title('Mean EdU signal Alsu stat_EdU_WT_20170731_1')

THRESH = max(EdUsig(strcmp(EDU,'no edu'))) ;

for I = 1:numel(unique_groups)
    data = EdUsig( strcmp(groups_cat,unique_groups{I}));
    x = get(bh(1,I),'XData'); x = x(1);
    txt = sprintf('%0.0f%% (%d/%d)' , ( sum(data>THRESH) / numel(data))*100 , sum(data>THRESH) , numel(data));
    text( x-0.25 , 45 , txt);
    [~,p] = fishertest( [ sum(data>THRESH)  sum(data<=THRESH) ; 0 sum(regexpcmp(Q.Label(1:2:end),'MAX_Control'))] ) ;
    text( x-0.25 , 42 , sprintf('p=%0.03f' , p) );
    plot( random('normal',x,0.05,size(data)) , data ,'ok');
    
end

line( xlim , [THRESH THRESH], 'LineStyle','--')
%% %% %% ecdfs %% %%
% %% ecdfs 
% fh = figure('units','centimeters','position',[5 5 10 10]);hold on ;
% us = unique(Stages);
% us = {'eG1' 'G1' 'M' 'A' 'T'};
% clrs = get(gca,'ColorOrder');
% for I = 1:numel(us)
%     mother_data = ( EdUsig( strcmp(MD,'M') & strcmp(Stages,us{I})));
%     daughter_data = ( EdUsig( strcmp(MD,'D') & strcmp(Stages,us{I})));
%    [f,x] = ecdf(  ( vertcat(mother_data,daughter_data) ) );
%    plot(x,f*100,'s-','LineWidth',3,'DisplayName',us{I});
% %     [f,x] = ecdf( log10(  mother_data  ) );
% %     plot(x,f*100,'s-','LineWidth',3,'Color',clrs(I,:),'DisplayName',sprintf('M %s' , us{I}));
% %     if ~isempty(daughter_data)
% %     [f,x] = ecdf( log10( daughter_data) ) ;
% %     plot(x,f*100,'o-.','LineWidth',3,'Color',clrs(I,:),'DisplayName',sprintf('D %s' , us{I}));
% %     [~,p] = ttest2(  log10(  mother_data  ) , log10( daughter_data));
% %     fprintf('%s M vs D p=%0.04f\n' , us{I} , p );
% %    end
% 
% end
% 
% xlabel('log10( EdU (mean intensity, nuc - bckgrnd) )')
% ylabel('% of cells')
% axis tight; 
% legend('location','nw')
% title('Mean')
% %set(gca,'xtick',log10([.05 .1 .2 .5 1 2 5 10 20 50 100 200]))
