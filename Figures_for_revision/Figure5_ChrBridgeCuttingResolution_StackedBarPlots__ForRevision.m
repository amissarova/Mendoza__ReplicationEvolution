%% cdc15-as vs top2-ts cdc15-as
%  the top2 plots differ from the DNA polymerase/BIR plots.
%    Because of the inefficient cytokinesis in the top2 double mutant, this 
%       graph plots all cells that resolve the bridge, whether or not they undergo 
%      cytokinesis during aquisition. If I remember correctly, Manuel wanted to mention 
%     this discrepancy and the percentage of top2 double mutants that undergo cytokinesis in the legend.
%  
% In short, looking at all cells that resolve the bridge, 
%    without contraction or with contraction, this graph is correct.
%
% The top2 cdc15 strain specifically is particularly bad at doing cytokinesis after inhibitor washout (~20%) compared to the other strains (60-80%). 
%        I speculate this could be due to background difference in inhibitor sensitivity, 
%       the top2 allele was crossed in. 
%

%% load data
DATADIR = '~/Downloads/' ; 
DATAFILE = [ DATADIR 'Ivanova source data.xlsx']  ; 
Fig5b = readtable( DATAFILE , 'Sheet' , 'fig5b','Range','A2:K4');
Fig5b.Var4 = []; 
Fig5b.Var8 = [] ; 
Fig5b.Var5 = [] ; 
Fig5b.Var9 = [] ; 
Fig5b.Properties.RowNames = Fig5b.Var1 ;
Fig5b.Var1 = [] ;



R = table(); 
R.experiment = regexprep( Fig5b.Properties.VariableNames , '_.*','')';
R.PercentResolved = (table2array(Fig5b(2,:)) ./ sum(table2array(Fig5b) )  * 100)' ; 
R.experiment = categorical(R.experiment) ; 

N = sum(reshape( sum(table2array(Fig5b) ) , 2 , []),2) ; 
notcut = sum( reshape( table2array(Fig5b(2,:)) , 2 , []) , 2) ; 
PctResolved = notcut ./ N * 100 ; 
%Fig5b = stack( Fig5b ,  Fig5b.Properties.VariableNames  , 'NewDataVariableName' , 'count'  , 'IndexVariableName' , 'experiment');
%Fig5b.replicateN = ([1 1 2 2 3 3 1 1 2 2 3 3 ])' ;
%Fig5b.cut = string(regexprep(Fig5b.Properties.RowNames,'_.*','')) ;
%Fig5b.experiment = regexprep(string(Fig5b.experiment),'_.*','') 
%% set constans for making the figures 
clrs = cbrewer('qual','Dark2',8) ; 
FA = 0.75 ; 

% %% B cdc15-as1 vs cdc15-as1 + top2-ts
% % Cells quantified here are only those that loose the bridge during the movie. Those cells are categorized as cut,
% %    loose the bridge when the ring has contracted, or as not cut, loose the bridge before the ring has contracted,
% %    regardless of whether or not it contracts at all during the movie.
% %
% % A suggestion would be to label the y-axis with something that describes a loss of the bridge due 
% %    to any reason. I would label the bars as "cut" and/or "not cut" but you could of course replace 
% %     "not cut" with "resolved".
% 
% B_resolve = [ 82.0779220779221 91.42857142857144 ] ; 
% 
% N_B_resolve = [ 49 55 ] ; 
% V_B_resolve = { zeros(N_B_resolve(1),1) zeros(N_B_resolve(2),1) } ; 
% V_B_resolve{1}(1:round(N_B_resolve(1)*(B_resolve(1)/100))) = 1; 
% V_B_resolve{2}(1:round(N_B_resolve(2)*(B_resolve(2)/100))) = 1; 
% 
% m1 = 100*bootstrp( 1e5 , @mean , V_B_resolve{1} );
% m2 = 100*bootstrp( 1e5 , @mean , V_B_resolve{2} );
% 
% fh = figure('units','centimeters','position',[5 5 4 8]);
% hold on ;
% bar( [ 100 100] , 'FaceColor' , [.95 .95 .95]);
% bar(1 , B_resolve(1) , 'FaceColor' ,clrs(1,:) ,'FaceAlpha',FA);
% bar(2 , B_resolexitve(2) , 'FaceColor' ,clrs(7,:) ,'FaceAlpha',FA);
% errorbar( 1 , mean(m1) , std(m1) , '.' ,'LineWidth',2 ,'Color',clrs(1,:)); 
% errorbar( 2 , mean(m2) , std(m2) , '.' ,'LineWidth',2 ,'Color',clrs(6,:)); 
% set(gca,'xtick',[])
% xlim([0.5 2.5])
% ylabel('% of cells with a bridge that goes away')

% plot showing results from all replicates for final version
fh = figure('units','centimeters','position',[5 5 4 8]);
hold on ;
bar( [ 100 100] , 'FaceColor' , [.95 .95 .95]);
bar(1 , PctResolved(2) , 'FaceColor' ,clrs(1,:) ,'FaceAlpha',FA);
bar(2 , PctResolved(1) , 'FaceColor' ,clrs(7,:) ,'FaceAlpha',FA);
m1 = (R.PercentResolved(R.experiment=='cdc15')) ; 
m2 = (R.PercentResolved(R.experiment=='cdc15Top2')) ; 
errorbar( 1 , mean(m1) , std(m1) , '.' ,'LineWidth',0.5 ,'Color','k'); 
errorbar( 2 , mean(m2) , std(m2) , '.' ,'LineWidth',0.5 ,'Color','k'); 
plot( [1 1 1]' , R.PercentResolved(R.experiment=='cdc15') ,'ok','LineWidth',1) ;
plot( [2 2 2]' , R.PercentResolved(R.experiment=='cdc15Top2') ,'ok','LineWidth',1  ) ;
set(gca,'xtick',[])
xlim([0.5 2.5])
ylabel('% of cells with a bridge that goes away')

%%
A_resolve = [12.468 9.351 ] 
hold on ;
bar( [ 100 100] , 'FaceColor' , [.95 .95 .95]);
bar(1 , B_resolve(1) , 'FaceColor' ,clrs(3,:) ,'FaceAlpha',0.5);
bar(2 , B_resolve(2) , 'FaceColor' ,clrs(1,:) ,'FaceAlpha',0.5);
set(gca,'xtick',[])
xlim([0.5 2.5])
ylabel('% of cells with a visible bridge')

x = [ sum(V_B_resolve{1}==0)  sum(V_B_resolve{1}==1) ; sum(V_B_resolve{2}==0)  sum(V_B_resolve{2}==1) ] ; 
[~,p,~] = fishertest(x)

%% cdc15-as1 vs cohesin-ts + cdc15-as1

A_resolved = [10.309 12.680  ] ; 

N_A = [ 114 206 ] ; 
V_A = { zeros(N_A(1),1) zeros(N_A(2),1) } ; 
V_A{1}(1:round(N_A(1)*(A_resolved(1)/100))) = 1; 
V_A{2}(1:round(N_A(2)*(A_resolved(2)/100))) = 1; 

fh = figure('units','centimeters','position',[5 5 5 6]);
hold on 
bar( [100 100] , 'FaceColor',[.95 .95 .95] );
bar( 1 , A_resolved(1) , 'FaceColor', clrs(1,:), 'FaceAlpha',FA);
bar( 2 , A_resolved(2) , 'FaceColor', clrs(6,:), 'FaceAlpha',FA);


ylabel('% of cells with a visible bridge')
set(gca,'xtick',[])
xlim([0.25 2.75])
ylim([0 40])
set(gca,'ytick',[0 10 20  40])
set(gca,'yticklabel',[0 10 20  100])


%x = [ sum(V_RFA{1}==0)  sum(V_RFA{1}==1) ; sum(V_RFA{2}==0)  sum(V_RFA{2}==1) ] ; 
%[~,p,~] = fishertest(x)

%% pol32 & BIR ; pol3-ts & pol2-ts
% Michael : In the polymerase/BIR plots I have only included cells that undergo cytokinesis 
%          and ignored those that resolve bridges in absence of cytokinesis. 
C_cut = [ 17.5824 18.813  36.7473 63.29670  49.934] ;
fh = figure('units','centimeters','position',[5 5 5 8]);
hold on 
bar( [100 100 100 100 100] , 'FaceColor',[.95 .95 .95] );
bar( 1 , 100-C_cut(1) , 'FaceColor', clrs(1,:), 'FaceAlpha',FA);
bar( 2 , 100-C_cut(2) , 'FaceColor', clrs(2,:),  'FaceAlpha',FA);
bar( 3 , 100-C_cut(3) , 'FaceColor', clrs(3,:),  'FaceAlpha',FA);
bar( 4 , 100-C_cut(4) , 'FaceColor', clrs(4,:),  'FaceAlpha',FA);
bar( 5 , 100-C_cut(5) , 'FaceColor', clrs(5,:),  'FaceAlpha',FA);

ylabel('% of cells w/bridges that do cytokinesis')
set(gca,'xtick',[])
xlim([0.5 5.5])


%% D :  Rfa2-foci
D_RFA = [ 32.325 74.76298 ] ;
N_RFA = [ 114 206 ] ; 
V_RFA = { zeros(N_RFA(1),1) zeros(N_RFA(2),1) } ; 
V_RFA{1}(1:round(N_RFA(1)*(D_RFA(1)/100))) = 1; 
V_RFA{2}(1:round(N_RFA(2)*(D_RFA(2)/100))) = 1; 


m1 = 100*bootstrp( 1e5 , @mean , V_RFA{1} );
m2 = 100*bootstrp( 1e5 , @mean , V_RFA{2} );

fh = figure('units','centimeters','position',[5 5 4 6]);
hold on 
bar( 2 , D_RFA(1) , 'FaceColor', 'k', 'FaceAlpha',FA);
errorbar( 2 , mean(m1) , std(m1) , '.k' ,'LineWidth',2 ); 
bar( 1 , D_RFA(2) , 'FaceColor', clrs(1,:),  'FaceAlpha',FA);
errorbar( 1 , mean(m2) , std(m2) , '.' ,'LineWidth',2 ,'Color',clrs(1,:)); 
set(gca,'xtick',[])
ylabel(' ')
xlim([0.5 2.5])
ylim([0 100])

fh = figure('units','centimeters','position',[5 5 6 7]);
ylabel('% of anaphase cells with RPA foci')

x = [ sum(V_RFA{1}==0)  sum(V_RFA{1}==1) ; sum(V_RFA{2}==0)  sum(V_RFA{2}==1) ] ; 
[~,p,~] = fishertest(x)

%% D from excel showing individual data points
Fig5d = readtable( DATAFILE , 'Sheet' , 'fig5d','Range','A1:C12','ReadVariableNames',true);
Fig5d = Fig5d( ~isnan(Fig5d.RPAFocus) , :);
Fig5d.repN = [1 1 1 2 2 2 3 3 3 ]' ;
Fig5d.Pct = Fig5d.RPAFocus ./ Fig5d.total * 100  ; 
Fig5d.mov1 = categorical( regexprep( Fig5d.mov1 , ' c.' ,'') ) ;

Fig5dG = grpstats( Fig5d , 'mov1' ,{'sum' 'mean' 'sem'}, 'DataVars' , {'RPAFocus' 'total' 'Pct'});

fh = figure('units','centimeters','position',[5 5 4 6]);
hold on 
bar( 1 , 100*Fig5dG.sum_RPAFocus(2)/Fig5dG.sum_total(2) , 'FaceColor', clrs(1,:), 'FaceAlpha',FA);
errorbar( 1 , Fig5dG.mean_Pct(2) , Fig5dG.sem_Pct(2) , '.k'  ); 

bar( 2 ,100*Fig5dG.sum_RPAFocus(1)/Fig5dG.sum_total(1)  , 'FaceColor', 'k',  'FaceAlpha',FA);
errorbar( 2 , Fig5dG.mean_Pct(1) , Fig5dG.sem_Pct(1) , '.k'  ); 

plot( 1 , Fig5d.Pct( Fig5d.mov1 == 'cdc15') , 'ok')
plot( 2 , Fig5d.Pct( Fig5d.mov1 == 'WT') , 'ok')

set(gca,'xtick',[])
ylabel(' ')
xlim([0.5 2.5])
ylim([0 100])

[~,p] = ttest2( Fig5d.Pct( Fig5d.mov1 == 'cdc15')  , Fig5d.Pct( Fig5d.mov1 == 'WT') )
fprintf( 'ttest across replicates, p = %0.05f\n' , p)

[~,p]=fishertest( [ Fig5dG.sum_RPAFocus(1) Fig5dG.sum_total(1)-Fig5dG.sum_RPAFocus(1) ...
    ; Fig5dG.sum_RPAFocus(2) Fig5dG.sum_total(2)-Fig5dG.sum_RPAFocus(2)]) ;
fprintf( 'Fisher''s exact test , p = %0.05e\n' , p)


%%

A_resolved = [10.309 12.680  ] ; 

fh = figure('units','centimeters','position',[5 5 4 7]);
hold on 
bar( [100 100] , 'FaceColor',[.95 .95 .95] );
bar( 1 , A_resolved(1) , 'FaceColor', 'k', 'FaceAlpha',0.5);
bar( 2 , A_resolved(2) , 'FaceAlpha',0.5);
ylabel('% of cells with bridges')
set(gca,'xtick',[])
xlim([0.5 2.5])
%%
C_cut = [ 17.5824 18.813  36.7473 63.29670  49.934] ;
fh = figure('units','centimeters','position',[5 5 4 7]);
clrs = get(gca,'ColorOrder');
hold on 
bar( [100 100 100 100 100] , 'FaceColor',[.95 .95 .95] );
bar( 1 , 100-C_cut(1) , 'FaceColor', clrs(1,:), 'FaceAlpha',0.5);
bar( 2 , 100-C_cut(2) , 'FaceColor', clrs(2,:),  'FaceAlpha',0.5);
bar( 3 , 100-C_cut(3) , 'FaceColor', clrs(3,:),  'FaceAlpha',0.5);
bar( 4 , 100-C_cut(4) , 'FaceColor', clrs(4,:),  'FaceAlpha',0.5);
bar( 5 , 100-C_cut(5) , 'FaceColor', clrs(5,:),  'FaceAlpha',0.5);

ylabel('% of cells with bridges')
set(gca,'xtick',[])
xlim([0.5 5.5])