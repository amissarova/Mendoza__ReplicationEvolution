%% 2019.10.17 -- new data on nuclear segregation and chr bridges
% Figure1_nuclear_segragation_and_bridges_new_data.m


%% load data
PROJECTDIR = '~/Develop/Mendoza__ReplicationEvolution/' ;
DATADIR = '~/Develop/Mendoza__ReplicationEvolution/figure1_data_analysis/2019.10.17 -- new data on nuclear segregation and chr bridges/' ; 
DATA_FILE_1 = [ PROJECTDIR 'Figures_for_revision/Fig2 old data/division pols123.txt' ] ;
DATA_FILE_2 = [ PROJECTDIR 'Figures_for_revision/Fig2 old data/division pol2 noco.txt' ] ;
DATA_FILE_3 = [ PROJECTDIR 'Figures_for_revision/Fig2 old data/bridges pols123.txt' ] ;
DATA_FILE_4 = [ PROJECTDIR 'Figures_for_revision/Fig2 old data/bridges pol2 noco.txt' ] ;

FIGUREOUTDIR = [ PROJECTDIR '/Figures_for_revision/' ] ;  system(['mkdir -p ' FIGUREOUTDIR]);

T1 = readtable( DATA_FILE_1 ) ; 
T1.Properties.VariableNames = {'WT' 'pol1' 'pol3' 'pol2' } ; 

T2 = readtable( DATA_FILE_2 ) ; 
T2.Properties.VariableNames = {'WT' 'pol2'} ; 

T3 = readtable( DATA_FILE_3 ) ; 
T3.Properties.VariableNames = {'WT' 'pol1' 'pol3' 'pol2' } ; 

T4 = readtable( DATA_FILE_4 ) ; 
T4.Properties.VariableNames = {'WT' 'pol2' } ; 

S3 = stack(T3,1:4) ; S3.Properties.VariableNames = {'genotype' 'time'} ;
S4 = stack(T4,1:2) ; 

%T1 = [T1 ;  arrayfun(@(X){X},[300 max(table2array(T1(:,2:end)))] )  ];
%T2 = T2(1:end-1,:);
%T2 = [T2 ;  arrayfun(@(X){X},[300 max(table2array(T2(:,2:end)))] )  ];
%% define colors for the lines
clrs1 = spring(12); clrs2 = winter(12); clrs3 = parula(12); clrs4 = summer(12); clrs5 = hot(12); clrs6 = lines(6);
clrs_to_plot = cell(4,1);
lw = 2.5 ; 

clrs_to_plot{1} = [.4 .4 .4; clrs6(4,:)];
clrs_to_plot{2} = [.4 .4 .4; clrs4(4,:); clrs2(7,:)];
clrs_to_plot{3} = [.4 .4 .4; clrs1(7,:)];
clrs_to_plot{4} = [.4 .4 .4; clrs2(7,:)];

figure; 
clrs = get(gca,'ColorOrder');
clrs = [ clrs(1:2,:) ; clrs(4,:)];
close;

clrs = cbrewer('qual','Dark2',8) ; 
lw = 2.5 ;


%% Main figure; 

%% nuclear division timecourse figures

% pol1 , pol2 , pol3 + CDC20
T1.pol3(T1.pol3==300)=1e9 ; T1.pol2(T1.pol2==300)=1e9 ; T1.pol1(T1.pol1==300)=1e9 ;  T1.WT(T1.WT==300)=1e9 ;
x1 = [sum(T1.WT<300) sum(T1.WT>=300) ; sum(T1.pol1<300) sum(T1.pol1>=300)  ] ; 
x2 = [sum(T1.WT<300) sum(T1.WT>=300) ; sum(T1.pol2<300) sum(T1.pol2>=300)  ] ; 
x3 = [sum(T1.WT<300) sum(T1.WT>=300) ; sum(T1.pol3<300) sum(T1.pol3>=300)  ] ; 

[~,p1,~] = fishertest(x1)
[~,p2,~] = fishertest(x2)
[~,p3,~] = fishertest(x3)

[mwu_p1,~,~] = ranksum(T1.WT , T1.pol1 )
[mwu_p2,~,~] = ranksum(T1.WT , T1.pol2 )
[mwu_p3,~,~] = ranksum(T1.WT , T1.pol3 )

fprintf('pol1 %0.09f\n' , p1);
fprintf('pol2 %0.09f\n' , p2);
fprintf('pol3 %0.09f\n' , p3);

figure('units','centimeters','position',[5 5 6 5]); 
hold on; 
grid on;
set(gca,'ColorOrder', clrs );
[f,x]=ecdf(T1.WT);plot( x , 100*f , 'LineWidth',lw,'Color','k','DisplayName','WT')
[f,x]=ecdf(T1.pol1);plot( x , 100*f , 'LineWidth',lw,'DisplayName','{\itpol1}','Color',clrs(1,:))
[f,x]=ecdf(T1.pol2);plot( x , 100*f , 'LineWidth',lw,'DisplayName','{\itpol2}','Color',clrs(2,:))
[f,x]=ecdf(T1.pol3);plot( x , 100*f , 'LineWidth',lw,'DisplayName','{\itpol3}','Color',clrs(3,:))

xlim([0 300]) 
legend('location','se')

% pol2 + nocodazole
T2.pol2(T2.pol2>=200)=1e9 ;  T2.WT(T2.WT>=200)=1e9 ;
x2 = [sum(T2.WT<200) sum(T2.WT>=200) ; sum(T2.pol2<200) sum(T2.pol2>=200)  ] ; 

[~,p2,~] = fishertest(x2)
[mwu_p2,~,~] = ranksum(T2.WT , T2.pol2 )

fprintf('pol2 + nocodazole %0.09f\n' , p2);

figure('units','centimeters','position',[5 5 6 5]); 
hold on; 
grid on;
set(gca,'ColorOrder', clrs );
[f,x]=ecdf(T2.WT);plot( x , 100*f , 'LineWidth',lw,'Color','k','DisplayName','WT')
[f,x]=ecdf(T2.pol2);plot( x , 100*f , 'LineWidth',lw,'DisplayName','{\itpol2}' ,'Color',clrs(4,:))

xlim([0 151]) 
legend('location','se')
%% chromosome bridge durations

%% boxplots

% pol2 + nocodazole
figure('units','centimeters','position',[5 5 4 6]); 
clrsTHIS = [0 0 0 ; clrs(4,:)];
bh = boxplot(S4.WT_pol2 , S4.WT_pol2_Indicator,'symbol','','notch','on','Color',[.2 .2 .2]);
set(bh,'LineWidth',1.3);
h = findobj(gca,'Tag','Box');
ylim([0 20])

K = numel(h);
for j = 1:K
        patch(get(h(K-j+1),'XData'),get(h(K-j+1),'YData'), clrsTHIS(j,:) ,'FaceAlpha',.8 );
end    
set(gca,'xticklabel',[])

[~,p]=ttest2(T4.WT,T4.pol2)
xlim([0.7 2.3])


%% pol1,2,3 + CDC20
figure('units','centimeters','position',[5 5 5 6]); 
clrsTHIS = [0 0 0 ; clrs(1:3,:)];
bh = boxplot(S3.time , S3.genotype,'symbol','','notch','off','Color',[.2 .2 .2],'GroupOrder',{'WT' 'pol1' 'pol2' 'pol3'});
set(bh,'LineWidth',1.3);
h = findobj(gca,'Tag','Box');
ylim([0 60])

K = numel(h);
for j = 1:K
        patch(get(h(K-j+1),'XData'),get(h(K-j+1),'YData'), clrsTHIS(j,:) ,'FaceAlpha',.8 );
end    
set(gca,'xticklabel',[])

[~,p1]=ttest2(T3.WT,T3.pol1) ; fprintf('pol1 + CDCD20 %0.09f\n' , p1);
[~,p2]=ttest2(T3.WT,T3.pol2) ; fprintf('pol2 + CDCD20 %0.09f\n' , p2);
[~,p3]=ttest2(T3.WT,T3.pol3) ; fprintf('pol3 + CDCD20 %0.09f\n' , p3);
