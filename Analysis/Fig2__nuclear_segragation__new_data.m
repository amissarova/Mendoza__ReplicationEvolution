%% 2019.10.17 -- new data on nuclear segregation and chr bridges
% Figure1_nuclear_segragation_and_bridges_new_data.m


%% load data
PROJECTDIR = '~/CareyLab/Projects/2017__MendozaReplication/' ;
DATADIR = '~/Develop/Mendoza__ReplicationEvolution/2019.10.17 -- new data on nuclear segregation and chr bridges/' ; 
DATA_FILE_1 = [ DATADIR 'Histogram of nuclear division pooled -- psf2.txt' ] ;
DATA_FILE_2 = [ DATADIR 'Histogram of nuclear division pooled MD -- HU.txt' ] ;

FIGUREOUTDIR = [ PROJECTDIR '/Figures_for_revision/' ] ;  system(['mkdir -p ' FIGUREOUTDIR]);

T1 = readtable( DATA_FILE_1 ) ; 
T1.Properties.VariableNames = {'Hours' 'WT' 'psf2' 'rad9' 'psf2rad9'} ; 

T2 = readtable( DATA_FILE_2 ) ; 
T2.Properties.VariableNames = {'Hours' 'WT' 'HU' 'rad9' 'rad9HU'} ; 


T1 = [T1 ;  arrayfun(@(X){X},[300 max(table2array(T1(:,2:end)))] )  ];
T2 = T2(1:end-1,:);
T2 = [T2 ;  arrayfun(@(X){X},[300 max(table2array(T2(:,2:end)))] )  ];
%% define colors for the lines
clrs1 = spring(12); clrs2 = winter(12); clrs3 = parula(12); clrs4 = summer(12); clrs5 = hot(12); clrs6 = lines(6);
clrs_to_plot = cell(4,1);
lw = 2.5 ; 

clrs_to_plot{1} = [.4 .4 .4; clrs6(4,:)];
clrs_to_plot{2} = [.4 .4 .4; clrs4(4,:); clrs2(7,:)];
clrs_to_plot{3} = [.4 .4 .4; clrs1(7,:)];
clrs_to_plot{4} = [.4 .4 .4; clrs2(7,:)];
clrs = get(gca,'ColorOrder');

%% Main figure; 

%% nuclear division timecourse figures

% psf2
x1 = [ 65 1 ; 68 4 ] ; % September 4th replicate
x2 = [ 84 0 ; 29 4 ] ; % July 16th  replicate
x = x1+x2 ;
[h,p,stats] = fishertest(x);

figure('units','centimeters','position',[5 5 6 5]); 
hold on; 
grid on;
plot(T1.Hours(~isnan(T1.WT)), T1.WT(~isnan(T1.WT)) , 'Color', 'k','LineWidth',2.5 ,'DisplayName','WT');
plot(T1.Hours(~isnan(T1.psf2)) , T1.psf2(~isnan(T1.psf2))  , 'Color', clrs(1,:) ,'LineWidth',2.5 ,'DisplayName','{\itpsf2}');
xlim([0 150]);
set(gca,'ytick',0:20:100)
set(gca,'xtick',0:60:200)
legend('location','se')
fprintf('FE test p=%0.05f\n' , p ); 
print( '-dpsc2' , [FIGUREOUTDIR 'psf2.eps'] ) ; 
close; 

% HU
x1 = [ 77 0 ; 45 10 ] ; % June 25 replicate
x2 = [ 63 0 ; 47 22 ] ; % June 13 replicate
x = x1+x2 ;
[h,p,stats] = fishertest(x);

figure('units','centimeters','position',[5 5 6 5]); 
hold on; 
grid on;
plot(T2.Hours(~isnan(T2.WT)), T2.WT(~isnan(T2.WT)) , 'Color', 'k','LineWidth',2.5 ,'DisplayName','WT');
plot(T2.Hours(~isnan(T2.HU)) , T2.HU(~isnan(T2.HU))  , 'Color', clrs(2,:),'LineWidth',2.5 ,'DisplayName','HU');
xlim([0 150]);
set(gca,'ytick',0:20:100)
set(gca,'xtick',0:60:200)
legend('location','se')
fprintf('FE test p=%0.05e\n' , p ); 
print( '-dpsc2' , [FIGUREOUTDIR 'HU.eps'] ) ; 
close; 


% rad9 + HU
x1 = [ 57 1 ; 81 0 ] ; % June 13 replicate
x2 = [ 102 0 ; 79 10 ] ; % June 25 replicate
x = x1+x2 ;
[h,p,stats] = fishertest(x);
fprintf('FE test p=%0.05f\n' , p ); 

figure('units','centimeters','position',[5 5 6 5]); 
hold on; 
grid on;
plot(T2.Hours(~isnan(T2.rad9)), T2.WT(~isnan(T2.rad9)) , 'Color', [.4 .4 .4],'LineWidth',2.5 ,'DisplayName','{\itrad9}');
plot(T2.Hours(~isnan(T2.rad9HU)) , T2.rad9HU(~isnan(T2.rad9HU))  , 'Color', clrs(3,:),'LineWidth',2.5 ,'DisplayName','{\itrad9}+HU');
xlim([0 150]);
set(gca,'ytick',0:20:100)
set(gca,'xtick',0:60:200)
legend('location','se')
print( '-dpsc2' , [FIGUREOUTDIR 'rad9HU.eps'] ) ; 
close; 

% psf2+rad9
x =  [ 174 2  ; 102 3  ]; % rad9 ; psf2,rad9
[h,p,stats] = fishertest(x);
fprintf('FE test p=%0.05f\n' , p ); 

figure('units','centimeters','position',[5 5 6 5]); 
hold on; 
grid on;
plot(T1.Hours(~isnan(T1.rad9)), T1.rad9(~isnan(T1.rad9)) , 'Color', [.4 .4 .4],'LineWidth',2.5 ,'DisplayName','{\itrad9}');
plot(T1.Hours(~isnan(T1.psf2rad9)) , T1.psf2rad9(~isnan(T1.psf2rad9))  , 'Color', clrs(4,:),'LineWidth',2.5 ,'DisplayName','{\itpsf2,rad9}');
xlim([0 150]);
set(gca,'ytick',0:20:100)
set(gca,'xtick',0:60:200)
legend('location','se')
print( '-dpsc2' , [FIGUREOUTDIR 'psf2,rad9.eps'] ) ; 
close; 



