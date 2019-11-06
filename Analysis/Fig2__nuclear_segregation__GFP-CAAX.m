%% 2019.10.17 -- new data on nuclear segregation and chr bridges
% Figure1_nuclear_segragation_and_bridges_new_data.m


%% load data
PROJECTDIR = '~/CareyLab/Projects/2017__MendozaReplication/' ;
DATADIR = '~/Develop/Mendoza__ReplicationEvolution/Data/' ; 
DATA_FILE_1 = [ DATADIR 'Figure2__NucSeg__GFP-CAAX__1335, July30-31 and August1.txt' ] ;
DATA_FILE_2 = [ DATADIR 'Figure2__NucSeg__GFP-CAAX__5066, July30-31 and August3.txt' ] ;

FIGUREOUTDIR = [ PROJECTDIR '/Figures_for_revision/' ] ;  system(['mkdir -p ' FIGUREOUTDIR]);

T1 = readtable( DATA_FILE_1 ) ; 

T2 = readtable( DATA_FILE_2 ) ; 

T1.noHUAll(T1.noHUAll==200) = 999 ; 
T1.withHUAll(T1.withHUAll==200) = 999 ; 

T2.noHUAll(T2.noHUAll==200) = 999 ; 
T2.withHUAll(T2.withHUAll==200) = 999 ; 


%% Main figure; 
clrHU = [0.7 0 0.7] ; 
lw = 2.5 ; 

figure('units','centimeters','position',[5 5 4.5 4]); 
hold on; 
grid on;
[f,x] = ecdf( T1.noHUAll ) ;
plot(x,100*f,'-','LineWidth',lw,'Color','k','DisplayName','WT');
[f,x] = ecdf( T1.withHUAll ) ;
plot(x,100*f,'-','LineWidth',lw,'Color',clrHU,'DisplayName','+HU');
xlim([0 160])
set(gca,'ytick',0:20:100)
set(gca,'xtick',0:50:1000)
lh = legend('location','se');
lh.Position = [0.6 0.1830 0.1 0.2] ;

n_no = sum(~isnan(T1.noHUAll))
n_HU = sum(~isnan(T1.withHUAll))

% Fisher's exact test
x = [ sum(T1.noHUAll<999) sum(T1.noHUAll==999) ; ...
    sum(T1.withHUAll<999) sum(T1.withHUAll==999) ]
[~,p,~] = fishertest( x ) 

% Mann?Whitney U test 
[p,~,~] = ranksum(T1.noHUAll , T1.withHUAll )

figure('units','centimeters','position',[5 5 4.5 4]); 
hold on; 
grid on;
[f,x] = ecdf( T2.noHUAll ) ;
plot(x,100*f,'-','LineWidth',lw,'Color','k','DisplayName', 'WT')
[f,x] = ecdf( T2.withHUAll ) ;
plot(x,100*f,'-','LineWidth',lw,'Color',clrHU,'DisplayName','+HU');
xlim([0 160])
set(gca,'ytick',0:20:100)
set(gca,'xtick',0:50:1000)
lh = legend('location','se');
lh.Position = [0.6 0.1830 0.1 0.2] ;

n_no = sum(~isnan(T2.noHUAll))
n_HU = sum(~isnan(T2.withHUAll))

x = [ sum(T2.noHUAll<999) sum(T2.noHUAll==999) ; ...
    sum(T2.withHUAll<999) sum(T2.withHUAll==999) ]
[~,p,~] = fishertest( x ) 
[p,~,~] = ranksum(T2.noHUAll , T2.withHUAll )

