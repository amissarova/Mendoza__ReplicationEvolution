%% Figure S7 showing individual replicates
%  S7 is YOYO bridges
%  
%

%% load data
DATADIR = '~/Develop/Mendoza__ReplicationEvolution/' ; 
DATAFILE = [ DATADIR 'Ivanova source data.xlsx']  ; 
A = readtable( DATAFILE , 'Sheet' , 'figS7a','Range','A7:M13','ReadRowNames',true);
B = readtable( DATAFILE , 'Sheet' , 'figS7b','Range','A7:M13','ReadRowNames',true);
clrs = cbrewer('qual','Dark2',8) ; 


%%
tem1 = table2array( A( '%YOYO in DAPI bridge' , regexpcmp(A.Properties.VariableNames,'Tem1')) )
cdc15 = table2array( A( '%YOYO in DAPI bridge' , regexpcmp(A.Properties.VariableNames,'Cdc15')) )
dbf2 = table2array( A( '%YOYO in DAPI bridge' , regexpcmp(A.Properties.VariableNames,'Dbf2')) )
mat = [tem1 ; cdc15 ; dbf2]; 
fh = figure('units','centimeters','position',[5 5 8 3]);
hold on ;
barh( mean(mat,2) , 'FaceColor', [.6 .6 .6] ) ;
herrorbar( mean(mat,2) , [1 2 3] , sem(mat,2),'.k'); 
plot( mat(1,:) , [1 1 1] ,'ok')
plot( mat(2,:) , [2 2 2] ,'ok')
plot( mat(3,:) , [3 3 3] ,'ok')
xlim([0 60])
ylim([0.5 3.5])
set(gca,'xtick',0:10:100)
xlabel('% YOYO-1 positive histone bridges')
set(gca,'yticklabel',{'\it{tem1-3}' '\it{cdc15-2}' '\it{dbf2-2}'})

print('-dpng','~/Downloads/S7A_barplot.png','-r600');


%%
%%
clear 'tem1' 'cdc15' 'dbf2' 'mat' ;

tem1 = table2array( B( '%YOYO  bridges' , regexpcmp(A.Properties.VariableNames,'Tem1')) )
cdc15 = table2array( B( '%YOYO  bridges' , regexpcmp(A.Properties.VariableNames,'Cdc15')) )
dbf2 = table2array( B( '%YOYO  bridges' , regexpcmp(A.Properties.VariableNames,'Dbf2')) )
mat = [tem1 ; cdc15 ; dbf2]; 
fh = figure('units','centimeters','position',[5 5 8 3]);
hold on ;
barh( mean(mat,2) , 'FaceColor', [.6 .6 .6] ) ;
herrorbar( mean(mat,2) , [1 2 3] , sem(mat,2),'.k'); 
plot( mat(1,:) , [1 1 1] ,'ok')
plot( mat(2,:) , [2 2 2] ,'ok')
plot( mat(3,:) , [3 3 3] ,'ok')
xlim([0 60])
ylim([0.5 3.5])
set(gca,'xtick',0:10:100)
xlabel('% YOYO-1 positive histone bridges')
set(gca,'yticklabel',{'\it{tem1-3}' '\it{cdc15-2}' '\it{dbf2-2}'})

print('-dpng','~/Downloads/S7B_barplot.png','-r600');
