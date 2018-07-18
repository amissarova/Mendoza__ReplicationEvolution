%% load data
PROJECTDIR = '~/Develop/Mendoza__ReplicationEvolution/' ; 
T = readtable([ PROJECTDIR 'Data/ExternalData/Lang11.xlsx'] );
load([ PROJECTDIR '/Data/DS_stat__features_new.mat'] );
T = innerjoin( T , dataset2table(DS) , 'LeftKey','Locus','RightKey','ORF') ; 
clear 'DS' ;
%%
Q = table();
Q.URA3_10_8 = T.URA3_10_8 ;
Q.x_RepTime = zscore(T.x_RepTime);
Q.x_GC = zscore(T.x_GC);
Q.percent_underreplicated_cdc20 = zscore(T.percent_underreplicated_cdc20);
Q.percent_underreplicated_cdc20_DM = zscore(T.percent_underreplicated_cdc20_DM);
Q.dist_to_ARS = zscore(T.dist_to_ARS);
Q.x_RecRate = zscore(T.x_RecRate);


GLM = stepwiseglm( Q , 'linear', 'ResponseVar' ,'URA3_10_8' , 'Distribution','normal' ) 