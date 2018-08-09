%% builds a single dataset w/rep timing from everyone's data
DATADIR = '~/Develop/Mendoza__ReplicationEvolution/Data/ExternalData/RepTiming/' ;

%% load our data
load( [DATADIR '../../DS_stat__200bp_new.mat']);

%% chr lengths
ChrLengths = readtable( [DATADIR 'genome.chrlengths'],'Delimiter','\t','FileType','text');
ChrLengths.Properties.VariableNames = {'chr_num' 'len'};

%% load .bed files made w/bedools map mean
K = readtable([ DATADIR 'Koren10_WT_200.txt'] );
K.Properties.VariableNames = {'chr_num' 'start_point' 'Koren10'} ; 
R = readtable([ DATADIR 'Raghuraman01_200.txt'] );
R.Properties.VariableNames = {'chr_num' 'start_point' 'Raghuraman01'} ; 
A = readtable([ DATADIR 'Alvino07_200.txt'] );
A.Properties.VariableNames = {'chr_num' 'start_point' 'Alvino07'} ; 


T = outerjoin(K,R,'Key',{'chr_num' 'start_point'},'MergeKeys',true);
T = outerjoin(T,A,'Key',{'chr_num' 'start_point'},'MergeKeys',true);
T = innerjoin( T , dataset2table(DS(:,{'chr_num' 'start_point' 'dist_to_the_end_kb' 'Trep_spline' 'G4' 'GC'...
    'percent_underreplicated_cdc20' 'percent_underreplicated_dbf2' })) ...
    );


%%
idx = ~isnan(T.Raghuraman01);
figure; 
scatter(  100*T.percent_underreplicated_cdc20(idx) , T.Raghuraman01(idx))
xlim([-20 60])
%%
xl = -30:60 ; 
ThreeClasses = zeros(height(T),1);
ThreeClasses( T.Raghuraman01 > prctile(T.Raghuraman01,80) ) = 1 ; 
ThreeClasses( T.Raghuraman01 > prctile(T.Raghuraman01,99) ) = 2 ; 
fh = figure('units','centimeters','position',[5 5 10 10]); 
hold on ;
for I = 0:2
    idx = R_ThreeClasses == I;
    [f,x] = ksdensity( T.percent_underreplicated_cdc20(idx)*100 , xl );
    plot(x,f,'LineWidth',3)
end
xlim([-10 60])
legend({'80% earliest' '20% latest' '1% latest'},'location','best')
title('Raghuraman 2001 (\alpha-factor & cdc7-1 block & release)')
xlabel('% under-rep MET3pr-CDC20')
ylabel('Fraction of genomic loci')
set(gca,'ytick',[0 max(get(gca,'ytick'))]);


ThreeClasses = zeros(height(T),1);
ThreeClasses( T.Koren10 < prctile(T.Koren10,20) ) = 1 ; 
ThreeClasses( T.Koren10 < prctile(T.Koren10,1) ) = 2 ; 
fh = figure('units','centimeters','position',[5 5 10 10]);
hold on ;
title('Koren 2010 (sort S & G1 cells)')
hold on ;
for I = 0:2
    idx = ThreeClasses == I;
    [f,x] = ksdensity( T.percent_underreplicated_cdc20(idx)*100 , xl );
    plot(x,f,'LineWidth',3)
end
xlim([-10 60])
legend({'80% earliest' '20% latest' '1% latest'},'location','best')
xlabel('% under-rep MET3pr-CDC20')
ylabel('Fraction of genomic loci')
set(gca,'ytick',[0 max(get(gca,'ytick'))]);
