%% generate DS with coverage info for all samples
% run once
cd ~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq
addpath(genpath('~/Develop/matlab'));
DS = dataset();

folder_names = {'~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/coverage_data_1000/G1_HU',...
    '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/coverage_data_1000/M_G1',...
    '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/coverage_data_1000/M_M',...
    '~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/coverage_data_1000/M_G1__noBrdU'};
folder_IDs = {'G1HU' , 'MG1' , 'MM' , 'MG1noBrdU'};
for Z = 1:length(folder_names)
    folder_name = folder_names{Z}; folder_ID = folder_IDs{Z};
    samples = dir(folder_name); samples = samples(3:end);
    cd(folder_name)
    for I = 1:length(samples)
        D = dataset('file' , samples(I).name);
        D = D(: , [1:4]);
        D.CN = NaN(length(D) , 1);
        D.pct = NaN(length(D) , 1);
        m = modefit( D.cov , 0 , [-5:5:1000] );
        if m > 0
            D.CN = D.cov/m;
            for J = 1:length(D)
                D.pct(J) = nanmean(D.cov < D.cov(J));
            end
        end
        
        D.Properties.VarNames{4} = strcat('cov_' , folder_ID , '_' , sprintf('%d',I));
        D.Properties.VarNames{5} = strcat('CN_' , folder_ID , '_' , sprintf('%d',I));
        D.Properties.VarNames{6} = strcat('pct_' , folder_ID , '_' , sprintf('%d',I));
        if isempty(DS)
            DS = D;
        else
            DS = join(DS , D , 'Type' , 'Left' , 'Keys' , {'chr' , 'start_point' , 'end_point'} , 'MergeKeys' , true);
        end
    end
end
% keep only gDNA chromosomes
idx = [];
for I = 1:length(DS)
    if findstr(DS.chr{I} , 'chr')
        idx = [idx I];
    end
end
DS = DS(idx,:);
% add info from underreplication dataset regarding each window (dist to the end, RT, dist to ARS)
D = DS;
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp_new.mat');
DS = DS(: , {'chr' , 'start_point' , 'end_point' , ...
    'Trep_spline' , 'dist_to_ARS' , 'dist_to_the_end_kb' , ...
    'percent_underreplicated_cdc20' , 'percent_underreplicated_dbf2'});
D.dist_to_the_end_kb = NaN(length(D) , 1);
D.Trep_spline = NaN(length(D) , 1);
D.dist_to_ARS = NaN(length(D) , 1);
D.percent_underreplicated_cdc20 = NaN(length(D) , 1);
D.percent_underreplicated_dbf2 = NaN(length(D) , 1);

for I = 1:length(D)
    idx = find(DS.start_point >= D.start_point(I) & DS.start_point < D.end_point(I) & strcmp(DS.chr , D.chr{I}));
    D.dist_to_the_end_kb(I) = nanmean(DS.dist_to_the_end_kb(idx));
    D.Trep_spline(I) = nanmean(DS.Trep_spline(idx));
    D.dist_to_ARS(I) = nanmean(DS.dist_to_ARS(idx));
    D.percent_underreplicated_cdc20(I) = nanmean(DS.percent_underreplicated_cdc20(idx));
    D.percent_underreplicated_dbf2(I) = nanmean(DS.percent_underreplicated_dbf2(idx));
end
DS = D;

% add rDNA region
DS.rDNA_bool = zeros(length(DS),1);
N1 = 459797; N2 = 468931;
idx = find(strcmp(DS.chr , 'chrXII'));
for J = 1:length(idx)
    if DS.start_point(idx(J)) >= N1 & DS.start_point(idx(J)) < N2
        DS.rDNA_bool(idx(J)) = 1;
    end
end
% save DS
cd ~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq
save('~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000.mat' , 'DS');

%% add wether region has transposables or tRNAs 
load('~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000.mat');
DS.tRNA_bool = zeros(length(DS) , 1);
DS.transposable_bool = zeros(length(DS) , 1);
D = DS;
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features_new.mat');
% tRNA
idx = find(strcmp(DS.TYPE , 'tRNA'));
for I = 1:length(idx)
    middle_point = (DS.start_point(idx(I)) + DS.end_point(idx(I)))/2;
    idx_point = find(strcmp(D.chr , DS.chr{idx(I)}) & middle_point >= D.start_point & middle_point < D.end_point);
    D.tRNA_bool(idx_point) = 1;
end
% transposables
idx = find(strcmp(DS.TYPE , 'transposable_element_gene') | strcmp(DS.TYPE , 'retrotransposon'));
for I = 1:length(idx)
    middle_point = (DS.start_point(idx(I)) + DS.end_point(idx(I)))/2;
    idx_point = find(strcmp(D.chr , DS.chr{idx(I)}) & middle_point >= D.start_point & middle_point < D.end_point);
    D.transposable_bool(idx_point) = 1;
end
DS = D;
save('~/Develop/Mendoza__ReplicationEvolution/BrdU_Seq/DS_1000.mat' , 'DS');
    

