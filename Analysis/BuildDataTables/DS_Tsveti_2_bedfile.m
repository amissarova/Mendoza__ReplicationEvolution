% convert matlab dataset to a bed file
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_Tsveti.mat')
WD = '~/Google Drive/CareyLab/Projects/2017__MendozaReplication/';
%load( [WD 'data_figures_final/G_Michi.mat' ] ) DS = G ; DS = DS( ~DS.weird_bool , :); % remove bad replicates
clear 'G';

SGD = loadSGDFeaturesDS();

GENOME = readtable('refNC_Genome.fasta.fai','FileType','text','Delimiter','\t','ReadVariableNames',false);
GENOME = GENOME( : , 1:2);
GENOME.Properties.VariableNames = {'chrom' 'nt'};

%% rewrite to get UCSC chr names
DS.UCSCchr = cell( size(DS,1)  , 1);
chr0N = 1132 ; 
basename = 'ref|NC_00' ;
for I = 1:size(DS,1)
    chr = DS.chr_num(I);
    DS.UCSCchr{I} = [ basename num2str(chr0N + DS.chr_num(I)) '|'];
end

%% calculate median across replicates
id = strcat(DS.MutantID , DS.chr); 
uids = unique(id);
R = table()
for I = 1:numel(uids)
    idx = strcmp( id , uids{I}) ;
    R.chr_name{I} = DS.UCSCchr{ find( idx , 1)} ; 
    R.id{I} = uids{I} ;
    R.chr_num(I) = DS.chr_num( find( idx , 1)) ; 
    R.chr{I} = DS.chr{ find( idx , 1)} ; 
    pct_unrep_cells_mat = cell2mat(DS.unreplicated_cells(idx)')' ;
    R.unreplicated_cells{I} = nanmedian(pct_unrep_cells_mat)';

    cn_mat = cell2mat(DS.CN(idx)')' ;
    R.cn{I} = nanmedian(cn_mat)';    
    
    R.start_point{I} = DS.start_point{ find( idx , 1)} ; 
    R.end_point{I} = DS.end_point{ find( idx , 1)} ; 
    
end
R = sortrows(R , 'chr_num' , 'ascend'); 
%% create a final table in BED format
clear 'T' ;
for I = 1:height(R)
    Q = table();
    N = numel(R.start_point{I});
    Q.chrom = repmat(R.chr_name(I) , N , 1);
    Q.chr = repmat(R.chr(I) , N , 1);
    Q.start_point = R.start_point{I};
    Q.end_point = R.end_point{I};
    Q.name = repmat(R.id(I) , N , 1);
    Q.cn = R.cn{I} ; 
    Q.unreplicated_cells = R.unreplicated_cells{I};
    if ~exist('T','var')
        T = Q ; 
    else
        T = vertcat( T , Q);
    end
end

%% save as .bed files
delete('dbf2.bed')
delete('cdc20.bed')
writetable( T( regexpcmp( T.name , 'dbf2') , {'chrom' 'start_point' 'end_point' 'name' 'cn'}) ...
    , 'dbf2.bed' ,'Delimiter', '\t'  , 'WriteVariableNames' , false , 'FileType','text' );
writetable( T( regexpcmp( T.name , 'cdc20') , {'chrom' 'start_point' 'end_point' 'name' 'cn'}) ...
    , 'cdc20.bed' ,'Delimiter', '\t'  , 'WriteVariableNames' , false , 'FileType','text' );

system('/usr/local/bin/bedtools map -a UCSC_sacCer3_all_refNC.bed -b dbf2.bed -o median | grep -v NC_001224 > dbf2_features.bed')
system('/usr/local/bin/bedtools map -a UCSC_sacCer3_all_refNC.bed -b cdc20.bed -o median| grep -v NC_001224 > cdc20_features.bed')


%% calc enrich
dbf2 = readtable('dbf2_features.bed','FileType','text','TreatAsEmpty',{'.'});
cdc20 = readtable('cdc20_features.bed','FileType','text','TreatAsEmpty',{'.'});

% optionally, truncate dataset at Xkb from ends
uchr = unique(cdc20.Var1);
TRUNCATE = 75 * 1000 ; 
cdc20 = cdc20( cdc20.Var3 < TRUNCATE ,:);
for I = 1:numel(uchr)
    maxpos = GENOME.nt( strcmp(GENOME.chrom , uchr{I}));
    idx1 = ~strcmp(cdc20.Var1 , uchr{I}) ; 
    idx2 = cdc20.Var2 < ( maxpos - TRUNCATE) ;
    cdc20 = cdc20( (idx1 | idx2) ,:);
end
%%
ut = unique(SGD.TYPE);
R = table();
R.TYPE = ut ;
X = cdc20.Var7 ; 
orfidx = find( ismember( cdc20.Var4 , SGD.ORF( strcmp(SGD.TYPE,'ORF'))));
for I = 1:numel(ut)
    R.idx{I} = find( ismember( cdc20.Var4 , SGD.ORF( strcmp(SGD.TYPE,ut{I}))));
    R.N(I) = numel(R.idx{I});
    R.mean(I) = nanmean( X( R.idx{I} ));
    [~ , R.p_all(I) ] = ttest2(  X( R.idx{I} ), X);
    [~ , R.p_orf(I) ] = ttest2(  X( R.idx{I} ), X(orfidx));

end
R = R( ~isnan(R.mean) & R.N>5,:);
R = sortrows( R , 'TYPE','ascend');
R = R( ismember( R.TYPE , {'telomere' 'transposable_element_gene' 'ARS' 'tRNA_gene'}) , :);
% R = R( R.p_orf < 1e-3 ,:);
R.TYPE = regexprep( R.TYPE , '_gene' ,'');
R.TYPE{ strcmp(R.TYPE,'transposable_element')} = 'transposon' ; 

%
NOT_SIG_CLR = [.8 .8 .8] ; 
PVAL_THRESH = 0.01 ; 
data = NaN( height(R) , max(R.N));
for I = 1:height(R)
    data(  I , 1:R.N(I)) = X( R.idx{I});
end
data  =  100*(2.*(1-data));

fh = figure('units','centimeters','position',[5 5 5 5 ]);
hold on ;
bh = boxplot( data' , 'PlotStyle','compact','symbol','' ,'Color','k','MedianStyle','target');
for I = 1:height(R)
    if R.p_orf(I) > PVAL_THRESH
        set(bh(1,I),'Color',NOT_SIG_CLR)
        set(bh(2,I),'Color',NOT_SIG_CLR)
        set(bh(3,I),'MarkerEdgeColor',NOT_SIG_CLR)
        set(bh(4,I),'MarkerEdgeColor',NOT_SIG_CLR)
        set(bh(5,I),'Color',NOT_SIG_CLR)
        set(bh(5,I),'Color',NOT_SIG_CLR)
    end
end
line(xlim,[0 0],'LineStyle','--','Color',[.7 .7 .7])
set(gca,'xtick',1:height(R))
set(gca,'xticklabel', regexprep( R.TYPE , '_' ,' '))
ylim([-10 70])
set(gca,'ytick',0:20:100)
rotateticklabel( gca , 8 , 45)
ylabel('% of cells unreplicated')

%% GO Enrichment
load('~/Google Drive/CareyLab//SingleCellFitness/Microscopy/data/YeastGO.mat');
GO = readtable('go_terms.tab','Delimiter','\t','FileType','text','ReadVariableNames',false);
cdc20 = readtable('cdc20_features.bed','FileType','text','TreatAsEmpty',{'.'});

%% GO term calculations
R = table();
R.FamilyName = {YeastGO.FamilyName}' ; 
R.p = NaN(height(R),1) ; 
R.mean = NaN(height(R),1) ;
R.FamilyID = [YeastGO.FamilyID]' ; 

orfidx = find( ismember( cdc20.Var4 , SGD.ORF( strcmp(SGD.TYPE,'ORF'))));
X = 100*(2.*(1-cdc20.Var7)); % conver to % unreplicated cells
for I = 1:numel(YeastGO)
    idx = ismember( cdc20.Var4 , YeastGO(I).ORF );
    [~,p] = ttest2( X(idx) ,  X(orfidx));
    R.p(I) = p ;
    R.mean(I) = nanmean(X(idx));
    R.N(I) = sum(idx);
    R.GTORFS(I) = nanmean(X(orfidx)) < nanmean(X(idx)) ; 
end

R = sortrows(R,'mean','ascend');
R = R(R.p < 0.001,:);
%%
R( R.GTORFS  , {'FamilyID' 'mean'})
%%
idx = R.N>5 ; 
figure; 
plot(log10(R.p(idx)), R.mean(idx),'.w')
text(log10(R.p(idx)), R.mean(idx),R.FamilyName(idx))