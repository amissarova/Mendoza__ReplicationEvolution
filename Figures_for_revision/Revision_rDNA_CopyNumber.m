% does rDNA copy number ( 2 * rDNA_coverage / gDNA_coverage ) increase from
% metaphase -> G1 or in metaphase under CDKinhibition? 
%
% for reviewer #2
% LBC 2020

% load original coverage txt files and annotation file
DATADIR = '/Volumes/MicroscopyAndSequencing/DataFromBarcelona/2017__Mendoza_ReplicationTiming/MendozaRevision_Coverage_rDNA/' ;

ANNO = readtable( [ DATADIR 'annotation.tab' ] , 'FileType','text' , 'ReadVariableNames',false);
ANNO.Properties.VariableNames = {'GSM' 'condition' 'file_ID'};
ANNO.condition_cat = regexprep( ANNO.condition , '_rep.*' ,'');
ANNO.cond2 = regexprep( regexprep(ANNO.condition_cat , 'arrest_' , '') , '_mutant.*' ,'');
ANNO.cond2 = regexprep( ANNO.cond2 , '_' ,' ');

T = readtable( [DATADIR 'chrXII_rDNA_coverage.txt']  , 'FileType','text' , 'ReadVariableNames',false , 'Delimiter','\t');
T.Properties.VariableNames = {'chr' 'start' 'stop' 'mean_coverage' 'bamfile'} ;
T.file_ID = regexprep( regexprep( T.bamfile , '.*BAM/' , '')  , '.bam' ,'');
T.exp_date = regexprep( regexprep( T.bamfile , '^../' , '')  , '/BAM.*' ,'');
T.exp_date = regexprep( T.exp_date , 'Data_' ,'');


T = innerjoin( T , ANNO , 'Key' , 'file_ID');


T.condition_cat = regexprep( T.condition , '_rep.*' ,'');

T.cond2 = regexprep( regexprep(T.condition_cat , 'arrest_' , '') , '_mutant.*' ,'');
T.cond2 = regexprep( T.cond2 , '_' ,' ');

T.exp_date = regexprep( T.exp_date , 'Data_' ,'');

%%

%% calculate coverage and make a scatterplot for each experiment
figname = '~/Downloads/coverage_plot.eps' ; delete(figname);
R = table();
file_ID = unique(T.file_ID);
R.file_ID = file_ID( ismember( file_ID , ANNO.file_ID(~regexpcmp(ANNO.condition_cat,'dbf2')))) ; % remove dbf2

for I = 1:height(R)
    idxThisExp = strcmp(T.file_ID , R.file_ID{I} ) ; 
    idxNoSubTel = T.start>75000 & T.start< (max(T.start)-75000) ;
    idxL = idxNoSubTel & T.start < 437500 ; 
    idxR = idxNoSubTel & T.start > 475000 ;
    idx = (idxL | idxR) & idxThisExp  ; 
    v = T.mean_coverage( idx ) ;
    R.cov(I) = modefit(v( v < prctile(v,98) )) ;
   % R.cov(I) = median(v( v < prctile(v,98) )) ;
    
    idx_rDNA = T.start > 451500 & T.start < 468000  & idxThisExp  ;
    r = T.mean_coverage( (idx_rDNA) & strcmp(T.file_ID , R.file_ID{I})) ;
    R.cov_rDNA(I) = modefit(r) ;
 %   R.cov_rDNA(I) = median(r) ;

    R.exp_date{I} = T.exp_date{ find( strcmp(T.file_ID,R.file_ID{I}),1) } ; 

   % if (false)
        fh = figure('units','centimeters','position',[5 5  15 8]) ;
        hold on ;
        plot( T.start(idx)./1000 , T.mean_coverage(idx) ,'.k');
        plot( T.start(idx_rDNA)./1000 , T.mean_coverage(idx_rDNA) ,'.-r');
        set(gca,'yscale','log');
        ylim( [prctile(v,5) , 1e4])
        line(xlim,[ R.cov_rDNA(I)  R.cov_rDNA(I)]);
        line(xlim,[ R.cov(I) R.cov(I)]);

        ylabel('# of reads or coverage in 200bp window')
        xlabel('position along chr XII')
        
        cond2 = ANNO.cond2{ strcmp(ANNO.file_ID,R.file_ID{I})};
        txt = sprintf('I=%d , %s,%s g=%0.01f rDNA=%0.01f CN_{rDNA}=%0.01f' , I , R.exp_date{I} , cond2  , R.cov(I) , R.cov_rDNA(I)  , 0.5*R.cov_rDNA(I)/R.cov(I) ); 
        title( txt ) ;
        
        print('-dpsc2' , figname , '-append');
        close ; 

                
        
  % end
end


R = innerjoin(R,ANNO,'Key','file_ID');
R.exp_date = categorical(R.exp_date);

R.rDNA_CN = 0.5 .* R.cov_rDNA ./ R.cov ; 



% summary plot of all three experiments
fh = figure('units','centimeters','position',[5 5  30 8]) ;

subplot(1,3,1)
idx = R.exp_date == '2017-09-15' ; 
gscatter(R.cov(idx) , R.rDNA_CN(idx) , { R.cond2(idx) T.cond2(idx)}  )
ylabel('rDNA coverage (norm)')
xlabel('chrXII (excluding rDNA) mode-fit coverage')
title('2017-09-15' )

subplot(1,3,2)
idx = R.exp_date == '2017-03-23' ; 
gscatter(R.cov(idx)  , R.rDNA_CN(idx) , { R.cond2(idx) T.cond2(idx)}  )
ylabel('rDNA coverage (norm)')
xlabel('chrXII (excluding rDNA) mode-fit coverage')
title('2017-03-23')

subplot(1,3,3)
idx = R.exp_date == '2017-10-06' ; 
gscatter(R.cov(idx) , R.rDNA_CN(idx) , { R.cond2(idx) T.cond2(idx)}  )
ylabel('rDNA coverage (norm)')
xlabel('chrXII (excluding rDNA) mode-fit coverage')
title('2017-10-06' )

orient(fh,'landscape');
 print('-dpsc2' , figname , '-append');
  close ; 

%% %% %% %% OLD %% %% 
%%


%
% %% old version with mean across the chromosome
% %T = readtable( [DATADIR 'chrXII_rDNA_coverage.txt']  , 'FileType','text' , 'ReadVariableNames',false , 'Delimiter','\t');
% %T.Properties.VariableNames = {'chr' 'start' 'stop'  'mean_coverage' 'bamfile'} ;
%
%
% T.file_ID = regexprep( regexprep( T.bamfile , '.*BAM/' , '')  , '.bam' ,'');
% T.exp_date = regexprep( regexprep( T.bamfile , '^../' , '')  , '/BAM.*' ,'');
%
% T.NormalizedCoverage  = NaN(height(T),1);
% for I = 2:3:height(T)
%     T.NormalizedCoverage(I) = T.mean_coverage(I) / mean( [ T.mean_coverage(I-1) T.mean_coverage(I+1)] ) ;
% end
% T.NormalizedCoverage = 2*T.NormalizedCoverage ;

%figure;
%boxplot( T.NormalizedCoverage , {T.exp_date  T.condition_cat }  )
%ylabel('rDNA copy number')
