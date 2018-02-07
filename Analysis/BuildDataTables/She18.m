%% look at coverage for all cerevisae gDNA sequencing data
% needs depth.tgz from the cluster
% uncompress into a folder and run this script
% the script does some filtering, makes some figures, and saves a .mat file w/the filtered experiments. 


%% load coverage data for all experiments
% download & uncompress depth.tgz as generated from the 'depth' target in 
%    /homes/users/lcarey/single_cell_behavior/Projects/2017__Mendoza_ReplicationTiming/AllYeastGDNAseqs/depth.tgz

DATADIR = '~/Desktop/She18/' ; 
DATASAVEDIR  = DATADIR ; 
%DATASAVEDIR = '~/Develop/Mendoza__ReplicationEvolution/Data/' ; 
EXPS = readtable([ DATADIR 'exps.tab' ] ,'FileType','text','ReadVariableNames',false);
DEPTH = dlmread( [ DATADIR 'depth200.tab' ] )' ; %note the transpose
LOCS = readtable([ DATADIR 'regions200.tab' ],'FileType','text','ReadVariableNames',false);

LOCS.Properties.VariableNames = {'chr' 'start' 'stop'};
LOCS.Chromosome = str2double( regexprep(  regexprep( LOCS.chr , 'ref.NC_','')  ,'\|','') ) - 1132 ; 
EXPS = EXPS.Var1' ; 

% remove non-gDNA locations
% DEPTH = DEPTH( regexpcmp( LOCS.chr , '^chr[IVX]') ,:); 
% LOCS = LOCS( regexpcmp( LOCS.chr , '^chr[IVX]') ,:); 


Npos = height(LOCS) ; 
Nexp = numel(EXPS)  ; 
%% remove experiments w/uneven coverage
%   arbitrary filters, but the results seem rather insensitive to reasonable values

% filter to get experiments w/reasonable coverage
THRESH_MEAN_VS_MEDIAN = 1.05 ; % 1.05 very stringent , 1.5 loose
THRESH_MEAN_VS_MEDIAN = 1.25 ; % 1.05 very stringent , 1.5 loose

THRESH_MAX_BINS_ZERO = 1000 ; % discard experiment if this many bins have zero reads


normalized_depth = DEPTH ./ repmat( sum(DEPTH) , Npos , 1);
medians = median(normalized_depth,1) ; 
means   = mean(normalized_depth,1) ; 

idx_reasonable_distribution = (medians > (means./THRESH_MEAN_VS_MEDIAN)) & (medians < (means.*THRESH_MEAN_VS_MEDIAN)) ;  % stringent threshold determined by visual inspection of imagesc()
idx_enough_few_reads = sum(DEPTH) > 10^5 ; 
%idx_not_missing_too_much = sum(normalized_depth==0) < THRESH_MAX_BINS_ZERO ; % no one that passes the 'reasonable distribution' threshold is missing > 200kb
idx =  idx_reasonable_distribution & idx_enough_few_reads ; 

EXPS_passed_thresh = EXPS(idx) ; 
normalized_depth = normalized_depth( : , idx ) ;
DEPTH = DEPTH( : , idx ) ; 

normalized_depth_vect = normalized_depth(:);

figure; imagesc( normalized_depth , [0 prctile(normalized_depth(:),95)] )

normalized_depth_nan = normalized_depth;
normalized_depth_nan( normalized_depth_nan == 0) = NaN ; 
idx_sufficient_coverage_when_present = nanmedian(normalized_depth_nan , 2) > 1e-5 ; 
%% many positions have low coverage across all reads (GC content?)
%   these are comming up as false positives
% only interested in regions that, when not zero, have decent coverage
%

%% for a few thresholds, which regions tend to have < THRESH reads
LOCS.SuffCov = idx_sufficient_coverage_when_present ; 
for THRESH = [ 0 1 2  ]
    LOCS.(['n_times_region_lost_' num2str(THRESH)]) =  sum( DEPTH <= THRESH , 2) ;
end
N = size(DEPTH,2) ; 

%%
LOCS.KB = round(LOCS.start ./ 1000) ;
G = grpstats( LOCS , {'chr' 'Chromosome' 'KB' },'sum' ,'DataVars' , {'n_times_region_lost_0' 'SuffCov'});
%%
uchrs = unique(LOCS.chr);
for I = 1:17
    figure; hold on ;
    idx =  strcmp(LOCS.chr,uchrs{I}) ; 
    idx2 = idx & LOCS.SuffCov ; 
    plot( LOCS.start(idx2) , 100*(LOCS.n_times_region_lost_0(idx2) ./ N) ,'.' ,'DisplayName','reads <= 0');
    xlim([-10 max(LOCS.stop(idx2))])
  title( uchrs{I} )
  xlabel(sprintf( 'position along chr %d' , I) )
  ylabel('% of strains lost')
 legend('location','best')
end
%% build 'chr truncated till x' table
THRESH = prctile(normalized_depth(:),1) ; %threshold for considering # reads == 0 

R = table();
warning('off','MATLAB:table:RowsAddedExistingVars') ; 
chrs = unique(LOCS.chr);
c=0;
for chri = 1:numel(chrs)
    idx = find(strcmp(LOCS.chr,chrs{chri}));
    for I = 1:min(1000,numel(idx))
        missing_this_much = find(all(normalized_depth( idx(1:I) , : )  <= THRESH, 1)) ; 
        if any(missing_this_much)
            c=c+1;
            R.chr{c} = chrs{chri} ;
            R.topos(c) = LOCS.stop(idx(I)) ;
            R.nmissing(c) =  numel(missing_this_much);
            R.whomissing{c} =  EXPS_passed_thresh(missing_this_much) ;
            fprintf( '%s\t0-%d\t%d\t' , chrs{chri} , LOCS.stop(idx(I)) ,  numel(missing_this_much) );
            for exp = missing_this_much
                fprintf('%s\t' , EXPS_passed_thresh{exp})
            end
            fprintf('\n');
            
        end
    end
end

%% plot results
figure;
hold on; 
clrs = parula(numel(chrs));
BINSIZE = 200 ; 
markers = 'o+svo+svo+svo+svo';
for I = 1:numel(chrs)
    idx = find(strcmp(R.chr,chrs{I}));
    if ~isempty(idx)
    X = [ R.topos(idx)' max(R.topos(idx))+BINSIZE]; 
    Y =  [ R.nmissing(idx)' 0] ; 
    Y = (Y ./ ncols(normalized_depth) ) * 100  ; % % of strains
    plot( X ./ 1000 , Y,'-' ,'Marker',markers(I), 'Color',clrs(I,:),'DisplayName',chrs{I} ,'LineWidth',3);
    end
end
ylabel('% of strains missing this')
xlabel('Kb from left end')
set(gca,'xscale','log')
axis tight;
set(gca,'xtick',[0.1 0.2 0.5 1 2 5 10])
legend('location','ne')

%% save results
save( [  DATASAVEDIR 'AllYeastGDNAseqs.mat' ] , 'R' , 'normalized_depth'  , 'EXPS_passed_thresh' , 'LOCS') ;
