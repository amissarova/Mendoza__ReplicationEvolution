%% look at coverage for all cerevisae gDNA sequencing data
% needs depth.tgz from the cluster
% uncompress into a folder and run this script
% the script does some filtering, makes some figures, and saves a .mat file w/the filtered experiments. 


%% load coverage data for all experiments
% download & uncompress depth.tgz as generated from the 'depth' target in 
%    /homes/users/lcarey/single_cell_behavior/Projects/2017__Mendoza_ReplicationTiming/AllYeastGDNAseqs

DATADIR = '~/Downloads/' ; 
EXPS = readtable([ DATADIR' exps.tab' ] ,'FileType','text','ReadVariableNames',false);
DEPTH = dlmread( [ DATADIR 'depth.tab' ] );
LOCS = readtable([ DATADIR' regions.tab' ],'FileType','text','ReadVariableNames',false);

LOCS.Properties.VariableNames = {'chr' 'start' 'stop'};

%remove experiments w/zero coverage at all positions (incomplete files)
EXPS  =  EXPS( ~all(DEPTH==0,2) ,:);
DEPTH = DEPTH( ~all(DEPTH==0,2) ,:);

% remove non-gDNA seqs
DEPTH = DEPTH( regexpcmp( LOCS.chr , '^chr[IVX]') ,:); 
LOCS = LOCS( regexpcmp( LOCS.chr , '^chr[IVX]') ,:); 

%% remove experiments w/uneven coverage
%   arbitrary filters, but the results seem rather insensitive to reasonable values

% filter to get experiments w/reasonable coverage
THRESH_MEAN_VS_MEDIAN = 1.05 ; % 1.05 very stringent , 1.5 loose
%THRESH_MEAN_VS_MEDIAN = 10.5 ; % 1.05 very stringent , 1.5 loose

THRESH_MAX_BINS_ZERO = 1000 ; % discard experiment if this many bins have zero reads


normalized_depth = DEPTH ./ repmat( sum(DEPTH) , nrows(DEPTH) , 1);
medians = median(normalized_depth,1) ; 
means   = mean(normalized_depth,1) ; 

idx_reasonable_distribution = (medians > (means./THRESH_MEAN_VS_MEDIAN)) & (medians < (means.*THRESH_MEAN_VS_MEDIAN)) ;  % stringent threshold determined by visual inspection of imagesc()

idx_not_missing_too_much = sum(normalized_depth==0) < THRESH_MAX_BINS_ZERO ; % no one that passes the 'reasonable distribution' threshold is missing > 200kb
idx = idx_not_missing_too_much & idx_reasonable_distribution ; 

EXPS_passed_thresh = EXPS.Var1(idx, :) ; 
normalized_depth = normalized_depth( : , idx ) ;
figure; imagesc( normalized_depth , [0 prctile(normalized_depth(:),95)] )

%%
R = table();
warning('off','MATLAB:table:RowsAddedExistingVars') ; 
THRESH = prctile(normalized_depth(:),1) ; 
chrs = unique(LOCS.chr);
c=0;
for chri = 1:numel(chrs)
    idx = find(strcmp(LOCS.chr,chrs{chri}));
    for I = 1:1000
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

% plot results
figure;
hold on; 
clrs = parula(numel(chrs));
BINSIZE = 200 ; 
markers = 'o+svo+svo+svo+sv';
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
save('AllYeastGDNAseqs.mat' , 'R' , 'normalized_depth'  , 'EXPS_passed_thresh' , 'LOCS') ;
