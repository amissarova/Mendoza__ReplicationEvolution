%% load data
% nocomma has only mutational positions with a single alternative allele
%   all.tab has all mutational positions
DATADIR = '~/Data/Peter18/' ;
%T = readtable([ DATADIR 'nocomma.tab' ] ,'FileType','text','Delimiter','\t','Format','%d%d%s%s');
T = readtable([ '~/Downloads/' 'all.tab' ] ,'FileType','text','Delimiter','\t','Format','%d%d%s%s');
T.Properties.VariableNames = {'chr_num' 'pos' 'ref' 'alt'};

YG = fastaread('~/Data/Yeast/genome.fasta');
%% keep only SNPs
% T = T( cellfun(@length,T.alt)==1 , :); % removes SNPs w/multiple alt  alleles
T = T( cellfun(@length,T.ref)==1 , :);
T = T( cellfun(@length,T.alt)==1 |  ( regexpcmp(T.alt,',')  & cellfun(@length,T.alt)==3 ), :); % removes mutations w/multiple alt  alleles or indels


%% split lines w/a comma
has_comma = find( regexpcmp(T.alt,','));
new_alt = cellfun( @(X)X(end) ,  T.alt(has_comma)) ; 
prim_alt = cellfun( @(X)X(1) ,  T.alt(has_comma)) ; 
T.alt(has_comma) = arrayfun(@(X){X},prim_alt) ; 
R = table();
R.chr_num = T.chr_num(has_comma);
R.pos = T.pos(has_comma);
R.ref = T.ref(has_comma);
R.alt = arrayfun(@(X){X},new_alt) ; 
T = vertcat( T , R) ; 
T = sortrows(T , {'chr_num' 'pos'});

%% Define Subtelomeric regions as SUBTEL kb from the end
SUBTEL = 25*1000 ; 
T.Subtel = zeros(height(T),1);
T.Subtel( T.pos < SUBTEL) = 1 ; 
T.Subtel( T.pos > (max(T.pos)-SUBTEL)) = 1 ; 


% %% summary stats. GroupCount will be the # of occurances of each mutation
% G=grpstats(T,{'chr_num' 'ref' 'alt' 'Subtel'}, 'sum','DataVars','chr_num');
% G.mut = arrayfun(@(I)sprintf('%s->%s',G.ref{I},G.alt{I}),1:height(G),'UniformOutput',false)' ; 
% G.sum_chr_num = []; 
% %% compute stats for each mutation & chr
% G.FE_p = NaN(height(G),1);
% G.FE_or = NaN(height(G),1);
% for I = 1:2:height(G)
%     st = sum(G.GroupCount( G.Subtel==1 & G.chr_num==G.chr_num(I)));
%     nonst = sum(G.GroupCount( G.Subtel==0 & G.chr_num==G.chr_num(I)));
%     [~,p,or] = fishertest( [ nonst st ; G.GroupCount(I) G.GroupCount(I+1)] );
%    G.FE_p(I) = log10(p);
%    G.FE_or(I) = or.OddsRatio ;
% end
% 
% 
% G2 = grpstats( G( ~isnan(G.FE_p) ,:) , {'ref' 'alt'},{'median' 'std'},'DataVars',{'FE_p' 'FE_or'});
% G2 = sortrows(G2,'median_FE_p')



%% Mutational context
% confirm that the ref allele corresponds to this genome version
%  find(~(YG(1).Sequence(T.pos(T.chr_num==1)) == [T.ref{T.chr_num==1}]))
Q2 = table();
for I = 1:16
    Q = T(T.chr_num==I,:);
    find(~(YG(I).Sequence(Q.pos) == [Q.ref{:}]))
    Q.context = arrayfun( @(X)YG(I).Sequence([X-1 X X+1]) , Q.pos ,'UniformOutput',false);
    Q2 = vertcat(Q2,Q);
end
T = Q2; 
clear 'Q' 'Q2' 
% %% compute stats for each mutation & chr
% G=grpstats(T,{'chr_num' 'context' 'alt' 'Subtel'}, 'sum','DataVars','chr_num');
% G.mut = arrayfun(@(I)sprintf('%s->%s',G.context{I},G.alt{I}),1:height(G),'UniformOutput',false)' ; 
% G.sum_chr_num = []; 
% G.FE_p = NaN(height(G),1);
% G.FE_or = NaN(height(G),1);
% for I = 1:2:height(G)
%     st = sum(G.GroupCount( G.Subtel==1 & G.chr_num==G.chr_num(I)));
%     nonst = sum(G.GroupCount( G.Subtel==0 & G.chr_num==G.chr_num(I)));
%     [~,p,or] = fishertest( [ nonst st ; G.GroupCount(I) G.GroupCount(I+1)] );
%    G.FE_p(I) = log10(p);
%    G.FE_or(I) = or.OddsRatio ;
% end
% 
% 
% G2 = grpstats( G( ~isnan(G.FE_p) ,:) , {'context' 'alt'},{'median' 'std'},'DataVars',{'FE_p' 'FE_or'});
% G2 = sortrows(G2,'median_FE_p')
%% not splitting by chr
G=grpstats(T,{ 'context' 'ref' 'alt' 'Subtel'}, 'sum','DataVars','chr_num');
G2=grpstats(T,{'Subtel'}, 'sum','DataVars','chr_num');

G.mut = arrayfun(@(I)sprintf('%s->%s',G.context{I},G.alt{I}),1:height(G),'UniformOutput',false)' ; 
G.sum_chr_num = []; 
G.FE_p = NaN(height(G),1);
G.FE_or = NaN(height(G),1);
for I = 1:2:height(G)
    [~,p,or] = fishertest( [ G2.GroupCount(2) G2.GroupCount(1) ; G.GroupCount(I+1) G.GroupCount(I)] );
   G.FE_p(I) = p;
   G.FE_or(I) = or.OddsRatio ;
end
G = G( ~isnan(G.FE_p),:);
%% plot using FDR calculations
alpha =  0.1 ; 
[G.FDR, G.FDR_Q ] = mafdr(G.FE_p) ;
G.FDR_Storey =  mafdr( G.FE_p ,'BHFDR',true);
G.p_Bon = (G.FE_p * height(G)) ;

G.sig = G.FDR_Storey < alpha ; 

G = sortrows(G,{'ref' 'context' 'alt'});
uc = unique(G.context);
[~,o]=sort(cellfun(@(X)X(2),uc));
uc = uc(o);
ltrs = 'ACTG';
fh = figure('units','centimeters','position',[5 5 25 25]); 
clrs = parula(4);
for I = 1:numel(uc)
    subplot(8,8,I)
    hold on ;
    X = G( strcmp(G.context, uc{I}) & ~isnan(G.FE_p),:) ; 
    p = find(X.sig) ; 
    bar( 1:3 , X.FE_or ,'FaceColor',(clrs( find(ltrs==uc{I}(2)) , :)).^0.3  );
    if ~isempty(p)
        bar( p ,  X.FE_or(p) ,'FaceColor',clrs( find(ltrs==uc{I}(2)) , :).^2);
        X(p,:)
    end
    set(gca,'xtick',1:3)
    set(gca,'xticklabel',X.alt')
    if isempty(p)
        text(1 , 1.4 , X.context{1} );
    else
        text(1 , 1.4 , X.context{1} ,'FontWeight','bold');
    end
    ylim([0.5 1.5]);
    set(gca,'ytick',1)
    set(gca,'yticklabel',[])
    if any(I == 1:8:64) 
        set(gca,'ytick',[0.5 1 1.5])
        set(gca,'yticklabel',[0.5 1 1.5])
    end
    line( xlim , [1 1] ,'LineStyle','--','Color',[.7 .7 .7])
end
    

