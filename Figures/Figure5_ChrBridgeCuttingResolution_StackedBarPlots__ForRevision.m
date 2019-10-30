
%% test
clrs = cbrewer('qual','Dark2',8) ; 

B_resolve = [ 82.0779220779221 91.42857142857144 ] ; 

fh = figure('units','centimeters','position',[5 5 4 7]);
hold on ;
bar( [ 100 100] , 'FaceColor' , [.99 .99 .99]);
bar(1 , B_resolve(1) , 'FaceColor' ,clrs(1,:) ,'FaceAlpha',0.5);
bar(2 , B_resolve(2) , 'FaceColor' ,clrs(2,:) ,'FaceAlpha',0.5);
set(gca,'xtick',[])
xlim([0.5 2.5])
ylabel('% of cells with a visible bridge')

%%
A_resolve = [12.468 9.351 ] 
hold on ;
bar( [ 100 100] , 'FaceColor' , [.95 .95 .95]);
bar(1 , B_resolve(1) , 'FaceColor' ,clrs(3,:) ,'FaceAlpha',0.5);
bar(2 , B_resolve(2) , 'FaceColor' ,clrs(1,:) ,'FaceAlpha',0.5);
set(gca,'xtick',[])
xlim([0.5 2.5])
ylabel('% of cells with a visible bridge')

%%

A_resolved = [10.309 12.680  ] ; 

fh = figure('units','centimeters','position',[5 5 4 7]);
hold on 
bar( [100 100] , 'FaceColor',[.95 .95 .95] );
bar( 1 , A_resolved(1) , 'FaceColor', 'k', 'FaceAlpha',0.5);
bar( 2 , A_resolved(2) , 'FaceAlpha',0.5);
ylabel('% of cells with bridges')
set(gca,'xtick',[])
xlim([0.5 2.5])
%%
C_cut = [ 17.5824 18.813  36.7473 63.29670  49.934] ;
fh = figure('units','centimeters','position',[5 5 4 7]);
clrs = get(gca,'ColorOrder');
hold on 
bar( [100 100 100 100 100] , 'FaceColor',[.95 .95 .95] );
bar( 1 , 100-C_cut(1) , 'FaceColor', clrs(1,:), 'FaceAlpha',0.5);
bar( 2 , 100-C_cut(2) , 'FaceColor', clrs(2,:),  'FaceAlpha',0.5);
bar( 3 , 100-C_cut(3) , 'FaceColor', clrs(3,:),  'FaceAlpha',0.5);
bar( 4 , 100-C_cut(4) , 'FaceColor', clrs(4,:),  'FaceAlpha',0.5);
bar( 5 , 100-C_cut(5) , 'FaceColor', clrs(5,:),  'FaceAlpha',0.5);

ylabel('% of cells with bridges')
set(gca,'xtick',[])
xlim([0.5 5.5])
