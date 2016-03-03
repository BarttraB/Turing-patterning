
L8n=load('/Volumes/Storage/FORTRAN/turing_AG/L_E4_2D_6.txt') ;

%to make the matrix of data
L8nm=matrix_maker2D_3(L8n,263);

%to solve the ubuntu dual screen issue D,
set(0, 'DefaultFigureRendererMode', 'manual')
set(0,'DefaultFigureRenderer','zbuffer')

set(gcf, 'PaperUnits','centimeters')
xSize=8;ySize=10;
set(gcf,'Position',[0 0 xSize*100 ySize*100]);


n=4;
m=4;
subhandles=panels(n,m);

    cmin=min(L8nm(20,40,:));
    cmax=max(L8nm(20,40,:));

%frame=1:13:300;    
    

for i=1:n*m;
hh=subplot(subhandles(i));
pic=pcolor(squeeze((L8nm(:,:,(1+13*(i-1)))))); 
colormap(jet); shading flat; caxis([cmin cmax]); axis square
set(gca,'XTick',[]) 
set(gca,'YTick',[])
p = get(hh,'pos'); 
set(hh,'pos',[p(1)-0.05 p(2:4)]) 
end

set(gca, 'LineWidth', 1.2)

h=colorbar;

set(h, 'Position', [.8314 .11 .0581 .8150])
set(gca, 'FontSize', 42)
annotation(gcf,'textbox',...
    [0.872389791183294 0.00832342449464922 0.138051044083526 0.0463733650416171],...
    'String',{'L [\muM]'},...
    'FontSize',42,...
    'FontName','Arial',...
    'FitBoxToText','off',...
    'LineStyle','none');

%saving figure to eps or pdf file
print -depsc panels_Turing2D_spots_4
print -dpdf panels_Turing2D_spots_4
