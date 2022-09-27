function nfigure(n,x,y)
%
% fuction nfigure(n,x,y)
% Questa funzione permette di collocare le figure in ordine sullo schermo
% del pc come avviene con il comando subplot, solo che si utilizza tutta la
% dimensione del video.
%
% n è il numero della figura
% x numero di colonne di figure
% y numero di righe di figure

% clc
% clear all
% close all
% format compact


% x=10;
% y=4;
screen = get(0,'ScreenSize');
screen(4)=screen(4)-45;
min_x=screen(1);
min_y=screen(2);
max_x=screen(3);
max_y=screen(4);

d_x=(max_x-min_x)/x;        % dimensione x della figura
d_y=(max_y-min_y)/y;        % dimensione y della figura
quanti=n/(y*x);             % numero di schermate richieste per plottare le n figure

if n>y*x
    quanti=(n-1)/(y*x);
    quanti=floor(quanti);
    n_new=n-quanti*y*x;
else
    n_new=n;
end

resto=mod(n_new-1,x);
intero=(n_new-1-resto)/x;
% intero=x-1-intero;
% resto=y-1-resto;
% disp(['n= ' num2str(n) '   riga= ' num2str(intero) '   colonna= ' num2str(resto)])
Y=floor(max_y/2)-intero*d_y;

if y == 1
    Y = 0;
end

Y = Y+ 45;
X=resto*d_x;

figure(n)
set(n,'OuterPosition',[X,Y,d_x,d_y])

end