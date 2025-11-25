

%vertical color scale
vumin = -3;
vumax = 3;
NC=64;
z1 = linspace(0,1,NC/2)';
z1f = flipud(z1);
cmap = [ones(NC/2,1) z1 z1; 
        z1f          z1f  ones(NC/2,1)];
cmap = flipud(cmap);


