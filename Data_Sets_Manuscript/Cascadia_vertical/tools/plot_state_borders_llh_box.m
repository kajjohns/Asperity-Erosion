function plot_state_borders_llh_box(color,bbox)

addpath tl_2017_us_state 
S = shaperead('tl_2017_us_state.shp');

for j=1:56
    
    ind= S(j).X<bbox(2) &  S(j).X>bbox(1) & S(j).Y<bbox(4) &  S(j).Y>bbox(3);
     
    if sum(ind)>0
       plot(S(j).X,S(j).Y,'color',color)
    end

 end