function plot_coast(bbox)
%bbox [minlon minlat;maxlon maxlat]

S = shaperead('ne_10m_coastline','BoundingBox',bbox);

for k=1:length(S)
    
    plot(S(k).X,S(k).Y,'k')
    
end
