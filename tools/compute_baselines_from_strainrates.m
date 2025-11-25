function Vbase = compute_baselines_from_strainrates(xysites,obs_coords_xy,Exx_bs,Exy_bs,Eyy_bs,Vec_unit,baselines,L)


alpha = atan2(Vec_unit(:,2),Vec_unit(:,1));
D = 50;

%integrate strainrates along line
for j=1:size(baselines,1)
  
    xs = xysites(baselines(j,2),1) + cos(alpha(j))*( L(j)/D/2:L(j)/D:L(j) );
    ys = xysites(baselines(j,2),2) + sin(alpha(j))*( L(j)/D/2:L(j)/D:L(j) );
  
    Exx = griddata(obs_coords_xy(:,1),obs_coords_xy(:,2),Exx_bs,xs,ys);
    Exx = sum(Exx);

    Exy = griddata(obs_coords_xy(:,1),obs_coords_xy(:,2),Exy_bs,xs,ys);
    Exy = sum(Exy);

    Eyy = griddata(obs_coords_xy(:,1),obs_coords_xy(:,2),Eyy_bs,xs,ys);
    Eyy = sum(Eyy);

    %baselin velocity differential from intergration of strain rates
    Vbase(j) = L(j)/D*(Exx*cos(alpha(j))^2 + Exy*2*cos(alpha(j))*sin(alpha(j)) + Eyy*sin(alpha(j))^2);  %equation 1.21, Segall text

   

end





