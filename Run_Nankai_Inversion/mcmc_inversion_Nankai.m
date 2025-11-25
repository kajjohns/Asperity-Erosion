

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  BEGIN SETUP


%number of boundary nodes along strike 
num_Dpts = 30;
num_Zpts = 8;

%starting value -- radii of asperities 
radii = 15*ones(num_Dpts*num_Zpts,1);

%MCMC step size
stepsize_radii = 10*ones(size(radii));


%minimum radius allowed
minrad = 10;

%starting value -- accumulated stress on ring around asperities (Pa)
ring_taus = 0.15e6*ones(size(radii));  


%MCMC step size
stepsize_ring_tau = 0.15e6*ones(size(radii));



%time step for propgation 
dT = 10;

D = 8;  %ringth width (km)


%elastic modulus
mu = 3e10;


%continue sampling, or start new?
%setting to false overwrites output files
%setting to true appecnds cexisting files
continuing = false;

%use starting value from a different folder?
startval = false;
starting_folder = 'Nankai_outputs';

%folder name for saving output files
folder_name = 'Nankai_outputs';

%%%%  END SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ---- Software configuration
% Make hmmvp available
addpath ../hmmvp0.16
addpath ../tools

[gamb.hmat.id nnz] = hm_mvp('init',gamb.hmat.savefn,4);

%reweight so that vertical is not ignored
Sigu = 1*ones(size(Vu));




patch_stuff = make_triangular_patch_stuff(el,nd);
centroids = patch_stuff.centroids_faces;
strikes = patch_stuff.strike_faces;

%rotate and shift interface so that average strike is aligned with y axis
%this is needed for computing boundary positions
str = mean(strikes)*pi/180;
R = [cos(str) -sin(str); sin(str) cos(str)];
centroids_rot = (R*centroids(:,1:2)')';
centroids_rot(:,3) = centroids(:,3);
centroids_rot(:,2) = centroids_rot(:,2)-min(centroids_rot(:,2));


%node points along strike
Dpts = linspace(0,max(centroids_rot(:,2)),num_Dpts);
Zpts = [4 7 10 13 15 20 25 30];

%Zpts = [4 7 10 13 15 20];
[ptsD ptsZ] = meshgrid(Dpts,Zpts);
ptsDZ = [ptsD(:) ptsZ(:)];

%find nearest mesh centroid to node point
node_centroid = dsearchn([centroids_rot(:,2) -centroids_rot(:,3)],ptsDZ);

%rates in mm/yr
srate = rates/1000;  %m/yr

%adjust for stiffness
nel = size(el,1);



%no confidene in Sigbase
Sigbase(Sigbase<1) = 1;

%Greeens function matrices  
if invert_vel

    G = [Ge;Gn];
    GG = [Ge./repmat(Sige,1,size(Ge,2)); Gn./repmat(Sign,1,size(Gn,2)) ];
    dd = [Ve./Sige; Vn./Sign];

else

    G = [Gbase;Gu];
    GG = [Gbase./repmat(Sigbase,1,size(Gbase,2));Gu./repmat(Sigu,1,size(Gu,2))];
    dd = [Vbase./Sigbase; Vu./Sigu];

end


num = 1:size(el,1);
maxit = 200;
tol = 1e-4;

%need to scale stresses
%GFs computed with km, convert to meters
%GFs computed for mu=1
scale = mu*10^-3;



if continuing
eval(['load ./' folder_name '/M_radii.txt'])
radii =  M_radii(end,:)';

eval(['load ./' folder_name '/M_ring_tau.txt'])
ring_taus =  M_ring_tau(end,:)';

end


if startval
eval(['load ./' starting_folder '/M_radii.txt'])
radii =  M_radii(end,:)';

eval(['load ./' starting_folder '/M_ring_tau.txt'])
ring_taus =  M_ring_tau(end,:)';

end


 


%vector of unknown values
%starting values
X = [radii; ring_taus];
stepsize = [stepsize_radii; stepsize_ring_tau ];
%%%%%%

 



%forward model calculation
i_locked = get_locked_radii_indices2(node_centroid,ptsDZ,radii,centroids_rot,patch_stuff.dip_faces,minrad);

%accumulated stress
ring_tau = get_ring_tau2(node_centroid,ring_taus,ptsDZ,D,radii,centroids_rot,i_locked,patch_stuff.dip_faces,minrad);


asp_index = i_locked;
bslip = zeros(size(el,1),1);
bslip(asp_index) = srate(asp_index);  %m/yr, backslip rate
dtau_rate = scale*hm_mvp('mvp',gamb.hmat.id,bslip);
dtau = dtau_rate*dT; %accumulated stress (stressing rate times T1)

%include ring stress
dtau = -dtau + ring_tau;

%compute creep rate
rs = num(~asp_index);
s = gmres(@(x)mvp(x,gamb.hmat.id,rs,rs,~asp_index),-dtau(~asp_index)/scale,[],tol,maxit);
inc_slip = zeros(size(el,1),1);
inc_slip(~asp_index)=s; %increment of forward slip
U = srate*dT + inc_slip;  %srate*T1 to convert to displacement over time T1
U(asp_index)=0;  %locked asperity
creep_rate = U/dT*1000;  %conver from m/yr to mm/yr;


dhat = G*(rates-creep_rate); %unweighted
ddhat = GG*(rates-creep_rate);  %weighted by uncertainties

%solve for vertical shift
%Gs = [zeros(length(Vbase),1);1./Sigu];
%shift = Gs\(dd-ddhat);

%ddhat = ddhat + shift*[zeros(length(Vbase),1);1./Sigu];
%dhat = dhat + shift*[zeros(length(Vbase),1);ones(size(Vu))];

%add long-term velocities
%dhat = dhat + [lt_vels(:,1);lt_vels(:,2)];    
%ddhat = ddhat + [lt_vels(:,1)./Sige; lt_vels(:,2)./Sign ];
  

logrho = -.5*(dd-ddhat)'*(dd-ddhat);



%%%%%%

Xprev = X;
logrhoprev = logrho;
dhatprev = dhat;
i_locked_prev = i_locked;
radii_prev = radii;
ring_tau_prev = ring_tau;
creep_rate_prev = creep_rate;


if ~continuing
fid = fopen(['./' folder_name '/M_radii.txt'],'w'); fclose(fid);
fid = fopen(['./' folder_name '/M_ring_tau.txt'],'w'); fclose(fid);
fid = fopen(['./' folder_name '/logrho.txt'],'w'); fclose(fid);
fid = fopen(['./' folder_name '/dhat.txt'],'w'); fclose(fid);
fid = fopen(['./' folder_name '/locked_index.txt'],'w'); fclose(fid);
fid = fopen(['./' folder_name '/creep_rates.txt'],'w'); fclose(fid);
end

fidM_radii = fopen(['./' folder_name '/M_radii.txt'],'a'); 
fidM_tau = fopen(['./' folder_name '/M_ring_tau.txt'],'a'); 
fidrho = fopen(['./' folder_name '/logrho.txt'],'a'); 
fiddhat = fopen(['./' folder_name '/dhat.txt'],'a');
fid_locked = fopen(['./' folder_name '/locked_index.txt'],'a');
fid_creeping = fopen(['./' folder_name '/creep_rates.txt'],'a');


numparams = length(X);
numaccept = 0;

rand('state', sum(100*clock));
opts =  optimset('display','off');
for iter=1:10^6

    %step through all params and vary one at a time
    count = mod(iter,numparams)+1; 
    r=(-1)^round(rand(1))*rand(1);
    r=r*stepsize(count); 
    X(count)=X(count)+r;


    %%%%%%


    radii = X(1:length(radii));
    ring_taus = X(1+length(radii):end);
      
   
  

    
  if sum(ring_taus<0 | ring_taus>7e5)==0    %upper depth < lower depth
          
    
      
        
        
       
        %forward model calculation
        i_locked = get_locked_radii_indices2(node_centroid,ptsDZ,radii,centroids_rot,patch_stuff.dip_faces,minrad);
        
        %accumulated stress
        ring_tau = get_ring_tau2(node_centroid,ring_taus,ptsDZ,D,radii,centroids_rot,i_locked,patch_stuff.dip_faces,minrad);

        
        asp_index = i_locked;
        bslip = zeros(size(el,1),1);
        bslip(asp_index) = srate(asp_index);  %mm/yr, backslip rate
        dtau_rate = scale*hm_mvp('mvp',gamb.hmat.id,bslip);
        dtau = dtau_rate*dT; %accumulated stress (stressing rate times T1)
        
        %include ring stress
        dtau = -dtau + ring_tau;
        
        %compute creep rate
        rs = num(~asp_index);
        s = gmres(@(x)mvp(x,gamb.hmat.id,rs,rs,~asp_index),-dtau(~asp_index)/scale,[],tol,maxit,[],[],inc_slip(~asp_index));
        inc_slip = zeros(size(el,1),1);
        inc_slip(~asp_index)=s; %increment of forward slip
        U = srate*dT + inc_slip;  %srate*T1 to convert to displacement over time T1
        U(asp_index)=0;  %locked asperity
        creep_rate = U/dT*1000;  %conver from m/yr to mm/yr
        
        
        dhat = G*(rates-creep_rate);  
        ddhat = GG*(rates-creep_rate);  

        %solve for vertical shift
        %Gs = [zeros(length(Vbase),1);1./Sigu];
        %shift = Gs\(dd-ddhat)
        
        %ddhat = ddhat + shift*[zeros(length(Vbase),1);1./Sigu];
        %dhat = dhat + shift*[zeros(length(Vbase),1);ones(size(Vu))];


        logrho2 = -.5*(dd-ddhat)'*(dd-ddhat);


      


        accept=metropolis_log(1,logrho,1,logrho2);

    else

        accept = 0;

    end


    

    if accept==1

       logrho=logrho2;
      

       Xprev=X;
       dhatprev = dhat;
       logrhoprev = logrho;
       
        i_locked_prev = i_locked;
         radii_prev = radii;
        ring_tau_prev = ring_taus;
    
        creep_rate_prev = creep_rate;


       numaccept = numaccept+1;

    else

        X=Xprev;
        dhat=dhatprev;
        logrho=logrhoprev;
        i_locked = i_locked_prev;
        radii = radii_prev;
        ring_taus = ring_tau_prev;
        creep_rate = creep_rate_prev;

    

    end

    
    if count==1

        %iter

        fprintf(fidM_radii,'\n',' ');fprintf(fidM_radii,'%6.8f\t',radii');

        fprintf(fidrho,'\n',' ');fprintf(fidrho,'%6.8f\t',logrho);
        fprintf(fiddhat,'\n',' ');fprintf(fiddhat,'%6.8f\t',dhat');
        
        fprintf(fidM_tau,'\n',' ');fprintf(fidM_tau,'%6.8f\t',ring_taus');
        fprintf(fid_locked,'\n',' ');fprintf(fid_locked,'%6.8f\t',i_locked');

        fprintf(fid_creeping,'\n',' ');fprintf(fid_creeping,'%6.8f\t',creep_rate');


        disp(['Completed sample number ' num2str(iter)])
        disp(['Acceptance rate: ' num2str(numaccept/iter*100)])

                
        disp(['Log likelihood: ' num2str(logrho)])
                
    
    end

    Logrhos(iter) = logrho;

    figure(100)
    plot(Logrhos)
    drawnow
    
end



