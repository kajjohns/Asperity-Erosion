

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  BEGIN SETUP



%number of nodes along strike 
num_Dpts = 30;

%number of nodes down dip 
num_Zpts = 9;
maxZ = 40;  %maximum depth of centroid of asperity

%starting value -- radii of asperities 
radii = 8*ones(num_Dpts*num_Zpts,1);

%MCMC step size
stepsize_radii = 10*ones(size(radii));

%starting value -- accumulated stress on ring around asperities (Pa)
ring_taus = 0.15e6*ones(size(radii));  

%MCMC step size
stepsize_ring_tau = 0.15e6*ones(size(radii));

%time step for propgation 
dT = 10;

D = 8;  %ringth width (km)

%minimum radius size
minrad = 5;

%elastic modulus
mu = 3e10;

%continue sampling, or start new?
%setting to false overwrites output files
%setting to true appends existing files
continuing = false;

%use starting value from a different folder?
startval = false;
starting_folder = 'Hikurangi_outputs';

%folder name for saving output files
folder_name = 'Hikurangi_outputs';
             
%%%%  END SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ---- Software configuration
% Make hmmvp available
addpath ../hmmvp0.16
addpath ../tools

[gamb.hmat.id nnz] = hm_mvp('init',gamb.hmat.savefn,4);



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
Zpts = [4 7 10 12.5 13 15 20 25 30];
[ptsD ptsZ] = meshgrid(Dpts,Zpts);
ptsDZ = [ptsD(:) ptsZ(:)];



%slip rate in mm/yr
srate = rates/1000;
nel = size(el,1);



%Greeens function matrices  
if invert_vel

    G = [Ge;Gn;Gu];
    GG = [Ge./repmat(Sige,1,size(Ge,2)); Gn./repmat(Sign,1,size(Gn,2)); Gu./repmat(Sigu,1,size(Gu,2)); ];
    dd = [Ve./Sige; Vn./Sign; Vu./Sigu];

else

    G = [Gbase; Gu];
    GG = [Gbase./repmat(Sigbase,1,size(Gbase,2)); Gu./repmat(Sigu,1,size(Gu,2)) ];
    dd = [Vbase./Sigbase;Vu./Sigu];

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

 

%%%% Precomputations for geometry (shared by subfunctions)
thetas   = linspace(0,2*pi,50);
cos_t    = cos(thetas(:));
sin_t    = sin(thetas(:));
tan_dips = tan(patch_stuff.dip_faces*pi/180);

% Build KDTree once
KDTree = KDTreeSearcher([centroids_rot(:,2), -centroids_rot(:,3)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Forward model (initial calc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[i_locked, ring_tau] = get_locked_and_ring_tau(...
    ptsDZ, radii, ring_taus, D, centroids_rot, minrad, ...
    cos_t, sin_t, tan_dips, KDTree);

bslip = zeros(nel,1);
bslip(i_locked) = srate(i_locked);  
dtau_rate = scale * hm_mvp('mvp', gamb.hmat.id, bslip);
dtau = dtau_rate*dT; 
dtau = -dtau + ring_tau;

rs = num(~i_locked);
inc_slip = zeros(nel,1);
s = gmres(@(x)mvp(x,gamb.hmat.id,rs,rs,~i_locked), ...
          -dtau(~i_locked)/scale, [], tol, maxit);
inc_slip(~i_locked) = s;

U = srate*dT + inc_slip;
U(i_locked) = 0;
creep_rate = U/dT*1000;  

dhat  = G*(rates-creep_rate);  
ddhat = GG*(rates-creep_rate);  

logrho = -.5*(dd-ddhat)'*(dd-ddhat);

Xprev          = X;
logrhoprev     = logrho;
dhatprev       = dhat;
i_locked_prev  = i_locked;
radii_prev     = radii;
ring_tau_prev  = ring_tau;
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

rng('shuffle');  % modern RNG seeding
opts =  optimset('display','off');

Logrhos = zeros(1,1e6);

for iter=1:1e6

   
    %-------------------------------------------------------
    % Propose update (unchanged Metropolis step)
    %-------------------------------------------------------
    count = mod(iter,numparams)+1; 
    r = (-1)^round(rand(1)) * rand(1);
    r = r * stepsize(count); 
    X(count) = X(count) + r;
    
    radii     = X(1:length(radii));
    ring_taus = X(1+length(radii):end);
    
    %-------------------------------------------------------
    % Check parameter bounds
    %-------------------------------------------------------
    if all(ring_taus >= 0 & ring_taus <= 7e5)
        
        %-------------------------------------------------------
        % Forward model calculation (using optimized subfunctions)
        %-------------------------------------------------------
       [i_locked, ring_tau] = get_locked_and_ring_tau(...
    ptsDZ, radii, ring_taus, D, centroids_rot, minrad, ...
    cos_t, sin_t, tan_dips, KDTree);

        % Backslip
        bslip(:) = 0; % reuse preallocated vector
        bslip(i_locked) = srate(i_locked);
        dtau_rate = scale * hm_mvp('mvp', gamb.hmat.id, bslip);
        dtau = dtau_rate*dT;
        
        % Include ring stress
        dtau = -dtau + ring_tau;

        % Solve for creep rate
        rs = num(~i_locked);
        s = gmres(@(x)mvp(x,gamb.hmat.id,rs,rs,~i_locked), ...
                  -dtau(~i_locked)/scale, [], tol, maxit, [], [], inc_slip(~i_locked));
        inc_slip(:) = 0;
        inc_slip(~i_locked) = s;
        
        U = srate*dT + inc_slip;
        U(i_locked) = 0;
        creep_rate = U/dT*1000; % mm/yr
        
        % Predictions
        dhat  = G*(rates-creep_rate);
        ddhat = GG*(rates-creep_rate);
        
        % Log-likelihood
        logrho2 = -.5*(dd-ddhat)'*(dd-ddhat);

        %-------------------------------------------------------
        % Metropolis acceptance step (UNCHANGED)
        %-------------------------------------------------------
        accept = metropolis_log(1, logrho, 1, logrho2);
        
    else
        accept = 0;
    end
    

    

    %-------------------------------------------------------
    % Update state (UNCHANGED)
    %-------------------------------------------------------
    if accept == 1
        logrho = logrho2;
        Xprev  = X;
        dhatprev = dhat;
        logrhoprev = logrho;
        i_locked_prev = i_locked;
        radii_prev = radii;
        ring_tau_prev = ring_taus;
        creep_rate_prev = creep_rate;
        numaccept = numaccept + 1;
    else
        X         = Xprev;
        dhat      = dhatprev;
        logrho    = logrhoprev;
        i_locked  = i_locked_prev;
        radii     = radii_prev;
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
 
 % Plot only every 100 iterations
    if mod(iter,100) == 0
        figure(100); clf
        plot(Logrhos(1:iter))
        drawnow
    end
    
end



