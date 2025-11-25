function bs_rate = compute_creeprate(id,idx_rbig,idx_rsmall,nel,i_locked,ring_tau,srate,bs_rate_prev)

dT = 10;
% Scale stresses (km → m, mu = 1 in GF calc)
scale = 3e10 * 1e-3;

%gmres specifications
maxit = 200;
tol   = 1e-4;


asp_index = logical(i_locked);
bslip = zeros(size(asp_index));
bslip(asp_index) = srate(asp_index);  %mm/yr, backslip rate

dtau_rate = scale*hm_mvp('mvp',id,bslip,idx_rbig,idx_rbig);
dtau = dtau_rate*dT; 

%include ring stress
dtau = -dtau + ring_tau*dT;

num = (1:length(i_locked))';

ind = ~asp_index & ismember(num,idx_rbig);
rs = num(ind);


s = gmres(@(x)mvp(x,id,rs,rs,ind),-dtau(ind)/scale,[],tol,maxit);


U = srate(ind)*dT + s;  %srate*T1 to convert to displacement over time T1
creep_rate = U/dT;  % m/yr 
cr = zeros(num(end),1);
cr(ind) = creep_rate;

bs = srate - cr; %m/yr

bs_rate = bs_rate_prev;
bs_rate(idx_rsmall) = bs(idx_rsmall)*1000;  %convert to mm/yr;

