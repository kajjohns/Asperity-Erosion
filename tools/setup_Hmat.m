

% ---- Software configuration
% Make hmmvp available
addpath ./hmmvp0.16;

global gamb;
gamb.use = 1;
% Optional H-matrix compression of the BEM matrix:
% Use H-matrix compression as implemented in hmmvp?
gamb.hmat.use = 1;
% Element-wise relative error of approximation. Should be ~ RelTol in the
% ODE integration.
gamb.hmat.rerr = 1e-3;
% Directory to which to save H-mat for future use:
dr = './hmat';


% ODE software parameters:
% File to which to save ODE solution
ne = size(el,1);
fn_dec = sprintf('ne%dhmat%drerr%3.2f',...
       ne,gamb.hmat.use,log10(gamb.hmat.rerr));


% H-matrix file name. 
gamb.hmat.savefn = sprintf('%s/Hmat_%s.dat',dr,fn_dec);

%compute new h-matrix, if requested 
if compute_hmat
    %make stress-slip matrix, stressG, such that stress = stressG*slip;
    CalcG_rake('Make',el,nd,rakes,gamb.hmat.rerr,gamb.hmat.savefn);
end

