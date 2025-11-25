
function [Vbase,Sigbase,BaseEnds,L,baselines,Vec_unit,crossind] = make_baseline_rate_changes(xy,Veast,Vnorth,Sigeast,Signorth,PatchEnds)



if nargin<6
    PatchEnds = [];
end


%plot data and triangles
 tri = delaunay(xy(:,1),xy(:,2));
 %remove triangles with small angles

for k=1:size(tri,1)
    v1 = [xy(tri(k,3),1) xy(tri(k,3),2)]-[xy(tri(k,1),1) xy(tri(k,1),2)]; 
    v2 = [xy(tri(k,2),1) xy(tri(k,2),2)]-[xy(tri(k,1),1) xy(tri(k,1),2)];
    a1 = 180/pi*acos( dot(v1/norm(v1),v2/norm(v2)) );
    v1 = [xy(tri(k,1),1) xy(tri(k,1),2)]-[xy(tri(k,2),1) xy(tri(k,2),2)]; 
    v2 = [xy(tri(k,3),1) xy(tri(k,3),2)]-[xy(tri(k,2),1) xy(tri(k,2),2)];
    a2 = 180/pi*acos( dot(v1/norm(v1), v2/norm(v2)) );
    v1 = [xy(tri(k,1),1) xy(tri(k,1),2)]-[xy(tri(k,3),1) xy(tri(k,3),2)]; 
    v2 = [xy(tri(k,2),1) xy(tri(k,2),2)]-[xy(tri(k,3),1) xy(tri(k,3),2)];
    a3 = 180/pi*acos( dot(v1/norm(v1),v2/norm(v2)) );
    
    if abs(a1)<10 | abs(a2)<10 | abs(a3)<10
        index(k)=1;
    else
        index(k)=0;
    end

end
index=logical(index);
tri(index,:)=[];


%form baselines
tri_num = (1:size(tri,1))';
tri_num = repmat(tri_num,3,1);
baselines = [ tri(:,1:2) ; tri(:,2:3) ; [tri(:,1) tri(:,3) ]];  %keep track of triangle number
baselines = sort(baselines,2);

%baseline endpoints
BaseEnds = [xy(baselines(:,1),:) xy(baselines(:,2),:)];

%identify baslines that cross patches (will use model strainrate value at center of baseline)



if ~isempty(PatchEnds)
    for k=1:size(BaseEnds,1)
        %int = intersectEdges(BaseEnds(k,:), SegEnds(30:34,:));
        int = intersectEdges(BaseEnds(k,:), PatchEnds);
        crossind(k) = sum(sum(~isnan(int)))>0;
    end
else
    crossind = [];
end



[baselines, ib, ic] = unique(baselines,'rows');  %baselines(ic) = original baselines
tri_num = tri_num(ib);
BaseEnds = BaseEnds(ib,:);
if ~isempty(crossind)
    crossind = crossind(ib);
end


%unit directional vectors
Vec_unit = xy(baselines(:,1),:)-xy(baselines(:,2),:);
L = sqrt(Vec_unit(:,1).^2+Vec_unit(:,2).^2);
Vec_unit = Vec_unit./[L L];  



%dot velocity vectors into unit direction 
Vbase2 = Veast(baselines(:,2)).*Vec_unit(:,1) + Vnorth(baselines(:,2)).*Vec_unit(:,2);
Vbase1 = Veast(baselines(:,1)).*Vec_unit(:,1) + Vnorth(baselines(:,1)).*Vec_unit(:,2);

%propagate errors
Sig_base2 = sqrt(Vec_unit(:,1).^2.*Sigeast(baselines(:,2)).^2 + Vec_unit(:,2).^2.*Signorth(baselines(:,2)).^2);
Sig_base1 = sqrt(Vec_unit(:,1).^2.*Sigeast(baselines(:,1)).^2 + Vec_unit(:,2).^2.*Signorth(baselines(:,1)).^2);


Vbase = Vbase1 - Vbase2;
%propagate error
Sigbase = sqrt(Sig_base2.^2 + Sig_base1.^2);

end

function point = intersectEdges(edge1, edge2)
%INTERSECTEDGES Return all intersections between two set of edges
%
%   P = intersectEdges(E1, E2);
%   returns the intersection point of lines L1 and L2. E1 and E2 are 1-by-4
%   arrays, containing parametric representation of each edge (in the form
%   [x1 y1 x2 y2], see 'createEdge' for details).
%   
%   In case of colinear edges, returns [Inf Inf].
%   In case of parallel but not colinear edges, returns [NaN NaN].
%
%   If each input is [N*4] array, the result is a [N*2] array containing
%   intersections of each couple of edges.
%   If one of the input has N rows and the other 1 row, the result is a
%   [N*2] array.
%
%   See also:
%   edges2d, intersectLines
%
%   ---------
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 31/10/2003.
%

%   HISTORY
%   19/02/2004 add support for multiple edges.
%   15/08/2005 rewrite a lot, and create unit tests
%   08/03/2007 update doc
%   21/09/2009 fix bug for edge arrays (thanks to Miquel Cubells)
%   03/08/2010 fix another bug for edge arrays (thanks to Reto Zingg)

%% Initialisations

% ensure input arrays are same size
N1  = size(edge1, 1);
N2  = size(edge2, 1);

% ensure input have same size
if N1~=N2
    if N1==1
        edge1 = repmat(edge1, [N2 1]);
        N1 = N2;
    elseif N2==1
        edge2 = repmat(edge2, [N1 1]);
    end
end

% tolerance for precision
tol = 1e-14;

% initialize result array
x0  = zeros(N1, 1);
y0  = zeros(N1, 1);


%% Detect parallel and colinear cases

% indices of parallel edges
%par = abs(dx1.*dy2 - dx1.*dy2)<tol;
par = isParallel(edge1(:,3:4)-edge1(:,1:2), edge2(:,3:4)-edge2(:,1:2));

% indices of colinear edges
% equivalent to: |(x2-x1)*dy1 - (y2-y1)*dx1| < eps
col = abs(  (edge2(:,1)-edge1(:,1)).*(edge1(:,4)-edge1(:,2)) - ...
            (edge2(:,2)-edge1(:,2)).*(edge1(:,3)-edge1(:,1)) )<tol & par;

% Parallel edges have no intersection -> return [NaN NaN]
x0(par & ~col) = NaN;
y0(par & ~col) = NaN;


%% Process colinear edges

% colinear edges may have 0, 1 or infinite intersection
% Discrimnation based on position of edge2 vertices on edge1
if sum(col)>0
    % array for storing results of colinear edges
    resCol = Inf*ones(size(col));

    % compute position of edge2 vertices wrt edge1
    t1 = edgePosition(edge2(col, 1:2), edge1(col, :));
    t2 = edgePosition(edge2(col, 3:4), edge1(col, :));
    
    % control location of vertices: we want t1<t2
    if t1>t2
        tmp = t1;
        t1  = t2;
        t2  = tmp;
    end
    
    % edge totally before first vertex or totally after last vertex
    resCol(col(t2<-eps))  = NaN;
    resCol(col(t1>1+eps)) = NaN;
        
    % set up result into point coordinate
    x0(col) = resCol;
    y0(col) = resCol;
    
    % touches on first point of first edge
    touch = col(abs(t2)<1e-14);
    x0(touch) = edge1(touch, 1);
    y0(touch) = edge1(touch, 2);

    % touches on second point of first edge
    touch = col(abs(t1-1)<1e-14);
    x0(touch) = edge1(touch, 3);
    y0(touch) = edge1(touch, 4);
end


%% Process non parallel cases

% process intersecting edges whose interecting lines intersect
i = find(~par);

% use a test to avoid process empty arrays
if sum(i)>0
    % extract base parameters of supporting lines for non-parallel edges
    x1  = edge1(i,1);
    y1  = edge1(i,2);
    dx1 = edge1(i,3)-x1;
    dy1 = edge1(i,4)-y1;
    x2  = edge2(i,1);
    y2  = edge2(i,2);
    dx2 = edge2(i,3)-x2;
    dy2 = edge2(i,4)-y2;

    % compute intersection points of supporting lines
    delta = (dx2.*dy1-dx1.*dy2);
    x0(i) = ((y2-y1).*dx1.*dx2 + x1.*dy1.*dx2 - x2.*dy2.*dx1) ./ delta;
    y0(i) = ((x2-x1).*dy1.*dy2 + y1.*dx1.*dy2 - y2.*dx2.*dy1) ./ -delta;
        
    % compute position of intersection points on each edge
    % t1 is position on edge1, t2 is position on edge2
    t1  = ((y0(i)-y1).*dy1 + (x0(i)-x1).*dx1) ./ (dx1.*dx1+dy1.*dy1);
    t2  = ((y0(i)-y2).*dy2 + (x0(i)-x2).*dx2) ./ (dx2.*dx2+dy2.*dy2);

    % check position of points on edges.
    % it should be comprised between 0 and 1 for both t1 and t2.
    % if not, the edges do not intersect
    out = t1<-tol | t1>1+tol | t2<-tol | t2>1+tol;
    x0(i(out)) = NaN;
    y0(i(out)) = NaN;
end


%% format output arguments

point = [x0 y0];

end




function b = isParallel(v1, v2, varargin)
%ISPARALLEL Check parallelism of two vectors
%
%   B = isParallel(V1, V2)
%   where V1 and V2 are two row vectors of length ND, ND being the
%   dimension, returns 1 if the vectors are parallel, and 0 otherwise.
%
%   Also works when V1 and V2 are two N-by-ND arrays with same number of
%   rows. In this case, return a N-by-1 array containing 1 at the positions
%   of parallel vectors.
%
%   Also works when one of V1 or V2 is N-by-1 and the other one is N-by-ND
%   array, in this case return N-by-1 results.
%
%   B = isParallel(V1, V2, ACCURACY)
%   specifies the accuracy for numerical computation. Default value is
%   1e-14. 
%   
%
%   Example
%   isParallel([1 2], [2 4])
%   ans =
%       1
%   isParallel([1 2], [1 3])
%   ans =
%       0
%
%   See also
%   vectors2d, isPerpendicular, lines2d
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2006-04-25
% Copyright 2006 INRA - CEPIA Nantes - MIAJ (Jouy-en-Josas).

%   HISTORY
%   2007-09-18 copy from isParallel3d, adapt to any dimension, and add psb
%       to specify precision
%   2007-01-16 fix bug
%   2009-09-21 fix bug for array of 3 vectors
%   2011-01-20 replace repmat by ones-indexing (faster)
%   2011-06-16 use direct computation (faster)

% default accuracy
acc = 1e-14;
if ~isempty(varargin)
    acc = abs(varargin{1});
end

% adapt size of inputs if needed
n1 = size(v1, 1);
n2 = size(v2, 1);
if n1 ~= n2
    if n1 == 1
        v1 = v1(ones(n2,1), :);
    elseif n2 == 1
        v2 = v2(ones(n1,1), :);
    end
end

% performs computation
if size(v1, 2) == 2
    % computation for plane vectors
    b = abs(v1(:, 1) .* v2(:, 2) - v1(:, 2) .* v2(:, 1)) < acc;
else
    % computation in greater dimensions 
    b = vectorNorm(cross(v1, v2, 2)) < acc;
end

end


