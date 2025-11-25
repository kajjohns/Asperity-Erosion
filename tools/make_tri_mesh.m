function [el,nd] = make_tri_mesh(contrkm,intv)


%% identify the points of the shallowest and deepest contours
% nodet: the shallowest points
% nodeb: the deepest points
ind = (contrkm(:,3) == contrkm(1,3));% hardwired since first element is shallowest
nodet0 = contrkm(ind,:);
ind = (contrkm(:,3) == contrkm(end,3));
nodeb0 = contrkm(ind,:);

% here I remove edge points to make sure that griddata works properly
nodet0(end,:) = [];
nodet0(1,:) = [];
nodeb0(end,:) = [];
nodet0(1,:) = [];


%% set the intervals of shallowest and deepest points as specified 
N = round(curvlength(nodet0)./intv);
nodet = curvspace(nodet0,N);
N = round(curvlength(nodeb0)./intv);
nodeb = curvspace(nodeb0,N);


%% mesh generation
[nd0,el] = meshfrac2(nodet,nodeb,intv);
% nd0  are coordinates of nodes on the mesh
% el   are the elements, the values specify the x,y,z indices of the nodes

%% shift the height of the nodes based on the contour info
% the mesh previously made is a straight-line interpolation from the top
% nodes to bottom nodes
[xi,yi,zi] = griddata(contrkm(:,1),contrkm(:,2),contrkm(:,3),...
    nd0(:,1),nd0(:,2),'linear');
ind = find(nd0(:,3) ~= nd0(1,3) & nd0(:,3) ~= nd0(end,3));
nd = nd0;
nd(ind,:) = [xi(ind),yi(ind),zi(ind)];

end

function q = curvspace(p,N)

% CURVSPACE Evenly spaced points along an existing curve in 2D or 3D.
%   CURVSPACE(P,N) generates N points that interpolates a curve
%   (represented by a set of points) with an equal spacing. Each
%   row of P defines a point, which means that P should be a n x 2
%   (2D) or a n x 3 (3D) matrix.
%
%   (Example)
%   x = -2*pi:0.5:2*pi;
%   y = 10*sin(x);
%   z = linspace(0,10,length(x));
%   N = 50;
%   p = [x',y',z'];
%   q = curvspace(p,N);
%   figure;
%   plot3(p(:,1),p(:,2),p(:,3),'*b',q(:,1),q(:,2),q(:,3),'.r');
%   axis equal;
%   legend('Original Points','Interpolated Points');
%
%   See also LINSPACE.
%

%   22 Mar 2005, Yo Fukushima


%% initial settings %%
currentpt = p(1,:); % current point
indfirst = 2; % index of the most closest point in p from curpt
len = size(p,1); % length of p
q = currentpt; % output point
k = 0;

%% distance between points in p %%
for k0 = 1:len-1
   dist_bet_pts(k0) = distance(p(k0,:),p(k0+1,:));
end
totaldist = sum(dist_bet_pts);

%% interval %%
intv = totaldist./(N-1);

%% iteration %%
for k = 1:N-1
   
   newpt = []; distsum = 0;
   ptnow = currentpt;
   kk = 0;
   pttarget = p(indfirst,:);
   remainder = intv; % remainder of distance that should be accumulated
   while isempty(newpt)
      % calculate the distance from active point to the most
      % closest point in p
      disttmp = distance(ptnow,pttarget);
      distsum = distsum + disttmp;
      % if distance is enough, generate newpt. else, accumulate
      % distance
      if distsum >= intv
         newpt = interpintv(ptnow,pttarget,remainder);
      else
         remainder = remainder - disttmp;
         ptnow = pttarget;
         kk = kk + 1;
         if indfirst+kk > len
            newpt = p(len,:);
         else
            pttarget = p(indfirst+kk,:);
         end
      end
   end
   
   % add to the output points
   q = [q; newpt];
   
   % update currentpt and indfirst
   currentpt = newpt;
   indfirst = indfirst + kk;
   
end

end


function l = distance(x,y)

% DISTANCE Calculate the distance.
%   DISTANCE(X,Y) calculates the distance between two
%   points X and Y. X should be a 1 x 2 (2D) or a 1 x 3 (3D)
%   vector. Y should be n x 2 matrix (for 2D), or n x 3 matrix
%   (for 3D), where n is the number of points. When n > 1,
%   distance between X and all the points in Y are returned.
%
%   (Example)
%   x = [1 1 1];
%   y = [1+sqrt(3) 2 1];
%   l = distance(x,y)
%

% 11 Mar 2005, Yo Fukushima

%% calculate distance %%
if size(x,2) == 2
   l = sqrt((x(1)-y(:,1)).^2+(x(2)-y(:,2)).^2);
elseif size(x,2) == 3
   l = sqrt((x(1)-y(:,1)).^2+(x(2)-y(:,2)).^2+(x(3)-y(:,3)).^2);
else
   error('Number of dimensions should be 2 or 3.');
end

end

function newpt = interpintv(pt1,pt2,intv)

% Generate a point between pt1 and pt2 in such a way that
% the distance between pt1 and new point is intv.
% pt1 and pt2 should be 1x3 or 1x2 vector.

dirvec = pt2 - pt1;
dirvec = dirvec./norm(dirvec);
l = dirvec(1); m = dirvec(2);
newpt = [intv*l+pt1(1),intv*m+pt1(2)];
if length(pt1) == 3
   n = dirvec(3);
   newpt = [newpt,intv*n+pt1(3)];
end

end

function totlen = curvlength(pts)

% CURVLENGTH Calculate the length of a curve.
%   CURVLENGTH(PTS) calculates the length of a curve defined
%   by number of points PTS. PTS should either be a nx2 (2D)
%   or a nx3 (3D) matrix.
%
%   (Algorithm)
%   The distance between each neighbouring points are summed.
%
%   (Example)
%   t = linspace(0,2*pi,100);
%   x = cos(t); y = sin(t);
%   pts = [x;y]';
%   l = curvlength(pts)
%   figure; plot(pts(:,1),pts(:,2)); axis equal;
%
%   See also DISTANCE.
%

totlen = 0;
for k = 1:size(pts,1)-1
    
       
        
     
    l = distancept(pts(k,:),pts(k+1,:));

    
    totlen = totlen + l;
end

end

function L = distancept(p1,p2)


D=p1-p2;
L = sqrt(D(:,1).^2+D(:,2).^2+D(:,3).^2);

end


function [nd,el,nodetnew] = meshfrac2(nodet,nodeb,intv,curvaxispts,vertcurv)

% MESHFRAC Mesh a fracture from specified top and bottom nodes
%
%   [nd,el] = meshfrac(nodet,nodeb,intv) generates a
%   mesh (nodes: nd, elements: el) by interpolating the
%   top and bottom nodes specified by nodet (top) and
%   nodeb (bottom), with average node interval intv.
%
%   nodet, nodeb should be n x 3 vectors where each row
%   represents a point (x,y,z). nodet and nodeb do not
%   have to have the same points.
%
%   The nodes specified by nodet and nodeb will be kept
%   unless vertcurv ~= 0 (in this case, nodet will be
%   displaced; see explanation of vertcurv.)
%
%   [nd,el,nodetnew] = meshfrac(...,curvaxispts,vertcurv) vertically
%   bends the mesh.
%
%   curvaxispts, a [m x 3] vector, where m is the number of
%   points, defines the mesh bending axis. The mesh is bent
%   perpendicular to the line connecting a point in curvaxispts
%   and the corresponding point in nodeb (for this purpose,
%   the points in curvaxispts are interpolated to have the
%   same number of points as nodeb).
%   When curvaxispts is identical to nodet, then the mesh
%   has the top and bottom points at the points specified by
%   nodet and nodeb; if not, then the top points do not
%   correspond to nodet.
%   In this case, nodetnew defines the new top nodes.
%
%   vertcurv defines the vertical curvature.
%   It is the angle between the line connecting a point in
%   curvaxispts and the corresponding point in nodeb and the
%   tangential line of the parabola at these points.
%   When positive and if nodet and nodeb are written from
%   south to north, it is convex downward.
%
%
%   (Example 1)
%   nodet = [1000,1000,-1000;1100,1050,-950;...
%            1180,1120,-900;1250,1220,-870];
%   nodeb = [777,635,-2000;963,789,-1960;1080,888,-1900;...
%            1186,955,-1870;1338,1076,-1850;...
%            1476,1154,-1800;1628,1250,-1770];
%   intv = 100;
%   [nd,el] = meshfrac(nodet,nodeb,intv);
%   figure;trisurf(el,nd(:,1),nd(:,2),nd(:,3));axis equal;
%   xlabel('x');ylabel('y');zlabel('z');
%
%   (Example 2: introducing vertical curvature)
%   nodet = [1000,1000,-1000;1100,1050,-950;...
%            1180,1120,-900;1250,1220,-870];
%   nodeb = [777,635,-2000;963,789,-1960;1080,888,-1900;...
%            1186,955,-1870;1338,1076,-1850;...
%            1476,1154,-1800;1628,1250,-1770];
%   intv = 100;
%   curvaxispts = [0,0,0;1000,1500,0];
%   vertcurv = 30;
%   [nd,el] = meshfrac(nodet,nodeb,intv,curvaxispts,vertcurv);
%   figure;trisurf(el,nd(:,1),nd(:,2),nd(:,3));axis equal;
%   xlabel('x');ylabel('y');zlabel('z');
%
%
%   See also meshdike, meshquadrangle.
%
%      12 Aug 2009, (c)Yo Fukushima
%

% 12 Aug 2009, made a modification so that now this is compatible with
% the cases where nodet and nodet are at the (nearly) same height.
%
% 18 Aug 2005, bug fixed and code cleaned. Now it does not do unnecessary
% calculations when vertcurv=0.
%
% 17 Aug 2005, the n values in the following lines are changed.
% Originally, they were one value larger (curvspace(ptnew,n+2) for example)
% if n <= 1, n = 1; end
% nodes = zeros(size(nodet_intpl,1),n,3);
%     foo = curvspace(ptnew,n+1);
%     xi = foo(1:n,1);
%     yi = foo(1:n,2);
%     zi = foo(1:n,3);
%
% 15 Aug 2005, bug fixed. changed to following around l. 182.
%     foo = squeeze(nodes(:,j,:));
%     totlen = curvlength(foo);


%% interpolate nodet (-> num. pts identical to nodeb) %%
nodet_intpl = curvspace(nodet,size(nodeb,1));

%% when nargin == 3 %%
if nargin == 3
    vertcurv = 0;
    curvaxispts = nodet_intpl;
else
    curvaxispts = curvspace(curvaxispts,size(nodeb,1));
end

%% allocate dimensions of nodes %%
% nodes(i,:,:) constitute of points along a line
% connecting a point in nodet_intpl and that in
% nodeb.
mheight = abs(mean(distancept(nodet_intpl,nodeb)));
% mheight = abs(mean(nodet_intpl(:,3) - nodeb(:,3)));
n = round(mheight./cos(vertcurv*pi/180)./intv);
if n <= 1, n = 1; end
nodes = zeros(size(nodet_intpl,1),n,3);

if vertcurv ~= 0
    %% calc. the bending direction %%
    if size(curvaxispts,1) >= 3
        middle = round(size(curvaxispts,1)/2);
        p1 = curvaxispts(middle,:);
        p2 = nodeb(middle-1,:);
        p3 = nodeb(middle+1,:);
        %         sizenodeb = size(nodeb,1)
        %         sizecurvaxispts = size(curvaxispts,1)
        %         middle
    else
        p1 = curvaxispts(1,:);
        p2 = nodeb(1,:);
        p3 = nodeb(2,:);
    end
    nvec = planenormvec(p1,p2,p3);
    nvec = nvec'./norm(nvec);
    if size(nvec,1)>size(nvec,2)
        nvec = nvec';
    end
    
    %% make it so that curve is convex downward %%
    crossprod = cross(p3-p2,nvec);
    if crossprod(3) > 0
        nvec = -nvec;
    end
end

%% filling nodes and move nodet %%
for i = 1:size(nodet_intpl,1)
    
    endpt1 = [nodeb(i,1),nodeb(i,2),nodeb(i,3)];
    endpt2 = [nodet_intpl(i,1),nodet_intpl(i,2),nodet_intpl(i,3)];
    
    % % when dip > 90, change the sign of nvec
    % % it allows to have the same curving direction
    % if dip > 90
    %     nvec = -nvec;
    % end
    
    if vertcurv ~= 0
        % move points between endpt1 and endpt2 so that
        % they are placed along a parabola
        N = 20;
        ptorig = [linspace(endpt1(1),endpt2(1),N)',...
            linspace(endpt1(2),endpt2(2),N)',...
            linspace(endpt1(3),endpt2(3),N)'];
        
        L = distancept(curvaxispts(i,:),endpt1);
        a = tan(vertcurv*pi/180)./L;
        midpt = mean([endpt1;curvaxispts(i,:)]);
        l = a.*((L/2)^2 - (distancept(midpt,ptorig)).^2);% length between old (lined) and new (parabola) pt
        ptnew = ptorig - l*nvec;
        % make points evenly spaced along each vertical curve
        foo = curvspace(ptnew,n+1);
        % update nodet (it has moved when vertcurv~=0)
        nodettmp(i,:) = foo(n+1,:);
    else
        foo = [linspace(endpt1(1),endpt2(1),n+1)',...
            linspace(endpt1(2),endpt2(2),n+1)',...
            linspace(endpt1(3),endpt2(3),n+1)'];
    end
    
    %% add to nodes %%
    xi = foo(1:n,1);
    yi = foo(1:n,2);
    zi = foo(1:n,3);
    nodes(i,:,:) = [xi, yi, zi];
    
end

if vertcurv ~= 0
    nodet = curvspace(nodettmp,size(nodet,1));
end


%% interpolate nodes evenly (horizontally) and create nd and nddummy %%
% nd is the resulting nodes, nddummy is used in delaunay.
% It has the same dimensions as nd.

nd = nodeb;
yidummy = linspace(0,10,size(nodeb,1));
zidummy = ones(1,size(nodeb,1)).*1.*10;
nddummy = [yidummy',zidummy'];

% if nodet(1,2) < nodet(size(nodet,1),2) % if nodes south -> north
%     st = 1;
%     ed = size(nodeb,1);
% else % if nodes north -> south
%     st = size(nodeb,1);
%     ed = 1;
% end

st = 1;
ed = size(nodeb,1);

for j = 2:size(nodes,2)
    %     totlen = 0;
    %     for i = 1:length(nodeb)-1
    %         try
    %             len = sqrt((nodes(i,j,1)-nodes(i+1,j,1)).^2+(nodes(i,j,2)-nodes(i+1,j,2)).^2);
    %         catch
    %             len = 0;
    %         end
    %         totlen = totlen + len;
    %     end
    foo = squeeze(nodes(:,j,:));
    totlen = curvlength(foo);
    n = ceil(totlen/intv);
    if n <= 1, n = 2; end
    if st > ed
        foo = nodes(ed:st,j,:);
        ndnew = curvspace(foo,n);
        ndnew = flipud(ndnew);
    else
        foo = nodes(st:ed,j,:);
        ndnew = curvspace(nodes(st:ed,j,:),n);
    end
    nd = [nd; ndnew];
    
    % yidummy, zidummy are used in delaunay
    yidummy = linspace(0,10,size(ndnew,1));
    zidummy = ones(1,size(ndnew,1)).*j.*10;
    nddummy = [nddummy; yidummy' zidummy'];
end

%% add nodet %%
% if nodet(1,2) < nodet(size(nodet,1),2) % if nodes south -> north
%     nd = [nd; nodet];
% else % if nodes north -> south
%     nd = [nd; flipud(nodet)];
% end

nd = [nd;nodet];
yidummy = linspace(0,10,size(nodet,1));
zidummy = ones(1,size(nodet,1)).*(size(nodes,2)+1).*10;
nddummy = [nddummy; yidummy' zidummy'];


%% delaunay %%
v = version;
if v(1:3) == '6.0'
    el = delaunay(nddummy(:,1),nddummy(:,2),0.00001);
else
    try
        el = delaunay(nddummy(:,1),nddummy(:,2));
    catch
        el = delaunay(nddummy(:,1),nddummy(:,2),{'Qt','Qbb','Qc','Qz'});
    end
end

%% return nodetnew %%
if nargout == 3
    nodetnew = nodet;
end

end
