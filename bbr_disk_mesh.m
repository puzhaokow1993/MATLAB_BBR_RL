function [p, t, bnd] = bbr_disk_mesh(nr)
%BBR_DISK_MESH  Concentric-ring triangulation of the unit disk B_1.
%   [p, t, bnd] = bbr_disk_mesh(nr) uses nr rings with 6*i nodes on ring i.
%   p   : N x 2 node coordinates
%   t   : T x 3 triangles (counter-clockwise)
%   bnd : indices of the boundary nodes (outermost ring)

pts = zeros(1 + 3*nr*(nr+1), 2);    % 1 + sum_{i=1}^{nr} 6i nodes
c = 1;                              % pts(1,:) = centre (0,0)
for i = 1:nr
    n  = 6*i;
    th = 2*pi*(0:n-1)'/n + (0.5*pi/n)*mod(i,2);
    pts(c+1:c+n, :) = (i/nr) * [cos(th), sin(th)];
    c = c + n;
end
p = pts;
t = delaunay(p(:,1), p(:,2));

% orient counter-clockwise
d = (p(t(:,2),1)-p(t(:,1),1)).*(p(t(:,3),2)-p(t(:,1),2)) ...
  - (p(t(:,3),1)-p(t(:,1),1)).*(p(t(:,2),2)-p(t(:,1),2));
t(d < 0, :) = t(d < 0, [1 3 2]);

N   = size(p,1);
bnd = (N - 6*nr + 1 : N)';
end
