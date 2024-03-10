function ph = plotline(vvec, color, lw, nv)
% PLOTLINE(VVEC,COLOR,LW,NV) plots a separating line
% into an existing figure
% INPUTS:
%        VVEC   - (M+1) augmented weight vector
%        COLOR  - character, color to use in the plot
%        LW   - optional scalar, line width for plotting symbols
%        NV   - optional logical, plot the normal vector
% OUTPUT:
%        PH   - plot handle for the current figure
% SIDE EFFECTS:
%        Plot into the current window. 

% Set the line width
if nargin >= 3 & ~isempty(lw)
  lwid = lw;
else
  lwid = 2;
end

% Set the normal vector
if nargin >= 4 & ~isempty(nv)
  do_normal = true;
else
  do_normal = false;
end

% Current axis settings
axin = axis();

% Scale factor for the normal vector
sval = 0.025*(axin(4) - axin(3));

% Four corners of the current axis
ll = [axin(1) ; axin(3)];
lr = [axin(2) ; axin(3)];
ul = [axin(1) ; axin(4)];
ur = [axin(2) ; axin(4)];

% Normal vector, direction vector, hyperplane scalar
nlen = norm(vvec(1:2));
uvec = vvec/nlen;
nvec = uvec(1:2);
dvec = [-uvec(2) ; uvec(1)];
bval = uvec(3);

% A point on the hyperplane
pvec = -bval*nvec;

% Projections of the axis corners on the separating line
clist = dvec'*([ll lr ul ur] - pvec);
cmin = min(clist);
cmax = max(clist);

% Start and end are outside the current plot axis, no problem
pmin = pvec +cmin*dvec;
pmax = pvec +cmax*dvec;

% Create X and Y coordinates of a box for the current axis
xbox = [axin(1) axin(2) axin(2) axin(1) axin(1)];
ybox = [axin(3) axin(3) axin(4) axin(4) axin(3)];

% Intersections of the line and the box
[xi, yi] = polyxpoly([pmin(1) pmax(1)], [pmin(2) pmax(2)], xbox, ybox);

% Point midway between the intersections
pmid = [mean(xi) ; mean(yi)];

% Range of the intersection line
ilen = 0.5*norm([(max(xi) - min(xi)) ; (max(yi) - min(yi))]);

% Plot the line according to the color specification
hold on;
if ischar(color)
    ph = plot([pmin(1) pmax(1)], [pmin(2) pmax(2)], ...
        [color '-'], 'LineWidth', lwid);
else
    ph = plot([pmin(1) pmax(1)], [pmin(2) pmax(2)], ...
        'Color', color, 'LineStyle', '-', 'LineWidth', lwid);
end
if do_normal
    quiver(pmid(1), pmid(2), nvec(1)*ilen*sval, nvec(2)*ilen*sval, ...
        'Color', color, 'LineWidth', lwid, ...
        'MaxHeadSize', ilen/2, 'AutoScale', 'off');
end
hold off;
