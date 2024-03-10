function plot271a1(Amat, cvec)
% PLOTCLUSTER(AMAT,CVEC) plots the adjacency matrix AMAT twice;
% first, as a Cartesian grid, and seconnd, by using binary clusters
% in CVEC to plot the graph of AMAT based on two circles
%
% INPUTS: 
%         Amat - NxN adjacency matrix, symmetric with binary entries
%         cvec - Nx1 vector of class labels, having 2 distinct values
% OUTPUTS:
%         none
% SIDE EFFECTS:
%         Plots into the current figure

% %
% % Part 1 of 2: plot the graph as a rectangle
% %

% Problem size
[m n] = size(Amat);

% Factor the size into primes and use the largest as the X size
nfact = factor(n);
nx = nfact(end);
ny = round(n/nx);

% Create a grid and pull apart into coordinates; offset Y by +2
[gx, gy] = meshgrid((1:nx) - round(nx/2), (1:ny) + 2);

% Offset the odd rows to diagram the connections a little better
for ix=1:2:ny
    gx(ix, :) = gx(ix, :) + 0.25*ix;
end

% The plot function needs simple vectors to create the graph
x = gx(:);
y = flipud(gy(:));

% Plot the graph of A using the Cartesian grid
plot(graph(tril(Amat, -1), 'lower'), 'XData', x, 'YData', y);
axis('equal');

% %
% % Part 2 of 2: plot the graph as pair of circles
% %
% Set up the X and Y coordinates of each graph vertex
xy = zeros(2, numel(cvec));

% Number of cluster to process
kset = unique(cvec);
nk = numel(kset);

% Base circle is radius 2, points are centers of clusters
bxy = 2*circlen(nk);

% Process each cluster
for ix = 1:nk
    jx = cvec==kset(ix);
    ni = sum(jx);
    xy(:, jx) = bxy(:, ix) + circlen(ni);
end

hold on;
plot(graph(Amat), 'XData', xy(1,:), 'YData', xy(2,:));
hold off;
title(sprintf('Clusters of (%d,%d) nodes', ...
    sum(cvec==kset(1)), sum(cvec==kset(2))));
end

function xy = circlen(n)
% XY=CIRCLEN(N) finds N 2D points on a unit circle
%
% INPUTS:
%         N  - positive integer, number of points
% OUTPUTS:
%         XY - 2xN array, each column is a 2D point

xy = [cos(2*pi*(0:(n-1))/n) ; sin(2*pi*(0:(n-1))/n)];
end
