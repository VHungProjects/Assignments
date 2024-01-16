elist = load("20vhjh.txt");
n = max(elist(:));
[nRow, nCol] = size(elist);
A = zeros(n);

for row = 1:nRow
    node1 = elist(row, 1);
    node2 = elist(row, 2);

    A(node1,node2) = 1;
    A(node2,node1) = 1;
end

%Form the degree matrix
D1Vec = sum(A, 2);
D1 = diag(D1Vec);

%Form Laplacian Matrix
L1 = D1 - A;

%Calculate Eigenvalues and Eigenvectors
[eigVec1, eigVal1] = eig(L1);

%Cluster The Graph
fiedlerVec1 = eigVec1(:,2); %(all rows, column 2)

clusterNum1 = fiedlerVec1 >= 0; %Assign value of 1 if the element value is > 0 and 0 otherwise (ASSIGNMENT USE >=0)
clusterNum1 = (clusterNum1 * 2) - 1; %Ao cluster numbers become -1 and 1 instead of 0 and 1

set12 = clusterNum1;
set1 = [];
set2 = [];
for idx = 1:n
    if set12(idx) == 1
        set1 = [set1 , idx];
    else
        set2 = [set2 , idx];
    end
end

disp(set1)

%Use the function provided in the assignment to plot the grah with
%associated clusters
%plot271a1(A, clusterNum1)