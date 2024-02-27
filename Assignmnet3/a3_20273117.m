function a3_20273117
% Function for CISC271, Winter 2024, Assignment #3

    %% Problem A
    A = csvread("wine.csv",0,1);
    A = A'; %Cultivar is now top row
    lvec = A(:,1);
    Xmat = A(:,2:end);
    % % Compute the pair of columns of Xmat with the lowest DB indexs
    dbscores = zeros(size(Xmat,2),size(Xmat,2));
    for i = 1:size(Xmat, 2)
        for j = 1:size(Xmat, 2)
            dbscores(i, j) = dbindex(Xmat(:, [i j]), lvec);
        end
    end
    [val, idx] = min(dbscores(:));
    [col1, col2] = ind2sub(size(dbscores), idx);
    % Scatter plot of the reduced data using gscatter
    %% part I
    disp("PROBLEM A")
    disp("Lowest DB score: " + val);
    disp("Columns " + col1 + " and " + col2);
    %% Part II
    % Plot the two columns with the lowest DB index
    %Zero mean
    figure; 
    gscatter(Xmat(:, col1), Xmat(:, col2), lvec);
    xlabel("Column 1");
    ylabel("Column 2");
    title('Dimensionality reduction with lowest DB index');
    
    %% Problem B
    % % Compute the PCA's of the data using the SVD; score the labelings
    XmatZM = Xmat - mean(Xmat);
    [~, S, V] = svd(XmatZM);
    Z2 = XmatZM*V(:,1:2); %Project onto first 2 components Week5Day3A
    score = dbindex(Z2,lvec);

    %% PArt I
    disp("PROBLEM B")
    disp("DB score of PCA: " + score);
    %% Part II
    figure;
    gscatter(Z2(:,1),Z2(:,2),lvec);
    xlabel("Component 1");
    ylabel("Component 2");
    title('Scatter Plot of PCA Components')
 
    %% Problem C
    xmatSTD = zscore(Xmat);
    xmatSTDZM = xmatSTD - mean(xmatSTD);
    [~, S, V] = svd(xmatSTDZM);
    ZSTD = xmatSTDZM*V(:,1:2);
    score = dbindex(ZSTD,lvec);
    %% Part I
    disp("PROBLEM C")
    disp("DB Score of standardized PCA: " + score);
    %% Part II
    figure;
    gscatter(ZSTD(:,1),ZSTD(:,2),lvec);
    xlabel("Component 1");
    ylabel("Component 2");
    title('Scatter Plot of Standardized PCA Components')



end
function score = dbindex(Xmat, lvec)
% SCORE=DBINDEX(XMAT,LVEC) computes the Davies-Bouldin index
% for a design matrix XMAT by using the values in LVEC as labels.
% The calculation implements a formula in their journal article.
%
% INPUTS:
%        XMAT  - MxN design matrix, each row is an observation and
%                each column is a variable
%        LVEC  - Mx1 label vector, each entry is an observation label
% OUTPUT:
%        SCORE - non-negative scalar, smaller is "better" separation

    % Anonymous function for Euclidean norm of observations
    rownorm = @(xmat) sqrt(sum(xmat.^2, 2));

    % Problem: unique labels and how many there are
    kset = unique(lvec);
    k = length(kset);

    % Loop over all indexes and accumulate the DB score of each label
    % gi is the label centroid
    % mi is the mean distance from the centroid
    % Di contains the distance ratios between IX and each other label
    D = [];
    for ix = 1:k
        Xi = Xmat(lvec==kset(ix), :);
        gi = mean(Xi);
        mi = mean(rownorm(Xi - gi));
        Di = [];
        for jx = 1:k
            if jx~=ix
                Xj = Xmat(lvec==kset(jx), :);
                gj = mean(Xj);
                mj = mean(rownorm(Xj - gj));
                Di(end+1) = (mi + mj)/norm(gi - gj);
            end
        end
        D(end+1) = max(Di);
    end

    % DB score is the mean of the scores of the labels
    score = mean(D);
end
