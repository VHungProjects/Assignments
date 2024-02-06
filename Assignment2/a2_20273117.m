function [rmsvars lowIndexPositive lowIndexNegative rmstrain rmstest] = a2_20273117
% [RMSVARS LOWNDX RMSTRAIN RMSTEST]=A3 finds the RMS errors of
% linear regression of the data in the associated CSV file. The
% individual RMS errors are returned in RMSVARS and the index of the
% smallest RMS error is returned in LOWNDX. For the variable that
% best explains the dependent variable, a 5-fold cross validation is
% computed. The RMS errors for the training of each fold are returned
% in RMSTEST and the RMS errors for the testing of each fold are
% returned in RMSTEST.
%
% INPUTS:
%         none
% OUTPUTS:
%         RMSVARS  - 1xN array of RMS errors of linear regression
%         LOWINDEXPOSITIVE - integer scalar, index into RMSVALS
%         LOWINDEXNEGATIVE - integer scalar, index into RMSVALS
%         RMSTRAIN - 1x5 array of RMS errors for 5-fold training
%         RMSTEST  - 1x5 array of RMS errors for 5-fold testing
 
    [rmsvars lowIndexPositive lowIndexNegative] = a2q1
    [rmstrain rmstest] = a2q2(lowIndexPositive)

end


function [rmsvars lowIndexPositive lowIndexNegative] = a2q1
% [RMSVARS LOWNDX]=A2Q1 finds the RMS errors of
% linear regression of the data in the CSV file. The
% individual RMS errors are returned in RMSVARS and the index of the
% smallest RMS error is returned in LOWNDX. 
%
% INPUTS:
%         none
% OUTPUTS:
%         RMSVARS  - 1xN array of RMS errors of linear regression
%         LOWNDX   - integer scalar, index into RMSVALS

    % Read the test data from a CSV file and find the size of the data
    [fragilityVector,dataMatrix,countries,ageMatrix] = fragilitydata;
    [rows cols] = size(dataMatrix);
    %fragilityvector is second column <-- This is the dependent variable
    %dataMatrix is all the data normalized without frag and countries
    %countries is countries
    %ageMatrix are a 2xN of the age ranges
    % Compute the RMS errors for linear regression
    datastand = zscore(dataMatrix);%standardize
    cstand = zscore(fragilityVector);
    
    
    %fragilityVector is cstand
    
    rmsvars = zeros(1, cols);
    for i = 1:cols
        astand = datastand(:,i);
        wstand = astand\cstand;
        cstandpred = astand*wstand;
        
        corrvars(i) = linsolve(astand, cstand);
        rmsvars(i) = rms(cstand - cstandpred);
    end
    
    

    negativeCorrIndex = corrvars < 0;
    positiveCorrIndex = corrvars > 0;
    negativeIndices = find(negativeCorrIndex);
    positiveIndices = find(positiveCorrIndex);
    %Find the closest to 0 for rmsvars
    [~, minIndexNegative] = min(rmsvars(negativeIndices));
    [~, minIndexPositive] = min(rmsvars(positiveIndices));
    lowIndexNegative = negativeIndices(minIndexNegative);
    lowIndexPositive = positiveIndices(minIndexPositive);


end

function [rmstrain rmstest] = a2q2(lowndx)
% [RMSTRAIN RMSTEST]=A3Q2(LOWNDX) finds the RMS errors of 5-fold
% cross-validation for the variable LOWNDX of the data in the CSV file.
% The RMS errors for the training of each fold are returned
% in RMSTEST and the RMS errors for the testing of each fold are
% returned in RMSTEST.
%
% INPUTS:
%         LOWNDX   - integer scalar, index into the data
% OUTPUTS:
%         RMSTRAIN - 1x5 array of RMS errors for 5-fold training
%         RMSTEST  - 1x5 array of RMS errors for 5-fold testing

    % Read the test data from a CSV file and find the size of the data
    [fragilityVector,dataMatrix,countries,ageMatrix]=fragilitydata;
    [m n] = size(dataMatrix);
    
    A = dataMatrix;

    % Create Xmat and yvec from the data and the input parameter,
    % accounting for no standardization of data
    Xmat = A(:,[1:lowndx-1, lowndx+1:n]); %Other data as independent
    yvec = A(:,lowndx); %Select new age group as proxy of fragility index

    [rmstrain, rmstest] = mykfold(Xmat, yvec, 5);

end
function [rmstrain,rmstest]=mykfold(Xmat, yvec, k_in)
% [RMSTRAIN,RMSTEST]=MYKFOLD(XMAT,yvec,K) performs a k-fold validation
% of the least-squares linear fit of yvec to XMAT. If K is omitted,
% the default is 5.
%
% INPUTS:
%         XMAT     - MxN data vector
%         yvec     - Mx1 data vector
%         K        - positive integer, number of folds to use
% OUTPUTS:
%         RMSTRAIN - 1xK vector of RMS error of the training fits
%         RMSTEST  - 1xK vector of RMS error of the testing fits

    % Problem size
    M = size(Xmat, 1);

    % Set the number of folds; must be 1<k<M
    if nargin >= 3 & ~isempty(k_in)
        k = max(min(round(k_in), M-1), 2);
    else
        k = 5;
    end
    
    rng('default') %makes sure get same resuts each time (consistent seed)
    randidx = randperm(linspace(1,M,1)); %random permutation of rows
    
    Xmatperm = Xmat(randidx, :); %data
    yvecperm = yvec(randidx, :); %(new) dependent

    % Initialize the return variables
    rmstrain = zeros(1, k);
    rmstest  = zeros(1, k);

    % Determine the number of rows per experiment
    numInFold = floor(M/k);%Changed to floor because loop was going above 178
    %Orignially round

    % Process each fold
    for ix=1:k
        
        % %
        % % STUDENT CODE GOES HERE: replace the next 5 lines with code to
        % % (1) set up the "train" and "test" indexing for "xmat" and "yvec"
        % % (2) use the indexing to set up the "train" and "test" data
        % % (3) compute "wvec" for the training data
        % %
        %Start and end index for test set of fold
        teststart = (ix-1)*numInFold + 1; %Initial 1. currentteststart + numInFold
        testend = teststart + numInFold - 1; %currentteststart+numInFold-1
        
        
        %Set up train and test indices             
        testindices(ix,:) = teststart:testend; %Each row(ix) different fold
        trainidices(ix,:) = [1:teststart-1, testend+1:M]; %Train on non test data
        
        xmat_train  = Xmatperm(trainidices(ix,:),:);%Everything but first 36
        yvec_train = yvecperm(trainidices(ix,:));%Everything but first 36

        xmat_test = Xmatperm(testindices(ix,:),:);
        yvec_test = yvecperm(testindices(ix,:));

        %wvec for training
        wvec = xmat_train\yvec_train;
        

        rmstrain(ix) = rms(xmat_train*wvec - yvec_train);
        rmstest(ix)  = rms(xmat_test*wvec  - yvec_test);
    end

end
function [fragilityVector, dataMatrix, countries, ageMatrix] = ...
    fragilitydata
% [fragilityVector,dataMatrix,countries,ageMatrix]=fragilitydata
% loads and separates the 2013 fragility data for assessed countries
% N.B. These are proportions of male populations; other data
% are not loaded here. UN population estimates were used
%
% INPUTS:
%         none
% OUTPUTS:
%         fragilityVector - Mx1 vector of fragility "index" values
%         dataMatrix      - MxN matrix, M countries and N age groups
%         countries       - Mx1 cell array, strings for country names
%         ageMatrix       - 2xN matrix, start end end of each age group

    % Load the table from a known CSV file
    mt = readtable('fragility2013male.csv');

    % Extract the values of the fragility "index" for each country
    fragilityVector = table2array(mt(:,2));

    % Extract and normalize data for male populations
    dataRaw = table2array(mt(:,3:end));
    dataMatrix = dataRaw./sum(dataRaw, 2);

    % Extract the country names
    countries = table2array(mt(:,1));

    allNames = mt.Properties.VariableNames;
    ages = allNames(3:end);
    N = length(ages);

    ageMatrix = zeros(2, N);
    for jx = 1:N
        ageString = cell2mat(ages(jx));
        lowStart = strfind(ageString, 'm');
        lowEnd = strfind(ageString, '_');
        ageMatrix(1, jx) = str2num(ageString((lowStart+1):(lowEnd-1)));
        endValue = str2num(ageString((lowEnd+1):end));
        if ~isempty(endValue)
            ageMatrix(2, jx) = endValue;
        else
            ageMatrix(2, jx) = inf;
        end
    end
end