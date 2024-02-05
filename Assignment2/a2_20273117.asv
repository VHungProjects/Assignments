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
 
    [rmsvars lowIndexPositive lowIndexNegative] = a2q1;
    [rmstrain rmstest] = a2q2(lowIndexPositive);

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
    % %
    % % STUDENT CODE GOES HERE: REMOVE THE NEXT 3 LINES AND THIS COMMENT
    % % THEN PERFORM THE COMPUTATIONS
    % %
    %One Vector for this data size
    onesVec = ones(size(fragilityVector,1), 1);
    %Loop through each age group
    for i = 1:cols
        atest = dataMatrix(:,1:col ~= i);
        wtest = atest\fragilityVector;

        rmsvars = (fragilityVector - atest*wtest);
    end
    %Find lowest RMS error of fit
    [val,lowIndex] = min(rmsvars);

    lowIndexPositive = 1;
    lowIndexNegative = cols;
    
    % Find the regression on your choice of standardized
    % or unstandardized variables
    % %
    % % STUDENT CODE GOES HERE: REMOVE THIS COMMENT
    % % THEN PERFORM THE COMPUTATIONS
    % %

    % % SAMPLE PLOT: REMOVE THIS COMMENT AND THE NEXT 14 LINES OF CODE
    % % THE "WFIT" VALUES ARE MADE UP FOR THIS EXAMPLE
    % % XVALS ARE DEDUCED FROM THE PLOT; YVALS ARE THE LINEAR FIT
    ageLo = ageMatrix(1, 1); % Where the age range starts
    ageHi = ageMatrix(2, 1);% Where the age range ends
    ageString = sprintf('Males in range %d:%d', ageLo, ageHi);
    plot(dataMatrix(:,1), fragilityVector, '.', 'MarkerSize', 10);
    xlabel('Proportion of Male Population');
    ylabel('Fragility Index');
    title(ageString);
    wfit = [400, 30]; % Slope/intercept are made up for this data
    axisVector = axis();
    xVals = axisVector(1:2); % X for the left and right sides of the plot
    yVals = wfit(1)*xVals + wfit(2); % Slope/intercept computation
    hold on;
    plot(xVals, yVals, 'k-');
    hold off;
    

    % Plot the results
    % %
    % % STUDENT CODE GOES HERE: REMOVE THIS COMMENT
    % % THEN PLOT THE RESULTS
    % %

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


    % Create Xmat and yvec from the data and the input parameter,
    % accounting for no standardization of data
    % %
    % % STUDENT CODE GOES HERE: REMOVE THIS COMMENT
    % % THEN ASSIGN THE VARIABLES FROM THE DATASET
    % %

    % Compute the RMS errors of 5-fold cross-validation
    % %
    % % STUDENT CODE GOES HERE: REMOVE THE NEXT 2 LINES AND THIS COMMENT
    % % THEN PERFORM THE COMPUTATIONS
    % %
    rmstrain = 0.5*ones(1,5);
    rmstest =  0.6*ones(1,5);

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

    % Initialize the return variables
    rmstrain = zeros(1, k);
    rmstest  = zeros(1, k);

    % Process each fold
    for ix=1:k
        % %
        % % STUDENT CODE GOES HERE: replace the next 5 lines with code to
        % % (1) set up the "train" and "test" indexing for "xmat" and "yvec"
        % % (2) use the indexing to set up the "train" and "test" data
        % % (3) compute "wvec" for the training data
        % %
        xmat_train  = [0 1];
        yvec_train  = 0;
        wvec = [0 0];
        xmat_test = [0 1];
        yvec_test = 0;

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