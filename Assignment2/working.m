[rmsvars lowIndexPositive lowIndexNegative] = a2q1;

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
    ainterp = [dataMatrix() ones([size(dataMatrix,1), 1])]

    lowIndexPositive = min(rmsvars);
    lowIndexNegative = min(rmsvars);
    
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