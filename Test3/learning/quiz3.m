iris = load('fisheriris.mat');
data = iris.meas;
yvec = iris.species;
data_std = zscore(data);
[u, s, v] = svd(data_std);
r = rank(data_std);
plot(sum(s)/sum(sum(s)));

%% 

species = categorical(yvec);
species = renamecats(species,{'setosa','versicolor','virginica'},{'1','2','3'});
species = str2double(string(species));


x = data(:, 3:4);
clusters = kmeans(x, 3, 'start', [1.5 0.3; 4.2 1.3; 5.9 2.1]);

correct = sum((clusters == species) == 1);
disp(correct / size(clusters, 1))




zerox = x - mean(x);
zerox = zscore(zerox);
[U, S, V] = svd(zerox);
Xmat = zerox*V(:,1);
[n, ~] = size(Xmat);
yvec = zeros(n, 1);
for ix = 51:100
    yvec(ix) = 1;
end
for ix = 101:150
    yvec(ix) = 2;
end
gscatter(Xmat ,yvec ,species)