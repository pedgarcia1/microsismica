function [cluster_labels, outlier_labels] = kmeans_clustering_with_outliers(data, num_clusters, k_lof, lof_threshold,stage_i)
% Performs k-means clustering on a matrix of size Nx3, where N is the number
% of points in 3D space. Columns are X, Y, and Z coordinates.
%
% Inputs:
%   - data: Nx3 matrix of data points
%   - num_clusters: number of clusters to identify
%   - k_lof: number of neighbors to use for LOF algorithm
%   - lof_threshold: threshold for LOF algorithm (higher values result in fewer outliers)
%
% Outputs:
%   - cluster_labels: N-length vector of cluster assignments (values 1 to num_clusters)
%   - outlier_labels: N-length logical vector identifying outliers (true = outlier)

% Perform k-means clustering
cluster_labels = kmeans(data, num_clusters);

% Compute LOF scores
[lof_scores] = LOF(data, k_lof);

% Identify outliers based on LOF threshold
outlier_labels = lof_scores > lof_threshold;

% Plot results with outliers in red
colors = lines(num_clusters);
figure('Name',sprintf('Clustering: Stage %d',stage_i));
hold on;
for i = 1:num_clusters
    scatter3(data(cluster_labels==i,1), data(cluster_labels==i,2), data(cluster_labels==i,3), 36, colors(i,:), 'filled');
    lgdsCell{i} = sprintf("Cluster %d",i);
end
scatter3(data(outlier_labels,1), data(outlier_labels,2), data(outlier_labels,3), 36, 'r', 'filled', 'MarkerEdgeColor', 'k');
lgdsCell{end+1} = 'Outliers'; legend(lgdsCell);
view(-45,30);
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');

end

