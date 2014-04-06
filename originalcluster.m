function originalcluster

%以2天为单位进行聚类
%  聚类的方式采样谱聚类

%  B506
%  Computer Science School
%  Wuhan University, Wuhan 430072 China
%  zhujiahui@whu.edu.cn
%  2013-12-31

%% 读写文件目录 %%
tic;

% 重构的对象为金字塔压缩后的数据
% 重构按片进行，依据为时间
read_directory1 = 'dataset/cluster/topics_data22/original_merge_vsm';

write_directory1 = 'dataset/cluster/o2/original_cluster_center';
write_directory2 = 'dataset/cluster/o2/original_each_cluster_number';
write_directory3 = 'dataset/cluster/o2/original_cluster_tag';

if ~isdir(write_directory1)
    mkdir(write_directory1);
end

if ~isdir(write_directory2)
    mkdir(write_directory2);
end

if ~isdir(write_directory3)
    mkdir(write_directory3);
end

% 待聚类的数据总片数
data_files = dir(fullfile(read_directory1, '*.txt'));

cluster_number = 2;
ev_number = 6;

for i  = 1 : length(data_files)
    
    fprintf('正在处理第%d片数据\n', i);
    
    data = load(strcat(strcat(read_directory1, '/'), strcat(num2str(i), '.txt')));
    data_length = size(data, 1);
    
    each_batch = 500;
    batch_number = floor(data_length / each_batch);
    row = 1;
    
    merge_center_data = zeros(size(data, 2), cluster_number * batch_number);
    merge_cluster_number = zeros(batch_number * cluster_number, 1);
    merge_cluster_tag = zeros(each_batch * batch_number, 1);
    
    center_count = 0;
    
    % 分布聚类
    for j = 1 : batch_number
        data2 = data(row : (row + each_batch - 1), :);
        
        % 聚类
        [cluster_tag, center, sum_to_center, each_to_center] = spectral_cluster(data2, cluster_number, ev_number);
        
        % 聚类分析
        for k = 1 : size(each_to_center, 2)
            [min_value, min_index] = min(each_to_center(:, k));
            merge_center_data(:, center_count + k) = data(min_index, :);
        end
        
        each_cluster_number = zeros(cluster_number, 1);
        for k = 1 : length(cluster_tag)
            for l = 1 : cluster_number
                if cluster_tag(k) == l
                    each_cluster_number(l, 1) = each_cluster_number(l, 1) + 1;
                end
            end
        end
        
        
        for k = 1 : cluster_number
            merge_cluster_number((2 * j + k - 2), 1) = each_cluster_number(k, 1);
        end
        
        merge_cluster_tag(row : (row + each_batch -1), 1) = cluster_tag;
        
        center_count = center_count + cluster_number;
        row = row + each_batch;
    end
    
    
    
    % 归总聚类
    [cluster_tag, center, sum_to_center, each_to_center] = spectral_cluster(merge_center_data', cluster_number, ev_number);
    
    % 聚类分析
    center_data = zeros(size(merge_center_data, 1), cluster_number);
    
    for k = 1 : size(each_to_center, 2)
        [min_value, min_index] = min(each_to_center(:, k));
        center_data(:, k) = data(min_index, :);
    end
    
    final_each_cluster_number = zeros(cluster_number, 1);
    for k = 1 : length(cluster_tag)
        for l = 1 : cluster_number
            if cluster_tag(k) == l
                final_each_cluster_number(l, 1) = final_each_cluster_number(l, 1) + merge_cluster_number(k, 1);
            end
        end
        
        % 只聚2个类，可枚举
        if mod(k, 2) == 1
            for kk = 1 : each_batch
                part = (k + 1) / 2;
                if merge_cluster_tag((part - 1) * each_batch + kk) == 1
                    merge_cluster_tag((part - 1) * each_batch + kk) = cluster_tag(k);
                end
            end
        else
            for kk = 1 : each_batch
                part = k / 2;
                if merge_cluster_tag((part - 1) * each_batch + kk) == 2
                    merge_cluster_tag((part - 1) * each_batch + kk) = cluster_tag(k);
                end
            end
        end
    end
    
    % 当前天的聚类数据（聚类中心）写入文件
    % 写入后每一列代表一条信息
    dlmwrite(strcat(strcat(write_directory1, '/'), strcat(num2str(i), '.txt')), center_data, ' ');
    dlmwrite(strcat(strcat(write_directory2, '/'), strcat(num2str(i), '.txt')), final_each_cluster_number, ' ');
    dlmwrite(strcat(strcat(write_directory3, '/'), strcat(num2str(i), '.txt')), merge_cluster_tag, ' ');
    
    fprintf('第%d片数据聚类完毕\n', i);
    
end

fprintf('所有数据聚类完毕\n');

time = toc;
fprintf('用时%f秒\n', time);

end


%% 谱聚类
function [cluster_tag, center, sum_to_center, each_to_center] = spectral_cluster(data, cluster_number, ev_number)

% 计算node之间的相似度矩阵
n = size(data, 1);  % 行数代表数据个数
node_matrix = zeros(n, n);
degree_matrix = zeros(n, n);  %顶点度矩阵

for i = 1 : n
    for j = i : n
        d1 = KL_divergence(data(i, :), data(j, :));
        d2 = KL_divergence(data(j, :), data(i, :));
        node_matrix(i, j) = max(d1, d2);
        % node_matrix(i, j) = pdist2(data(i, :), data(j, :), 'Euclidean');
        node_matrix(j, i) = node_matrix(i, j);
    end
    degree_matrix(i, i) = sum(node_matrix(i, :));
end

disp('finish the node similarity computing!!!');

% 基于相似度矩阵的NJW谱聚类
L_matrix = degree_matrix - node_matrix;  % 构建拉普拉斯矩阵
for i = 1 : n
    degree_matrix(i, i) = degree_matrix(i, i) ^ (-1 / 2);
end

L_matrix = degree_matrix * L_matrix * degree_matrix;  % 拉普拉斯矩阵规范化

[E_vectors, E_values] = eig(L_matrix);
k = ev_number;  % 取特征向量的个数
if k > n
    fpfintf('\n数据个数太少！\n');
else
    
    select_E_vectors = E_vectors(:, 1 : k);
    norm_select_E_vectors = zeros(n, k);
    
    % 按行单位化
    for i = 1 : n
        for j = 1 : k
            norm_select_E_vectors(i, j) = select_E_vectors(i, j) / (sum(select_E_vectors(i, :) .^ 2) ^ (0.5));
        end
    end
    
    %if size(norm_select_E_vectors, 1) >= cluster_number
    
    % K-Means聚类
    % 注意K-Means聚类输入的数据按行来
    [cluster_tag, center, sum_to_center, each_to_center] = kmeans(norm_select_E_vectors, cluster_number, 'emptyaction','singleton');
    disp('finish the clustering!!!');
end

end
