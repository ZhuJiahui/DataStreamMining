function postcluster

%将粗粒度的聚类中心按天数变为细粒度
%  即将属于同一天但不属于同一片内的聚类中心何为一体
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
read_directory1 = 'dataset/cluster/topics_data22/merge_cluster_center';
read_directory2 = 'dataset/cluster/topics_data22/each_cluster_number';
date_filename = 'dataset/cluster/topics_data22/date.txt';

write_directory1 = 'dataset/cluster/topics_data22/final_cluster_center';
write_directory2 = 'dataset/cluster/topics_data22/final_each_cluster_number';
write_date_filename = 'dataset/cluster/topics_data22/new_date.txt';

if ~isdir(write_directory1)
    mkdir(write_directory1);
end

if ~isdir(write_directory2)
    mkdir(write_directory2);
end

all_date = load(date_filename);
previous_date = all_date(1);
previous_center_data = load(strcat(read_directory1, '/1.txt'));
previous_center_data = previous_center_data';
previous_support = load(strcat(read_directory2, '/1.txt'));

file_number = 1;
fid_d = fopen(write_date_filename, 'w+');

for i = 1 : length(all_date)
    
    this_date = all_date(i);
    this_center_data = load(strcat(strcat(read_directory1, '/'), strcat(num2str(i), '.txt')));
    this_center_data = this_center_data';
    this_support = load(strcat(strcat(read_directory2, '/'), strcat(num2str(i), '.txt')));
    
    if (i == 1)
        fprintf('The first Day!\n');
        
    elseif this_date - previous_date >= 1
        % 不是同一天
        
        % 聚类
        cluster_number = 2;
        ev_number = 6;
        if size(previous_center_data, 2) <= 2
            % 当前天的聚类数据（聚类中心）写入文件
            % 写入后每一列代表一条信息
            dlmwrite(strcat(strcat(write_directory1, '/'), strcat(num2str(file_number), '.txt')), previous_center_data, ' ');
            dlmwrite(strcat(strcat(write_directory2, '/'), strcat(num2str(file_number), '.txt')), previous_support, ' ');
            file_number = file_number + 1;
            % 将时间信息等写入文件，此处只有一个日期
            fprintf(fid_d, '%s\n', num2str(previous_date));
        else
            [cluster_tag, center, sum_to_center, each_to_center] = spectral_cluster(previous_center_data', cluster_number, ev_number);
            % 聚类分析
            center_data = zeros(size(previous_center_data, 1), cluster_number);

            for k = 1 : size(each_to_center, 2)
                [min_value, min_index] = min(each_to_center(:, k));
                center_data(:, k) = previous_center_data(:, min_index);
            end
            
            each_cluster_number = zeros(cluster_number, 1);
            for k = 1 : length(cluster_tag)
                for l = 1 : cluster_number
                    if cluster_tag(k) == l
                        each_cluster_number(l, 1) = each_cluster_number(l, 1) + previous_support(k);
                    end
                end
            end
            
            % 当前天的聚类数据（聚类中心）写入文件
            % 写入后每一列代表一条信息
            dlmwrite(strcat(strcat(write_directory1, '/'), strcat(num2str(file_number), '.txt')), center_data, ' ');
            dlmwrite(strcat(strcat(write_directory2, '/'), strcat(num2str(file_number), '.txt')), each_cluster_number, ' ');
            file_number = file_number + 1;
            % 将时间信息等写入文件，此处只有一个日期
            fprintf(fid_d, '%s\n', num2str(previous_date));
            
        end
        
        previous_date = all_date(i);
        previous_center_data = this_center_data;
        previous_support = this_support;
        
    else
        % 是同一天
        fprintf('\n与上面是同一天\n');
        
        previous_center_data = [previous_center_data, this_center_data];
        previous_support = [previous_support', this_support']';
        
    end
    
end


%% 最后部分处理
% 聚类
cluster_number = 2;
ev_number = 6;
if size(previous_center_data, 2) <= 2
    % 当前天的聚类数据（聚类中心）写入文件
    % 写入后每一列代表一条信息
    dlmwrite(strcat(strcat(write_directory1, '/'), strcat(num2str(file_number), '.txt')), previous_center_data, ' ');
    dlmwrite(strcat(strcat(write_directory2, '/'), strcat(num2str(file_number), '.txt')), previous_support, ' ');

    % 将时间信息等写入文件，此处只有一个日期
    fprintf(fid_d, '%s\n', num2str(previous_date));
else
    [cluster_tag, center, sum_to_center, each_to_center] = spectral_cluster(previous_center_data', cluster_number, ev_number);
    % 聚类分析
    center_data = zeros(size(previous_center_data, 1), cluster_number);
    for k = 1 : size(each_to_center, 2)
        [min_value, min_index] = min(each_to_center(:, k));
        center_data(:, k) = previous_center_data(:, min_index);
    end
    
    each_cluster_number = zeros(cluster_number, 1);
    for k = 1 : length(cluster_tag)
        for l = 1 : cluster_number
            if cluster_tag(k) == l
                each_cluster_number(l, 1) = each_cluster_number(l, 1) + previous_support(k);
            end
        end
    end
    
    % 当前天的聚类数据（聚类中心）写入文件
    % 写入后每一列代表一条信息
    dlmwrite(strcat(strcat(write_directory1, '/'), strcat(num2str(file_number), '.txt')), center_data, ' ');
    dlmwrite(strcat(strcat(write_directory2, '/'), strcat(num2str(file_number), '.txt')), each_cluster_number, ' ');
    % 将时间信息等写入文件，此处只有一个日期
    fprintf(fid_d, '%s\n', num2str(previous_date));
    
end

fclose(fid_d);


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
