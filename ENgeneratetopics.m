function ENgeneratetopics

%挖掘话题
%  根据时间选择数据
%  数据先用OMP算法重构
%  再在重构的空间内聚类
%  聚类的方法采用谱聚类
%  根据聚类中心挖掘话题

%  B506
%  Computer Science School
%  Wuhan University, Wuhan 430072 China
%  zhujiahui@whu.edu.cn
%  2013-12-31

date_start = '2013-02-15';
date_end = '2013-02-25';

% 压缩数据中时间戳对应关系 1970-01-01 0:00:00 = -28800
start_unix_time = (datenum(date_start) - datenum('1970-01-01')) * 86400;
end_unix_time = (datenum(date_end) - datenum('1970-01-01')) * 86400;

%% 读写文件目录 %%
tic;

% 重构的对象为金字塔压缩后的数据
% 重构按片进行，依据为时间
read_directory1 = 'dataset/pyramid/Music5/data';
read_directory2 = 'dataset/pyramid/Music5/phst';
read_directory3 = 'dataset/pyramid/Music5/original_index';

Q_filename = 'dataset/non_orthogonal/Music5/Q.txt';
dictionary_filename = 'dataset/non_orthogonal/Music5/字典.txt';

write_directory1 = 'dataset/cluster/Music5/cluster_center';
write_directory2 = 'dataset/cluster/Music5/each_cluster_number';

write_date_filename = 'dataset/cluster/Music5/date.txt';
batch_index_filename = 'dataset/cluster/Music5/batch_index.txt';

if ~isdir(write_directory1)
    mkdir(write_directory1);
end
if ~isdir(write_directory2)
    mkdir(write_directory2);
end


Q = load(Q_filename);
D = load(dictionary_filename);
D1 = Q * D;

% 稀疏度25
thresh = 25;

% 可能的待聚类的数据总片数
data_files = dir(fullfile(read_directory1, '*.txt'));
file_number = 1;

fid_d = fopen(write_date_filename, 'w+');

for i = 1 : length(data_files)
    % 注意此处有字符串数据
    fid = fopen(strcat(strcat(read_directory2, '/'), strcat(num2str(i), '.txt')));
    phst = textscan(fid, '%s %s %s %f');
    fclose(fid);
    
    if phst{1, 4}(end) < start_unix_time
        fprintf('Skip\n');
    elseif phst{1, 4}(1) > end_unix_time
        break
    else
        % 压缩的金字塔数据按列存储，每一列代表一条数据
        pyramid_data = load(strcat(strcat(read_directory1, '/'), strcat(num2str(i), '.txt')));
        % 金字塔压缩的索引，2列，第一列表示原始数据片编号，第二列表示在相应数据片中的位置，二者的编号都是从1开始
        original_index = load(strcat(strcat(read_directory3, '/'), strcat(num2str(i), '.txt')));
        
        candidate_j = intersect(find(phst{1, 4} >= start_unix_time), find(phst{1, 4} <= end_unix_time));
        
        j = candidate_j(1);
        previous_batch_id = original_index(j, 1);
        previous_time = start_unix_time;
        batch_data = [];
        
        while j <= candidate_j(end)
            this_batch_id = original_index(j, 1);
            this_unix_time = phst{1, 4}(j);
            
            %% while内第一层if %%
            % 既要满足在原始数据中是同一片，又要满足是同一天
            if this_batch_id == previous_batch_id && (this_unix_time - previous_time) < 86400
                batch_data = [batch_data, pyramid_data(:, j)];
                
                % 若已至末尾，则聚类
                if j == candidate_j(end)
                    %%%%% 聚类分析 %%%%%
                    % 当前数据片重构
                    S = OMP(D1, batch_data, thresh);
                    % 每一列代表一条数据
                    reconstruct_data = round(D * S);
                
                    % 聚类
                    cluster_number = 3;
                    ev_number = 5;
                    if size(reconstruct_data, 2) <= ev_number
                        fprintf('\n数据太少，无法聚类！\n');
                    else                
                        [cluster_tag, center, sum_to_center, each_to_center] = spectral_cluster(reconstruct_data', cluster_number, ev_number);
                
                        % 聚类分析
                        center_data = zeros(size(reconstruct_data, 1), cluster_number);
                        for k = 1 : size(each_to_center, 2)
                            [min_value, min_index] = min(each_to_center(:, k));
                            center_data(:, k) = reconstruct_data(:, min_index);
                        end
                    
                        each_cluster_number = zeros(cluster_number, 1);
                        for k = 1 : length(cluster_tag)
                            for l = 1 : cluster_number
                                if cluster_tag(k) == l
                                    each_cluster_number(l, 1) = each_cluster_number(l, 1) + 1;
                                end
                            end
                        end
                                
                
                        % 当前天的聚类数据（聚类中心）写入文件
                        % 写入后每一列代表一条信息
                        dlmwrite(strcat(strcat(write_directory1, '/'), strcat(num2str(file_number), '.txt')), center_data, ' ');
                        dlmwrite(strcat(strcat(write_directory2, '/'), strcat(num2str(file_number), '.txt')), each_cluster_number, ' ');
                        file_number = file_number + 1;
                        % 将时间信息等写入文件，此处只有一个日期
                        fprintf(fid_d, '%s\n', num2str(previous_time));
                        dlmwrite(batch_index_filename, this_batch_id, '-append', 'delimiter', ' ');                       
                    end
                    %%%%% 聚类分析结束 %%%%%
                end
                
            %% while内第一层else %%
            else
                %%%%% 聚类分析 %%%%%
                % 当前数据片重构
                S = OMP(D1, batch_data, thresh);
                % 每一列代表一条数据
                reconstruct_data = round(D * S);
                
                % 聚类
                cluster_number = 3;
                ev_number = 5;
                if size(reconstruct_data, 2) <= ev_number
                    fprintf('\n数据太少，无法聚类！\n');
                else                
                    [cluster_tag, center, sum_to_center, each_to_center] = spectral_cluster(reconstruct_data', cluster_number, ev_number);
                
                    % 聚类分析
                    center_data = zeros(size(reconstruct_data, 1), cluster_number);
                    for k = 1 : size(each_to_center, 2)
                        [min_value, min_index] = min(each_to_center(:, k));
                        center_data(:, k) = reconstruct_data(:, min_index);
                    end
                    
                    each_cluster_number = zeros(cluster_number, 1);
                    for k = 1 : length(cluster_tag)
                        for l = 1 : cluster_number
                            if cluster_tag(k) == l
                                each_cluster_number(l, 1) = each_cluster_number(l, 1) + 1;
                            end
                        end
                    end
                                
                
                    % 当前天的聚类数据（聚类中心）写入文件
                    % 写入后每一列代表一条信息
                    dlmwrite(strcat(strcat(write_directory1, '/'), strcat(num2str(file_number), '.txt')), center_data, ' ');
                    dlmwrite(strcat(strcat(write_directory2, '/'), strcat(num2str(file_number), '.txt')), each_cluster_number, ' ');
                    file_number = file_number + 1;
                    % 将时间信息等写入文件，此处只有一个日期
                    fprintf(fid_d, '%s\n', num2str(previous_time));
                    dlmwrite(batch_index_filename, previous_batch_id, '-append', 'delimiter', ' ');
                end
                %%%%% 聚类分析结束 %%%%%
                
                if this_batch_id == previous_batch_id && (this_unix_time - previous_time) >= 86400
                    previous_time = previous_time + 86400;
                elseif this_batch_id ~= previous_batch_id && (this_unix_time - previous_time) < 86400
                    previous_batch_id = this_batch_id;
                else
                    previous_time = previous_time + 86400;
                    previous_batch_id = this_batch_id;
                end
                
                batch_data = pyramid_data(:, j);
            end
            %% while内第一层end %%
            
            j = j + 1;
        end
    end
end

fclose(fid_d);

fprintf('\n所有数据聚类完毕\n');
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

    % K-Means聚类
    % 注意K-Means聚类输入的数据按行来
    [cluster_tag, center, sum_to_center, each_to_center] = kmeans(norm_select_E_vectors, cluster_number);
    disp('finish the clustering!!!');
end

end
