function samplermse
%% 读写文件目录 %%
tic;
read_directory1 = 'D:/Local/workspace/MicroblogDataStreamCompress/dataset/batch_data_segment/topics_data1/update_vsm';
dictionary_filename = 'dataset/non_orthogonal/topics_data1/字典.txt';

thresh = 20;

write_directory = 'dataset/non_orthogonal/topics_data1/error';
write_directory2 = 'dataset/non_orthogonal/topics_data1/重构数据';

if ~isdir(write_directory)
    mkdir(write_directory);
end

if ~isdir(write_directory2)
    mkdir(write_directory2);
end

data_files = dir(fullfile(read_directory1, '*.txt'));
D = load(dictionary_filename);

error_matrix = zeros(length(data_files), 1);

for i = 1 : length(data_files)
    data = load(strcat(strcat(read_directory1, '/'), strcat(num2str(i), '.txt')));
    printf('正在处理第%d片数据', i);
    
    M_sample = 120;  % 采样维数
    L = size(data, 2);  % 原始数据维度
    Q = normrnd(0, sqrt(1 / sqrt(M_sample)), M_sample, L);
    printf('测量矩阵的维度：%d * %d', size(Q, 1), size(Q, 2));
    
    data1 = Q * data';
    printf('数据采样后的维度：%d * %d',size(data1, 1), size(data1, 2));

    data = data';
    
    D1 = Q * D;
    Gamma = OMP(D1, data1, thresh);    % 这个OMP比较慢
    data2 = (round(D * Gamma))';
    dlmwrite(strcat(strcat(write_directory2, '/'), strcat(num2str(i), '.txt')), full(data2), ' ');

    error_matrix(i, 1) = compute_err(D, Gamma, data);
    
    printf('第%d片数据处理完毕\n', i);
end

write_filename = strcat(write_directory, '/采样字典重构误差.xlsx');
xlswrite(write_filename, error_matrix);

toc;
end



%% 计算残差 %%
function err = compute_err(D, Gamma, data)
% 基于稀疏度限制的残差计算
err = sqrt(sum(reperror2(data, D, Gamma)) / numel(data));
end


%% 分块计算残差的平方和
function err2 = reperror2(X, D, Gamma)

err2 = zeros(1, size(X, 2));
blocksize = 2000;
for i = 1 : blocksize : size(X, 2)
    blockids = i : min(i + blocksize - 1, size(X, 2));
    err2(blockids) = sum((X(:, blockids) - round(D * Gamma(:, blockids))) .^ 2);
end
end


