function ENsamplermse
%% 读写文件目录 %%
tic;
read_directory1 = 'D:/Local/workspace/MicroblogDataStreamCompress/dataset/batch_data_segment/Music/Music2/update_vsm';
dictionary_filename = 'dataset/non_orthogonal/Music2/字典.txt';

thresh = 25;

write_directory = 'dataset/non_orthogonal/Music2/error';

if ~isdir(write_directory)
    mkdir(write_directory);
end

data_files = dir(fullfile(read_directory1, '*.txt'));
D = load(dictionary_filename);

Q_filename = 'dataset/non_orthogonal/Music2/Q.txt';
Q = load(Q_filename);
D1 = Q * D;

error_matrix = zeros(length(data_files), 1);

for i = 1 : length(data_files)
    data = load(strcat(strcat(read_directory1, '/'), strcat(num2str(i), '.txt')));
    fprintf('\n正在处理第%d片数据\n', i);
   
    data1 = Q * data';

    data = data';
    
    Gamma = OMP(D1, data1, thresh);    % 这两个OMP比较慢
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


