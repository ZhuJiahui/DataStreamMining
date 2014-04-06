function gettrainingdata
%% 读写文件目录 %%
tic;
read_directory = 'D:/Local/workspace/MicroblogDataStreamCompress/dataset/batch_data_segment/Music/Music2/update_vsm';
%read_directory = 'D:/Local/workspace/MicroblogDataStreamCompress/dataset/batch_data_segment/topics_data1/update_vsm';

%write_directory = 'dataset/general_ksvd/topics_data1';
write_directory = 'dataset/general_ksvd/Music2';

if ~isdir(write_directory)
    mkdir(write_directory);
end

vsm_files = dir(fullfile(read_directory, '*txt'));

%% 采样准备 %%
primitive_data = load(strcat(strcat(read_directory, '/'), vsm_files(1).name));
data_demension = size(primitive_data, 2);
select_data_number = 10000;
global_data = zeros(data_demension, select_data_number);
global_index = 1;

each_data_number = floor(select_data_number / length(vsm_files));
important_index = round(length(vsm_files) / 2);
rest_data_number = select_data_number - each_data_number * (length(vsm_files) - 1);

%each_data_number = floor(select_data_number / (length(vsm_files) / 5));
%important_index = 4 * length(vsm_files) / 5 + round(length(vsm_files) / 5 / 2);
%rest_data_number = select_data_number - each_data_number * (length(vsm_files) / 5 - 1);

%% 循环采样数据 %%
for i = 1 : length(vsm_files)
    % 此时行数为信号数目，列为数据项维度
    X = load(strcat(strcat(read_directory, '/'), strcat(num2str(i), '.txt')));
    
    L = size(X, 1);  % 原始矩阵行数，即数据项的个数;
    
    perm = randperm(L);
    if i == important_index
        perm = perm(1 : rest_data_number);
    else
        perm = perm(1 : each_data_number);
    end
    perm = sort(perm);
    printf('第%d片数据选取%d项数据', i, length(perm));
    for j = 1 : length(perm)
        global_data(:, global_index) = X(perm(j), :);
        global_index = global_index + 1;
    end
end
time1 = toc;
fprintf('\n数据采集完毕，用时%f秒\n', time1);
dlmwrite_each_line(strcat(write_directory, '/训练数据.txt'), global_data);
end

function dlmwrite_each_line(write_filename, X)
row = size(X, 1);
for i = 1 : row
    dlmwrite(write_filename, full(X(i, :)), '-append', 'delimiter', ' ')
end
end

