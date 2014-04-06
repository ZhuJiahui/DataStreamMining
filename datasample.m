function datasample
%% 读写文件目录 %%

read_directory1 = 'D:/Local/workspace/MicroblogDataStreamCompress/dataset/batch_data_segment/topics_data2/update_vsm';
%read_directory1 = 'D:/Local/workspace/MicroblogDataStreamCompress/dataset/batch_data_segment/Music/Music2/update_vsm';

write_directory = 'dataset/non_orthogonal/topics_data2';
write_directory2 = 'dataset/non_orthogonal/topics_data2/采样信号';
%write_directory = 'dataset/non_orthogonal/Music2';
%write_directory2 = 'dataset/non_orthogonal/Music2/采样信号';

write_filename1 = 'dataset/non_orthogonal/topics_data2/Q.txt';
%write_filename1 = 'dataset/non_orthogonal/Music2/Q.txt';

if ~isdir(write_directory)
    mkdir(write_directory);
end
if ~isdir(write_directory2)
    mkdir(write_directory2);
end

data_files = dir(fullfile(read_directory1, '*.txt'));
sample_time = zeros(length(data_files), 1);

tic;
% 随便读取当前目录下的一个文件
data_first = load(strcat(read_directory1, '/1.txt'));

M_sample = 120;  % 中文微博采样维数
%M_sample = 250;  % 亚马逊采样维数
L = size(data_first, 2);  % 原始数据维度
Q = normrnd(1, sqrt(1 / sqrt(M_sample)), M_sample, L);  % 高斯随机矩阵
fprintf('测量矩阵的维度：%d * %d\n', size(Q, 1), size(Q, 2));
base_time = toc;

dlmwrite(write_filename1, full(Q), ' ');  % 写入文件

for i = 1 : length(data_files)
    data = load(strcat(strcat(read_directory1, '/'), strcat(num2str(i), '.txt')));
    fprintf('正在处理第%d片数据\n', i);
    
    tic;
    data1 = Q * data';
    this_time = toc;
    sample_time(i, 1) = base_time + this_time;
    
    fprintf('数据采样后的维度：%d * %d\n', size(data1, 1), size(data1, 2));
    
    % 若写入文件有错误，或数据是稀疏的，则使用full(data1)
    dlmwrite(strcat(strcat(write_directory2, '/'), strcat(num2str(i), '.txt')), data1, ' ');
    fprintf('第%d片数据处理完毕\n\n', i);
end

write_filename = strcat(write_directory, '/采样时间.xlsx');
xlswrite(write_filename, sample_time);

end
