function ENreconstruct

%金字塔压缩数据重构
%  数据先用OMP算法重构
%  按金字塔层来重构

%  B506
%  Computer Science School
%  Wuhan University, Wuhan 430072 China
%  zhujiahui@whu.edu.cn
%  2013-12-31

%% 读写文件目录 %%
tic;
read_directory1 = 'dataset/pyramid2/Music5/data';

Q_filename = 'dataset/non_orthogonal/Music5/Q.txt';
dictionary_filename = 'dataset/non_orthogonal/Music5/字典.txt';

write_directory = 'dataset/pyramid2/Music5/reconstruct_data';

if ~isdir(write_directory)
    mkdir(write_directory);
end

Q = load(Q_filename);
D = load(dictionary_filename);
D1 = Q * D;

% 稀疏度25
thresh = 25;

data_files = dir(fullfile(read_directory1, '*.txt'));

for i = 1 : length(data_files)
    data = load(strcat(strcat(read_directory1, '/'), strcat(num2str(i), '.txt')));
    printf('正在处理第%d片数据', i);
    
    % 当前数据片重构
    S = OMP(D1, data, thresh);
    % 每一列代表一条数据
    reconstruct_data = round(D * S);
    
    dlmwrite(strcat(strcat(write_directory, '/'), strcat(num2str(i), '.txt')), full(reconstruct_data), ' ');
    printf('第%d片数据处理完毕\n', i);
end

time = toc;
fprintf('用时%f秒\n', time);

end
