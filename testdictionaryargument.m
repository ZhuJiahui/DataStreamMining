function testdictionaryargument
%获取通用字典
%  从全局数据中训练字典
%  字典训练采用KSVD训练方法
%  字典原子的选择采用Batch-OMP算法
%  训练所得的字典为非正交字典
%  详见ksvd

%  B506
%  Computer Science School
%  Wuhan University, Wuhan 430072 China
%  zhujiahui@whu.edu.cn
%  2013-12-08

disp(' ');
disp('  **********  获取通用字典  **********');
disp(' ');

%% 读写文件目录 %%

%read_filename = 'dataset/general_ksvd/topics_data2/训练数据.txt';
read_filename = 'dataset/general_ksvd/Music5/训练数据.txt';

%write_directory = 'dataset/non_orthogonal/topics_data2/800';
write_directory = 'dataset/non_orthogonal/Music5/3500';

if ~isdir(write_directory)
    mkdir(write_directory);
end

% 训练数据
global_data = load(read_filename);
printf('数据读取完毕\n');

tic;
%% KSVD训练过程 %%
params.data = global_data;  % 训练数据
%params.Tdata = 20;  % 中文微博数据稀疏度
params.Tdata = 25;  % 英文亚马逊数据稀疏度
params.dictsize = 3500;  % 字典原子个数
params.iternum = 30;  % KSVD训练迭代次数
params.memusage = 'high';  % 内存使用为高级状态较好

printf('正在训练，请等待\n');

% KSVD字典训练
[Dksvd, g, err] = ksvd(params, '');

time1 = toc;
printf('\n训练字典完毕，共用时%f秒\n', time1);

%% 结果展示 %%
final_error = err(length(err));  % 最终收敛残差
printf('\n均方误差：%f\n', final_error);
figure;
plot(err);
title('K-SVD error convergence');
xlabel('Iteration'); ylabel('RMSE');

printf('\n字典大小 %d * %d', size(global_data, 1), params.dictsize);
printf('数据个数: %d\n', size(global_data, 2));

tic;
%% 写入文件 %%

write_filename1 = strcat(write_directory, '/字典.txt');
dlmwrite(write_filename1, Dksvd, ' ');

write_filename2 = strcat(write_directory, '/稀疏系数.txt');
dlmwrite_each_line(write_filename2, g);

write_filename3 = strcat(write_directory, '/RMSE.xlsx');
xlswrite(write_filename3, err);

time2 = toc;
printf('写入文件用时%f秒', time2);


end

function dlmwrite_each_line(write_filename, X)
row = size(X, 1);
for i = 1 : row
    dlmwrite(write_filename, full(X(i, :)), '-append', 'delimiter', ' ')
end
end