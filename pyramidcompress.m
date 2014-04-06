function pyramidcompress

%将采样后的微博数据按金字塔形进行压缩
%  压缩的对象为信号和时间信息等
%  保留数据在原始数据中的索引信息

%  B506
%  Computer Science School
%  Wuhan University, Wuhan 430072 China
%  zhujiahui@whu.edu.cn
%  2013-12-14

%% 读写文件目录 %%
tic;

% 压缩的对象为信号和时间信息等
% 压缩的依据为信息熵值
read_directory1 = 'dataset/non_orthogonal/topics_data1/采样信号';
read_directory2 = 'D:/Local/workspace/MicroblogDataStreamCompress/dataset/batch_data_segment/topics_data1/entropy';
read_directory3 = 'D:/Local/workspace/MicroblogDataStreamCompress/dataset/batch_data_segment/topics_data1/update_id_time';

write_directory1 = 'dataset/pyramid/topics_data1/data';
write_directory2 = 'dataset/pyramid/topics_data1/id_time';
write_directory3 = 'dataset/pyramid/topics_data1/original_index';
write_directory4 = 'dataset/pyramid/topics_data1/entropy';

if ~isdir(write_directory1)
    mkdir(write_directory1);
end
if ~isdir(write_directory2)
    mkdir(write_directory2);
end
if ~isdir(write_directory3)
    mkdir(write_directory3);
end
if ~isdir(write_directory4)
    mkdir(write_directory4);
end

% 待压缩的数据总片数
data_files = dir(fullfile(read_directory1, '*.txt'));
%data_files = 100;

% 产生金字塔压缩索引
if length(data_files) < 2
    compress_index = cell({1});
    fprintf('数据片数 < 2，无需压缩');
else
    compress_index = pyramid_index(length(data_files));
    %compress_index = pyramid_index(data_files);
    
    % 设置每一层窗口的大小
    gram = 5000;  % 中文微博数据为2 * 2500

    % 循环遍历金字塔索引的每一层
    for i = 1 : length(compress_index)
        
        fprintf('正在压缩第%d层数据\n', i);
        % 根据压缩的规律选择该层的每个数据片的信息条数
        each_select = floor(gram / length(compress_index{i}));
        rest_select = gram - each_select * (length(compress_index{i}) - 1);  
        
        % 每一层的所有数据片压缩后的数据
        data_level_result = zeros(120, gram);
        
        % 压缩的时间信息等
        id_level_result = cell(gram, 1);
        time_level_result = cell(gram, 1);
        
        % 每一层的所有数据熵值压缩后的数据
        entropy_level_result = zeros(gram, 1);

        column_count = 1;
        
        % 遍历每一层的每一个数据片编号
        for j = 1 : length(compress_index{i})
            if j == length(compress_index{i})
                select = rest_select;
            else
                select = each_select;
            end
            
            % 读取熵值和行数索引编号，作为待排序的数据项
            entropy = load(strcat(strcat(read_directory2, '/'), strcat(num2str(compress_index{i}(j)), '.txt')));
            line_index = 1 : length(entropy);
            line_index = line_index';
            
            % 按熵值降序排序
            el = [entropy, line_index];
            el1 = sortrows(el, -1);
            
            % 选择前select项行号索引
            % 进行升序排序，确保压缩后的信息的相对位置不变
            update_index = el1(:, 2);
            update_index = sort(update_index(1 : select));
            
            % 压缩后保留的数据在原先数据中的索引位置
            mix_index = [(compress_index{i}(j) * ones(select, 1)), update_index];
            dlmwrite(strcat(strcat(write_directory3, '/'), strcat(num2str(i), '.txt')), mix_index, '-append', 'delimiter', ' ');
            
            %注意，读入的采样数据每一列代表一条信息
            sample_data = load(strcat(strcat(read_directory1, '/'), strcat(num2str(compress_index{i}(j)), '.txt')));
            
            % 注意此处有字符串数据
            fid = fopen(strcat(strcat(read_directory3, '/'), strcat(num2str(compress_index{i}(j)), '.txt')));
            phst = textscan(fid, '%s %s');
            fclose(fid);
            
            % 赋值给压缩结果
            for k = 1 : length(update_index)
                data_level_result(:, column_count) = sample_data(:, update_index(k));
                
                id_level_result(column_count) = phst{1, 1}(update_index(k));
                time_level_result(column_count) = phst{1, 2}(update_index(k));
                
                entropy_level_result(column_count) = entropy(update_index(k), 1);
                
                column_count = column_count + 1;
            end
            
        end
        
        % 一层遍历完后，将该层的所有数据片压缩后的数据写入文件
        % 写入后每一列代表一条信息（未压缩的采样数据也是每一列代表一条信息）
        dlmwrite(strcat(strcat(write_directory1, '/'), strcat(num2str(i), '.txt')), full(data_level_result), ' ');
        % 将该层的所有数据的时间信息等写入文件
        write_it_to_text(strcat(strcat(write_directory2, '/'), strcat(num2str(i), '.txt')),  id_level_result,  time_level_result);
        % 将该层的所有数据的熵值信息等写入文件
        dlmwrite(strcat(strcat(write_directory4, '/'), strcat(num2str(i), '.txt')), entropy_level_result, ' ');

        fprintf('第%d层数据压缩完毕\n', i);
    end
end

fprintf('所有数据压缩完毕\n');

% 金字塔压缩索引
write_filename = 'dataset/pyramid/topics_data1/pyramid_index.txt';
dlmwrite_cell(write_filename, compress_index);

toc;
end


%% 产生合并压缩的索引
function pyramid_all = pyramid_index(batch_number)

% batch_number:文件总个数
% 金字塔初始设置
level_element = 2;
pyramid_all = cell({[1,2]});

%以2片数据为单位进行压缩
for i = 4 : 2 : batch_number
    flag = 0;
    level = length(pyramid_all);
    
    for j = level - 1 : -1 : 1
        if (length(pyramid_all{j}) < level_element(j))
            %若某一层未满，则合并压缩
            pyramid_all{j} = [pyramid_all{j}, pyramid_all{j + 1}];
            pyramid_all(j + 1) = [];
            pyramid_all{end + 1} = [i - 1, i];
            flag = 1;
            break;
        end
    end
    
    %若所有层均满，则新建一层
    if flag == 0
        pyramid_all{end + 1} = [i - 1, i];
        level_element = [2 ^ (level + 1), level_element];
    end
end
end


%% 金字塔压缩索引写入文件
function dlmwrite_cell(write_filename, X)
row = length(X);
for i = 1 : row
    dlmwrite(write_filename, full(X{i}), '-append', 'delimiter', ' ')
end
end


%% 字符串元胞矩阵写入文件 %%
function write_it_to_text(write_filename, id, time)

row = length(id);

fid = fopen(write_filename, 'w+');
for i = 1 : row
    fprintf(fid, '%s %s\n', id{i}, time{i});
end

fclose(fid); 
end