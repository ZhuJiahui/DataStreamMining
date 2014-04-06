function ENpyramidmft

%将采样后的英文亚马逊数据按金字塔形进行压缩
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
read_directory1 = 'D:/Local/workspace/MicroblogDataStreamCompress/dataset/spur/Music5/compress';
read_directory2 = 'D:/Local/workspace/MicroblogDataStreamCompress/dataset/batch_data_segment/Music/Music5/entropy';
read_directory3 = 'D:/Local/workspace/MicroblogDataStreamCompress/dataset/batch_data_segment/Music/Music5/update_phst';

write_directory1 = 'dataset/dspurmft/Music5/60/data';
write_directory2 = 'dataset/dspurmft/Music5/60/phst';
write_directory3 = 'dataset/dspurmft/Music5/60/original_index';
write_directory4 = 'dataset/dspurmft/Music5/60/entropy';

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
%data_files = dir(fullfile(read_directory1, '*.txt'));
file_number = 60;
%original_count = 0;
original_count = 9196110;
compress_count = 0;

% 产生金字塔压缩索引
if file_number < 2
    compress_index = cell({1});
    fprintf('数据片数 < 2，无需压缩');
else
    compress_index = pyramid_index(file_number);
    
    % 设置每一层窗口的大小
    gram = 9000;  % 英文亚马逊数据为2 * 4500

    % 循环遍历金字塔索引的每一层
    for i = 1 : length(compress_index)
        
        fprintf('正在压缩第%d层数据\n', i);
        % 根据压缩的规律选择该层的每个数据片的信息条数
        each_select = floor(gram / length(compress_index{i}));
        rest_select = gram - each_select * (length(compress_index{i}) - 1);
        %select = gram / length(compress_index{i});    
        
        % 每一层的所有数据片压缩后的数据
        data_level_result = cell(gram, 1);
        
        % 压缩的时间信息等
        p_level_result = cell(gram, 1);
        h_level_result = cell(gram, 1);
        s_level_result = cell(gram, 1);
        t_level_result = cell(gram, 1);
        
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

            %注意，读入的数据不规则，是一个cell
            [spur_data, this_count] = get_text_to_cell(strcat(strcat(read_directory1, '/'), strcat(num2str(compress_index{i}(j)), '.txt'))); 
            %original_count = original_count + this_count;
 
            % 注意此处有字符串数据
            fid = fopen(strcat(strcat(read_directory3, '/'), strcat(num2str(compress_index{i}(j)), '.txt')));
            phst = textscan(fid, '%s %s %s %s');
            fclose(fid);
            
            % 赋值给压缩结果
            for k = 1 : length(update_index)
                data_level_result{column_count} = spur_data{update_index(k)};
                this_number = 0;
                for l = 1 : length(data_level_result{column_count})
                    if data_level_result{column_count}{l}(1) == 'z'
                        this_number = this_number + 1;
                    end
                end
                compress_count = compress_count + this_number;
                
                p_level_result(column_count) = phst{1, 1}(update_index(k));
                h_level_result(column_count) = phst{1, 2}(update_index(k));
                s_level_result(column_count) = phst{1, 3}(update_index(k));
                t_level_result(column_count) = phst{1, 4}(update_index(k));
                
                entropy_level_result(column_count) = entropy(update_index(k), 1);
                
                column_count = column_count + 1;
            end
            
        end
        
        % 一层遍历完后，将该层的所有数据片压缩后的数据写入文件
        % 写入后每一列代表一条信息（未压缩的采样数据也是每一列代表一条信息）
        dlmwrite_cell2(strcat(strcat(write_directory1, '/'), strcat(num2str(i), '.txt')), data_level_result);
        
        % 将该层的所有数据的时间信息等写入文件
        write_phst_to_text(strcat(strcat(write_directory2, '/'), strcat(num2str(i), '.txt')),  p_level_result,  h_level_result,  s_level_result,  t_level_result);
        % 将该层的所有数据的熵值信息等写入文件
        dlmwrite(strcat(strcat(write_directory4, '/'), strcat(num2str(i), '.txt')), entropy_level_result, ' ');

        fprintf('第%d层数据压缩完毕\n', i);
    end
end

fprintf('所有数据压缩完毕\n');
fprintf('压缩总量： %d\n', original_count);

ratio = compress_count / original_count;
fprintf('当前的压缩占总数的比重： %f\n', ratio);

time = toc;
fprintf('用时%f秒\n', time);

% 金字塔压缩索引
write_filename = 'dataset/dspurmft/Music5/60/pyramid_index.txt';
dlmwrite_cell(write_filename, compress_index);


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
    dlmwrite(write_filename, X{i}, '-append', 'delimiter', ' ')
end
end

%% 压缩数据写入文件
function dlmwrite_cell2(write_filename, X)

row = length(X);

fid = fopen(write_filename, 'w+');
for i = 1 : row
    for j = 1 : length(X{i})
        fprintf(fid, '%s', X{i}{j});
        if j ~= length(X{i})
            fprintf(fid, ' ');
        end
    end
    fprintf(fid, '\n');
end

fclose(fid); 
end

%% 读取文件数据到元胞矩阵 %%
function [data_cell, total_count]= get_text_to_cell(read_filename)
total_count = 0;
fid = fopen(read_filename);
i = 1;
while ~feof(fid)
    line_cell{i} = fgetl(fid);
    data_cell{i} = regexp(line_cell{i}, '\s+', 'split');
    total_count = total_count + length(data_cell{i});
    i = i + 1;
end
fclose(fid); 
end

%% 字符串元胞矩阵写入文件 %%
function write_phst_to_text(write_filename, p, h, s, t)

row = length(p);

fid = fopen(write_filename, 'w+');
for i = 1 : row
    fprintf(fid, '%s %s %s %s\n', p{i}, h{i}, s{i}, t{i});
end

fclose(fid); 
end