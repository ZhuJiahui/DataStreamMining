function compressquery

%在金字塔形压缩数据进行查询

%  B506
%  Computer Science School
%  Wuhan University, Wuhan 430072 China
%  zhujiahui@whu.edu.cn
%  2013-12-14

%% 读写文件目录 %%
tic;

% 压缩的对象为信号和时间信息等
% 查询的依据为时间和关键字
read_directory1 = 'dataset/pyramid/topics_data1/data';
read_directory2 = 'dataset/pyramid/topics_data1/id_time';
read_directory3 = 'dataset/pyramid/topics_data1/original_index';

write_directory1 = 'dataset/pyramid/topics_data1/data';
write_directory2 = 'dataset/pyramid/topics_data1/id_time';
write_directory3 = 'dataset/pyramid/topics_data1/original_index';

if ~isdir(write_directory1)
    mkdir(write_directory1);
end
if ~isdir(write_directory2)
    mkdir(write_directory2);
end
if ~isdir(write_directory3)
    mkdir(write_directory3);
end

keyword_list = {'李', '天一'};
mode = 'AND';
time_interval = {'2013-7-1', '2013-7-5'};
time_interval2 = query_time_format(time_interval);
select = 30;

write_filename = write_directory + u'/q2.txt'

read_filename1 = 'dataset/non_orthogonal/topics_data1/Q.txt';
Q = load(read_filename1);

% 待压缩的数据总片数
data_files = dir(fullfile(read_directory2, '*.txt'));

% 选定待重构的数据片
% 读取时间信息
for i = 1 : length(data_files)
    id_time = load(strcat(strcat(read_directory2, '/'), strcat(num2str(i), '.txt')));
    time_lines = id_time(:, 2);

    % 当前片的最晚时间比查询设定的开始时间还早，则跳过该片
    if time_lines(end) < time_interval2(1)
         continue;
    %当前片的最早时间比查询设定的结束时间晚，则结束
    elseif time_lines(1) > time_interval2(2)
          break;
    else
        % 压缩后的数据项对应原始数据的索引
        % 每一个数据片中逐行遍历
        for j = 1 : length(time_lines)
            % 当前遍历时的时间
                now_t = float(time_lines[j].strip().split()[-1])
                
                if (now_t >= start) and (now_t <= end) and (j in data_index):
                    if mode == "OR":
                        flag = 0
                        for each1 in keyword_list:
                            for k in range(len(word_list)):
                                if (each1 in word_list[k]) and (float(each_weibo_vsm[j][k]) > 0.000001):
                                    this_message = " ".join(vsm_map_word(each_weibo_vsm[j], word_list))
                                    if this_message not in query_result:
                                        query_result.append(this_message)
                                        entropy_result.append(entropy_list[j])
                                    flag = 1
                                    break
                            if flag == 1:
                                break
                    else:
                        flag = 0
                        for each1 in keyword_list:
                            for k in range(len(word_list)):
                                if (each1 in word_list[k]) and (float(each_weibo_vsm[j][k]) > 0.000001):
                                    flag += 1
                                    break
                                
                        if flag == len(keyword_list):
                            this_message = " ".join(vsm_map_word(each_weibo_vsm[j], word_list))
                            if this_message not in query_result:
                                query_result.append(this_message)
                                entropy_result.append(entropy_list[j])
            

% 产生金字塔压缩索引
if length(data_files) < 2
    compress_index = cell({1});
    fprintf('数据片数 < 2，无需压缩');
else
    compress_index = pyramid_index(length(data_files));
    
    % 设置每一层窗口的大小
    gram = 4096;  % 中文微博数据为2 * 2048

    % 循环遍历金字塔索引的每一层
    for i = 1 : length(compress_index)
        
        fprintf('正在压缩第%d层数据\n', i);
        % 根据压缩的规律选择该层的每个数据片的信息条数
        select = gram / length(compress_index{i});    
        
        % 每一层的所有数据片压缩后的数据
        data_level_result = zeros(120, gram);
        
        % 压缩的时间信息等
        id_level_result = cell(gram, 1);
        time_level_result = cell(gram, 1);

        column_count = 1;
        
        % 遍历每一层的每一个数据片编号
        for j = 1 : length(compress_index{i})
            
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
                
                column_count = column_count + 1;
            end
            
        end
        
        % 一层遍历完后，将该层的所有数据片压缩后的数据写入文件
        % 写入后每一列代表一条信息（未压缩的采样数据也是每一列代表一条信息）
        dlmwrite(strcat(strcat(write_directory1, '/'), strcat(num2str(i), '.txt')), full(data_level_result), ' ');
         % 将该层的所有数据的时间信息等写入文件
        write_it_to_text(strcat(strcat(write_directory2, '/'), strcat(num2str(i), '.txt')),  id_level_result,  time_level_result);

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
function result_time = query_time_format(time_interval)

%2013-01-01----41275
    
result_time = zeros(1, 2);
base_time = datenum('2013-1-1');
base_number = 41274;

for i = 1 : length(time_interval)
    first_time = datenum(time_interval{i});
    result_time{i} = base_number + first_time - base_time;
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