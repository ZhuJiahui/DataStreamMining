function topicevolution

%挖掘话题
%  根据聚类中心挖掘话题
%  进行简要的话题演化分析

%  B506
%  Computer Science School
%  Wuhan University, Wuhan 430072 China
%  zhujiahui@whu.edu.cn
%  2013-12-31

%% 读写文件目录 %%
tic;

% 重构的对象为金字塔压缩后的数据
% 重构按片进行，依据为时间
read_directory1 = 'dataset/cluster/topics_data23/final_cluster_center';
read_directory2 = 'dataset/cluster/topics_data23/final_each_cluster_number';
date_filename = 'dataset/cluster/topics_data23/new_date.txt';

write_filename = 'dataset/cluster/topics_data23/evolution.xlsx';

relation_matrix = -100 * ones(100, 100);

all_date = load(date_filename);
previous_center_data = load(strcat(read_directory1, '/1.txt'));
previous_support = load(strcat(read_directory2, '/1.txt'));

t_gamma = 13.5;
t_delta = 12;
t_eps = 0.3;

for i = 1 : length(all_date)

    if (i == 1)
        %%%%% 无演化 %%%%%
        fprintf('The first day!');
        % 全部算作出现
        for j = 1 : 2
            relation_matrix(j, j) = 0;
        end
      
    else
        % 不是第一天
        center_data = load(strcat(strcat(read_directory1, '/'), strcat(num2str(i), '.txt')));
        this_support = load(strcat(strcat(read_directory2, '/'), strcat(num2str(i), '.txt')));
        
        %%%%% 有演化 %%%%%
        KL_matrix1 = zeros(size(previous_center_data, 2), 2);
        KL_matrix2 = zeros(2, size(previous_center_data, 2));
        for k = 1 : size(previous_center_data, 2)
            for l = 1 : 2
                % 非对称的KL距离
                KL_matrix1(k, l) = KL_divergence(center_data(:, l), previous_center_data(:, k));                
            end
        end
        
        for k = 1 : 2
            for l = 1 : size(previous_center_data, 2)
                KL_matrix2(k, l) = KL_divergence(previous_center_data(:, l), center_data(:, k));
            end
        end
        
        %%%%%% 出现 %%%%%
        % 全部算作出现
        for j = 1 : 2
            relation_matrix(size(previous_center_data, 2) + j, size(previous_center_data, 2) + j) = 0;
        end
        
        %for k = 1 : 3
            %if min(KL_matrix1(:, k)) > t_gamma
                % 出现用0表示
                %relation_matrix(size(previous_center_data, 2) + k, size(previous_center_data, 2) + k) = 0;
            %else
                % 不出现用-2表示
                %relation_matrix(size(previous_center_data, 2) + k, size(previous_center_data, 2) + k) = -2;
            %end
        %end
        
        %%%%%% 持续 %%%%%
        for k = 1 : 2
            for l = 1 : length(KL_matrix1(:, k))
                if KL_matrix1(l, k) < t_delta && (this_support(k) / previous_support(l) < (1 + t_eps)) && (this_support(k) / previous_support(l) > (1 - t_eps)) && (relation_matrix(l, l) >= 0)
                    % 持续用3表示
                    relation_matrix(l, size(previous_center_data, 2) + k) = 3;
                    relation_matrix(size(previous_center_data, 2) + k, l) = 3;
                    %relation_matrix(size(previous_center_data, 2) + k, size(previous_center_data, 2) + k) = 0;
                end
            end
        end
        
        %%%%%% 合并 %%%%%
        for k = 1 : 2     
            indexes = [];
            for l = 1 : length(KL_matrix1(:, k))
                if KL_matrix1(l, k) > t_delta && KL_matrix1(l, k) < t_gamma && (relation_matrix(l, l) >= 0)
                    indexes = [indexes, l];
                end
            end
            if length(indexes) >= 2
                f1 = zeros(size(center_data, 1), 1);
                f2 = 0;
                for l = 1 : length(indexes)
                    f1 = f1 + previous_center_data(:, indexes(l)) * previous_support(indexes(l));
                    f2 = f2 + previous_support(indexes(l));
                end
                Z = f1 / f2;
                
                if KL_divergence(center_data(:, k), Z) <= t_delta - 0.23
                    for s = 1 : length(indexes)
                        relation_matrix(indexes(s), size(previous_center_data, 2) + k) = 1;
                        relation_matrix(size(previous_center_data, 2) + k, indexes(s)) = 1;
                        % relation_matrix(indexes(s), indexes(s)) = 1;
                        %relation_matrix(size(previous_center_data, 2) + k, size(previous_center_data, 2) + k) = 0;
                    end
                end
            end
        end
            
        %%%%%% 分割 %%%%%
        for k = 1 : size(previous_center_data, 2)
            indexes2 = [];
            for l = 1 : 2
                if KL_matrix1(k, l) > t_delta && KL_matrix1(k, l) < t_gamma && (relation_matrix(k, k) >= 0)
                    indexes2 = [indexes2, l];
                end
            end
            if length(indexes2) >= 2
                ff1 = zeros(size(center_data, 1), 1);
                ff2 = 0;
                for l = 1 : length(indexes2)
                    ff1 = ff1 + center_data(:, indexes2(l)) * this_support(indexes2(l));
                    ff2 = ff2 + this_support(indexes2(l));
                end
                ZZ = ff1 / ff2;
                
                if KL_divergence(ZZ, previous_center_data(:, k)) <= t_delta
                    for s = 1 : length(indexes2)
                        relation_matrix(k, size(previous_center_data, 2) + indexes2(s)) = 2;
                        relation_matrix(size(previous_center_data, 2) + indexes2(s), k) = 2;
                        %relation_matrix(size(previous_center_data, 2) + k, size(previous_center_data, 2) + k) = 0;
                    end
                end
            end
        end
        
        %%%%%% 消失 %%%%%
        for k = 1 : size(previous_center_data, 2)
            if min(KL_matrix2(:, k)) > t_gamma               
                relation_matrix(k, k) = -1;
            end
        end
        
        previous_center_data = [previous_center_data, center_data];
        previous_support = [previous_support', this_support']';
    end
    
    if i == 5
        break;
    end
end

relation_matrix = relation_matrix(1 : size(previous_center_data, 2), 1 : size(previous_center_data, 2));

xlswrite(write_filename, relation_matrix);

time = toc;
fprintf('用时%f秒\n', time);

end