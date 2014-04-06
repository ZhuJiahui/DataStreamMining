function elementcompare

read_filename1 = 'dataset/batch_data_segment/topics_data1/update_vsm/1.txt';
read_filename2 = 'dataset/general_ksvd/topics_data1/error/XXX.xlsx';
dictionary_filename = 'dataset/general_ksvd/topics_data1/字典.txt';

thresh = 10;

D = load(dictionary_filename);
data = load(read_filename1);

M_sample = 50;  % 采样维数
L = size(data, 2);  % 原始数据维度
Q = sqrt(M_sample) * normrnd(0,1,M_sample, L);
printf('测量矩阵的维度：%d * %d', size(Q, 1), size(Q, 2));

data1 = Q * data';
data = data';

D1 = Q * D;
Gamma = OMP(D1, data1, thresh);    % 这两个OMP效果是一样的，但是这个慢2个数量级
%Gamma = omp(D1' * data1, D1' * D1, thresh, ompparams{:});
this_error = compute_err(D, Gamma, data);

printf('处理完毕，残差%f\n', this_error);


%% 比较
X1 = xlsread(read_filename2);
compare_matrix = zeros(size(data, 1), size(data, 2));
count = 0;
for i = 1 : size(data, 1)
    for j = i : size(data, 2)
        if data(i, j) ~= X1(i, j)
            compare_matrix(i, j) = 1;
            count = count + 1;
        end
    end
end
printf('错误元素个数：%d', count);
xlswrite('dataset/general_ksvd/topics_data1/error/1.xlsx', data);
xlswrite('dataset/general_ksvd/topics_data1/error/compare.xlsx', compare_matrix);

end

%% 计算残差 %%
function err = compute_err(D, Gamma, data)
% 基于稀疏度限制的残差计算
err = sqrt(sum(reperror2(data, D, Gamma)) / numel(data));
%err = sum(reperror2(data, D, Gamma)) / numel(data);
end


%% 分块计算残差的平方和
function err2 = reperror2(X, D, Gamma)

err2 = zeros(1, size(X, 2));
XXX = zeros(size(X, 1), size(X,2));

for i = 1 : size(X, 2)
    element_X = D * Gamma(:, i);
    for j = 1 : length(element_X)
        if element_X(j) <= 0
            element_X(j) = 0;
        else
            element_X(j) = round(element_X(j));
        end
    end
        
    XXX(:, i) = element_X;
    err2(i) = sum((X(:, i) - element_X) .^ 2);
end
xlswrite('dataset/general_ksvd/topics_data1/error/XXX.xlsx', XXX);
end
