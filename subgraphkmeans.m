clear,clc
tic;
readFilename = 'dataset/KSVDProcess3/王菲薛蛮子/稀疏系数1200.xlsx';
S = xlsread(readFilename);
S=S';
L=size(S,2);% 原始矩阵列数，即数据项的个数;
n=size(S,1);% 原始矩阵行数，即字典原子个数;
disp('finish the file reading!!!');
toc;

tic;
%% 计算node之间的相似度矩阵
%tabulate(S_index(:));  %%% 各个原子的使用频率相差不大
node_matrix = zeros(n, n);
N_matrix = zeros(n, n);  %顶点度矩阵

for i = 1 : n
    for j = i : n
        node_matrix(i, j) = pdist2(S(i, :), S(j, :), 'Euclidean');
        node_matrix(j, i) = node_matrix(i, j);
    end
    N_matrix(i,i)= sum(node_matrix(i,:));
end

disp('finish the node similarity computing!!!');
toc;
tic;
%% 基于相似度矩阵的NJW谱聚类
L_matrix = N_matrix - node_matrix; %构建拉普拉斯矩阵
for i=1:n
    N_matrix(i,i)= N_matrix(i,i)^(-1/2);
end
L_matrix=N_matrix*L_matrix*N_matrix; %拉普拉斯矩阵规范化
[D,V]=eig(L_matrix);
k=3; % 取特征值的个数
new_D=D(:,1:10);

for i=1:n
    for j=1:10
        nor_new_D(i,j)=new_D(i,j) / (sum(new_D(i,:).^2)^(0.5));
    end
end

% K-Means聚类
[cluster_tag, center, sum_to_center, each_to_center] = kmeans(nor_new_D, 2);

disp('finish the clustering!!!');
toc;

tic;
%% 找寻聚类中心在原始数据中的标注
min_distance = 10000000;
center_index = zeros(1, size(center, 1));
for i=1:size(center, 1)
    for j = 1:size(nor_new_D, 1)
        distance_to_center = pdist2(center(i, :), nor_new_D(j, :), 'Euclidean');
        if distance_to_center < min_distance
            min_distance = distance_to_center;
            center_index(i) = j;
        end
    end
    min_distance = 10000000;
end
disp('Find the center position!');

readFilename1 = 'dataset/ClassTag/王菲薛蛮子.txt';
tag = load(readFilename1);
exact_tag = tag(12001 : 13200);
center_tag = tag(12000 + center_index);
if (center_tag(2) == 2 && center_tag(1) == 2)
    center_tag(2) = 1;
else
    if center_tag(2) == 1 && center_tag(1) == 1
        center_tag(2) = 2;
    end
end

hold on;

right_tag = 0;
for i = 1 : length(cluster_tag)
    cluster_tag(i) = center_tag(cluster_tag(i));
    if cluster_tag(i) == exact_tag(i)
        right_tag = right_tag + 1;
    end
end

purity = right_tag / 1200;
printf('聚类纯度：%f', purity);
toc;

figure(1);
hold on;
xlim([-1, 1201]);
ylim([0, 4]);
for k = 1 : length(cluster_tag)
    if (cluster_tag(k) == 1)
        plot(k, cluster_tag(k), 'r*', 'MarkerSize', 3);
    else
        if (cluster_tag(k) == 2)
            plot(k, cluster_tag(k), 'go', 'MarkerSize', 3);
        else
            plot(k, cluster_tag(k), 'bp', 'MarkerSize', 3);
        end
    end
end
hold off;

figure(2);
%line(L_matrix(index1,1),L_matrix(index1,2),'marker','*','color','g');
%line(L_matrix(index2,1),L_matrix(index2,2),'marker','*','color','r');
plot([center([1 2],1)],[center([1 2],2)],'*','color','k');

