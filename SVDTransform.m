function SVDTransform
%% 矩阵的奇异值分解，然后重构

%readFilename = 'dataset/TopNWordsTFIDFVSM/第一夫人彭丽媛.txt';
readFilename = 'dataset/TopNWordsTFIDFVSM/王菲薛蛮子.txt';
X = load(readFilename);  % 读取txt数据
disp('Data Load OK!');
selectKeys = 50;

[U, S, V] = svds(X, selectKeys);
clear X;
%S = S(:,1:selectKeys);
%V = V(:,1:selectKeys);
X2 = U * S * V';
X2 = X2';
col = size(X2, 2);  % 原始矩阵列数，即数据项的个数;
row = size(X2, 1);  % 原始矩阵行数，即数据的维度;
%writeFilename = 'dataset/TopNWordsTFIDFVSM/第一夫人彭丽媛2.txt';
writeFilename = 'dataset/TopNWordsTFIDFVSM/王菲薛蛮子SVD.txt';
fid = fopen(writeFilename, 'w+');
for i = 1:row
    for j = 1:col
        fprintf(fid, '%f', X2(i,j));
        if j~=col
            fprintf(fid, ' ');
        end
    end
    fprintf(fid, '\n');
end
fclose(fid); 
