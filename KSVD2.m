function [D,Gamma,err,gerr] = KSVD2(params,varargin)

%%
%KSVD K-SVD dictionary training.
%  [D,GAMMA] = KSVD(PARAMS) runs the K-SVD dictionary training algorithm on
%  the specified set of signals, returning the trained dictionary D and the
%  signal representation matrix GAMMA.
%
%  KSVD has two modes of operation: sparsity-based and error-based. For
%  sparsity-based minimization, the optimization problem is given by
%
%      min  |X-D*GAMMA|_F^2      s.t.  |Gamma_i|_0 <= T
%    D,Gamma
%
%  where X is the set of training signals, Gamma_i is the i-th column of
%  Gamma, and T is the target sparsity. For error-based minimization, the
%  optimization problem is given by
%
%      min  |Gamma|_0      s.t.  |X_i - D*Gamma_i|_2 <= EPSILON
%    D,Gamma
%
%  where X_i is the i-th training signal, and EPSILON is the target error.
%
%  [D,GAMMA,ERR] = KSVD(PARAMS) also returns the target function values
%  after each algorithm iteration. For sparsity-constrained minimization,
%  the returned values are given by
%
%    ERR(D,GAMMA) = RMSE(X,D*GAMMA) = sqrt( |X-D*GAMMA|_F^2 / numel(X) ) .
%
%  For error-constrained minimization, the returned values are given by
%
%    ERR(D,GAMMA) = mean{ |Gamma_i|_0 } = |Gamma|_0 / size(X,2) .
%
%  Error computation slightly increases function runtime.
%
%  [D,GAMMA,ERR,GERR] = KSVD(PARAMS) computes the target function values on
%  the specified set of test signals as well, usually for the purpose of
%  validation (testing the generalization of the dictionary). This requires
%  that the field 'testdata' be present in PARAMS (see below). The length
%  of ERR and GERR is identical.
%
%  [...] = KSVD(...,VERBOSE) where VERBOSE is a character string, specifies
%  messages to be printed during the training iterations. VERBOSE should
%  contain one or more of the characters 'i', 'r' and 't', each of which
%  corresponds to a certain piece of information:
%
%    i - iteration number
%    r - number of replaced atoms
%    t - target function value (and its value on the test data if provided)
%
%  Specifying either 'r', 't' or both, also implies 'i' automatically. For
%  example, KSVD(PARAMS,'tr') prints the iteration number, number of
%  replaced atoms, and target function value, at the end of each iteration.
%  The default value for VERBOSE is 't'. Specifying VERBOSE='' invokes
%  silent mode, and cancels all messages.
%
%  [...] = KSVD(...,MSGDELTA) specifies additional messages to be printed
%  within each iteration. MSGDELTA should be a positive number representing
%  the interval in seconds between messages. A zero or negative value
%  indicates no such messages (default). Note that specifying VERBOSE=''
%  causes KSVD to run in silent mode, ignoring the value of MSGDELTA.
%
%
%  Required fields in PARAMS:
%  --------------------------
%
%    'data' - Training data.
%      A matrix containing the training signals as its columns.
%
%    'Tdata' / 'Edata' - Sparse coding target.
%      Specifies the number of coefficients (Tdata) or the target error in
%      L2-norm (Edata) for coding each signal. If only one is present, that
%      value is used. If both are present, Tdata is used, unless the field
%      'codemode' is specified (below).
%
%    'initdict' / 'dictsize' - Initial dictionary / no. of atoms to train.
%      At least one of these two should be present in PARAMS.
%
%      'dictsize' specifies the number of dictionary atoms to train. If it
%      is specified without the parameter 'initdict', the dictionary is
%      initialized with dictsize randomly selected training signals.
%
%      'initdict' specifies the initial dictionary for the training. It
%      should be either a matrix of size NxL, where N=size(data,1), or an
%      index vector of length L, specifying the indices of the examples to
%      use as initial atoms. If 'dictsize' and 'initdict' are both present,
%      L must be >= dictsize, and in this case the dictionary is
%      initialized using the first dictsize columns from initdict. If only
%      'initdict' is specified, dictsize is set to L.
%
%
%  Optional fields in PARAMS:
%  --------------------------
%
%    'testdata' - Validation data.
%      If present, specifies data on which to compute generalization error.
%      Should be a matrix containing the validation signals as its columns.
%
%    'iternum' - Number of training iterations.
%      Specifies the number of K-SVD iterations to perform. If not
%      specified, the default is 10.
%
%    'memusage' - Memory usage.
%      This parameter controls memory usage of the function. 'memusage'
%      should be one of the strings 'high', 'normal' (default) or 'low'.
%      When set to 'high', the fastest implementation of OMP is used, which
%      involves precomputing both G=D'*D and DtX=D'*X. This increasese
%      speed but also requires a significant amount of memory. When set to
%      'normal', only the matrix G is precomputed, which requires much less
%      memory but slightly decreases performance. Finally, when set to
%      'low', neither matrix is precomputed. This should only be used when
%      the trained dictionary is highly redundant and memory resources are
%      very low, as this will dramatically increase runtime. See function
%      OMP for more details.
%
%    'codemode' - Sparse-coding target mode.
%      Specifies whether the 'Tdata' or 'Edata' fields should be used for
%      the sparse-coding stopping criterion. This is useful when both
%      fields are present in PARAMS. 'codemode' should be one of the
%      strings 'sparsity' or 'error'. If it is not present, and both fields
%      are specified, sparsity-based coding takes place.
%
%    'exact' - Exact K-SVD update.
%      Specifies whether the exact or approximate dictionary update
%      should be used. By default, the approximate computation is used,
%      which is significantly faster and requires less memory. Specifying a
%      nonzero value for 'exact' causes the exact computation to be used
%      instead, which slows down the method but provides slightly improved
%      results. The exact update uses SVD to solve the rank-1 minimization
%      problem, while the approximate upate performs alternate-optimization
%      to solve this problem.
%
%
%  Optional fields in PARAMS - advanced:
%  -------------------------------------
%
%    'maxatoms' - Maximal number of atoms in signal representation.
%      When error-based sparse coding is used, this parameter can be used
%      to specify a hard limit on the number of atoms in each signal
%      representation (see parameter 'maxatoms' in OMP2 for more details).
%
%    'muthresh' - Mutual incoherence threshold.
%      This parameter can be used to control the mutual incoherence of the
%      trained dictionary, and is typically between 0.9 and 1. At the end
%      of each iteration, the trained dictionary is "cleaned" by discarding
%      atoms with correlation > muthresh. The default value for muthresh is
%      0.99. Specifying a value of 1 or higher cancels this type of
%      cleaning completely. Note: the trained dictionary is not guaranteed
%      to have a mutual incoherence less than muthresh. However, a method
%      to track this is using the VERBOSE parameter to print the number of
%      replaced atoms each iteration; when this number drops near zero, it
%      is more likely that the mutual incoherence of the dictionary is
%      below muthresh.
%
%
%   Summary of all fields in PARAMS:
%   --------------------------------
%
%   Required:
%     'data'                   training data
%     'Tdata' / 'Edata'        sparse-coding target
%     'initdict' / 'dictsize'  initial dictionary / dictionary size
%
%   Optional (default values in parentheses):
%     'testdata'               validation data (none)
%     'iternum'                number of training iterations (10)
%     'memusage'               'low, 'normal' or 'high' ('normal')
%     'codemode'               'sparsity' or 'error' ('sparsity')
%     'exact'                  exact update instead of approximate (0)
%     'maxatoms'               max # of atoms in error sparse-coding (none)
%     'muthresh'               mutual incoherence threshold (0.99)
%
%
%  References:
%  [1] M. Aharon, M. Elad, and A.M. Bruckstein, "The K-SVD: An Algorithm
%      for Designing of Overcomplete Dictionaries for Sparse
%      Representation", the IEEE Trans. On Signal Processing, Vol. 54, no.
%      11, pp. 4311-4322, November 2006.
%  [2] M. Elad, R. Rubinstein, and M. Zibulevsky, "Efficient Implementation
%      of the K-SVD Algorithm using Batch Orthogonal Matching Pursuit",
%      Technical Report - CS, Technion, April 2008.
%
%  See also KSVDDENOISE, OMPDENOISE, OMP, OMP2.


%%
%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  August 2009


%% 全局变量定义

% 迭代条件限制性模式，分基于稀疏度限制和基于残差限制2种
global CODE_SPARSITY CODE_ERROR codemode

% 内存使用状态，分low，normal，high3种
global MEM_LOW MEM_NORMAL MEM_HIGH memusage

% omp函数，omp参数，精确svd计算标识定义
global ompfunc ompparams exactsvd  

% 以下定义各种情况的对应标识
CODE_SPARSITY = 1;
CODE_ERROR = 2;

MEM_LOW = 1;
MEM_NORMAL = 2;
MEM_HIGH = 3;


%%
%%%%% 以下解析输入参数 %%%%%

%%
data = params.data;  % 原始信号数据
ompparams = {'checkdict', 'off'};

%% 迭代条件限制性模式 %%
if (isfield(params, 'codemode'))
    switch lower(params.codemode)
        case 'sparsity'
            codemode = CODE_SPARSITY;
            thresh = params.Tdata;
        case 'error'
            codemode = CODE_ERROR;
            thresh = params.Edata;
        otherwise
            error('Invalid coding mode specified');
    end
elseif (isfield(params, 'Tdata'))
    codemode = CODE_SPARSITY;
    thresh = params.Tdata;
elseif (isfield(params, 'Edata'))
    codemode = CODE_ERROR;
    thresh = params.Edata;
    
else
    error('Data sparse-coding target not specified');
end

if (codemode == CODE_ERROR && isfield(params, 'maxatoms'))
    ompparams{end + 1} = 'maxatoms';  % 元胞矩阵的end+1维赋值一个标签
    ompparams{end + 1} = params.maxatoms; % 用于表示单个信号的字典原子的最大个数
end


%% 内存使用 %%
% 分low，normal，high3种状态
if (isfield(params, 'memusage'))
    switch lower(params.memusage)
        case 'low'
            memusage = MEM_LOW;
        case 'normal'
            memusage = MEM_NORMAL;
        case 'high'
            memusage = MEM_HIGH;
        otherwise
            error('Invalid memory usage mode');
    end
else
    % 默认为normal
    memusage = MEM_NORMAL;
end


%% 迭代计数 %%
if (isfield(params, 'iternum'))
    iternum = params.iternum;
else
    % 默认迭代次数为10
    iternum = 10;
end


%% OMP函数 %%
if (codemode == CODE_SPARSITY)
    % 如果是基于稀疏度限制，那么用omp
    ompfunc = @omp;
else
    % 如果是基于残差限制，那么用omp2
    ompfunc = @omp2;
end


%% 状态消息标识 %%
printiter = 0;
printreplaced = 0;
printerr = 0;
printgerr = 0;

verbose = 't';
msgdelta = -1;  % 输出显示的的延迟停顿时间

for i = 1:length(varargin)
    if (ischar(varargin{i}))
        verbose = varargin{i};
    elseif (isnumeric(varargin{i}))
        msgdelta = varargin{i};
    else
        error('Invalid call syntax');
    end
end

for i = 1:length(verbose)
    switch lower(verbose(i))
        case 'i'
            printiter = 1;
        case 'r'
            printiter = 1;
            printreplaced = 1;
        case 't'
            printiter = 1;
            printerr = 1;
            if (isfield(params,'testdata'))
                printgerr = 1;
            end
    end
end

if (msgdelta<=0 || isempty(verbose))
    msgdelta = -1;
end

ompparams{end + 1} = 'messages';
ompparams{end + 1} = msgdelta;

% 当输出参数大于等于3，或者需要输出残差变化向量的时候，需要使计算残差标识为1
comperr = (nargout >= 3 || printerr);  % 计算残差标识


%% 验证标识 %%
testgen = 0;  % 默认无数据验证
if (isfield(params, 'testdata'))
    testdata = params.testdata;
    if (nargout >= 4 || printgerr)
        % 若输出有gerr，则需进行验证
        testgen = 1;
    end
end


%% 残差向量归一化 %%
XtX = []; 
XtXg = [];
if (codemode == CODE_ERROR && memusage == MEM_HIGH)
    % 如果迭代终止准则是基于残差限制的模式以及内存使用状态为high
    % 找出原始信号数据按列分块计算的平方和向量，向量维度为分块的数目
    XtX = colnorms_squared(data);
    if (testgen)
        % 若需数据验证，则计算用于验证的数据的平方和向量
        XtXg = colnorms_squared(testdata);
    end
end


%% 互相关系数阈值设定 %%
if (isfield(params, 'muthresh'))
    muthresh = params.muthresh;
else
    muthresh = 0.99;
end
if (muthresh < 0)
    error('invalid muthresh value, must be non-negative');
end


%% 精确svd计算 %
exactsvd = 0;  % 精确svd计算标识
if (isfield(params, 'exact') && params.exact ~= 0)
    exactsvd = 1;
end


%% 待训练字典原子数的确定 %%
if (isfield(params, 'initdict'))
    % any判断元素是否为0，非零返回1；all里面向量都为非0则返回逻辑值1.否则返回逻辑值0
    % 其实该语句就是判断矩阵是一维还是二维
    if (any(size(params.initdict) == 1) && all(iswhole(params.initdict(:))))
        % 如果给定了一维索引向量，索引向量的每一个值表示对应的原始信号数据的位置
        % 即用某几个原始信号数据作为待训练字典的原子总数
        % 那么待训练字典的原子数就跟索引向量的长度一样
        dictsize = length(params.initdict);  
    else
        % 其他情况下使用给定的初始字典的列数作为训练字典的原子数
        dictsize = size(params.initdict, 2);  
    end
end

% 如果输入参数直接给出待训练字典的原子数，那么就直接使用
% 此项相比上面具有最高优先权
if (isfield(params, 'dictsize')) 
    dictsize = params.dictsize;
end

% 如果原始信号数据个数小于字典原子总数，则抛出错误
if (size(data, 2) < dictsize)
    error('Number of training signals is smaller than number of atoms to train');
end


%% 初始化待训练字典 %%
if (isfield(params, 'initdict'))
    if (any(size(params.initdict) == 1) && all(iswhole(params.initdict(:))))
        % 使用某几个原始信号数据作为初始字典
        D = data(:, params.initdict(1:dictsize));
    else
        if (size(params.initdict,1) ~= size(data, 1) || size(params.initdict, 2) < dictsize)
            % 如果给定的初始字典的行与原始信号数据的维度不一样
            % 或者给定的出使矩阵的列数小于待训练字典原子个数，则抛出错误
            error('Invalid initial dictionary');
        end
        % 按照设定的原子数截取给定的初始化字典矩阵
        % 因为参数里给定的初始字典的原子数目可能比设定的字典的原子数目要多
        D = params.initdict(:, 1:dictsize);  
    end
else
    % 找出原始信号数据列分块后的平方和大于0的块
    % colnorms_squared(data)得到的是一个1x块数的向量
    data_ids = find(colnorms_squared(data) > 1e-6);
    % 返回1到length(data_ids)之间的length(data_ids)个随机整数
    perm = randperm(length(data_ids));
    % 随机选取原始信号数据的一些列初始化字典
    D = data(:, data_ids(perm(1 : dictsize)));
end


%% 归一化字典 %%
D = normcols(D);  % L2范数归一化
err = zeros(1, iternum);  % 残差向量
gerr = zeros(1, iternum);  % 测试误差（泛化误差）向量，是在与训练数据同分布的独立测试样本上的平均损失

if (codemode == CODE_SPARSITY)
    % 基于稀疏度限制的情形，默认是该选项
    errstr = 'RMSE';
else
    % 基于残差限制的情形
    errstr = 'mean atomnum';
end


%%  核心主程序  %%
for iter = 1 : iternum
    printf('KSVD第%d次迭代', iter);
    G = [];
    if (memusage >= MEM_NORMAL)
        G = D' * D;  % 内存使用状态为normal或high时，G需要预计算
    end
    
    %%%%%  稀疏表示  %%%%%
    % 输入参数
    % data：原始信号矩阵； D：字典；
    % XtX：残差向量； G：D'*D； thresh：阈值
    Gamma = sparsecode(data, D, XtX, G, thresh);  % Gamma为稀疏系数矩阵
    
    %%%%%  字典更新  %%%%%
    replaced_atoms = zeros(1, dictsize);  % 每次更新的原子的索引向量初始化为0向量
    unused_sigs = 1 : size(data, 2);  % 未表示的原始信号索引向量
    
    % 确保每一个信号只被处理一次
    p = randperm(dictsize);  % 在1到dictsize之间随机产生dictsize个随机整数
    tid = timerinit('updating atoms', dictsize);  % 初始化计时器
    for j = 1 : dictsize
        % 更新原子
        % 可改动此处，使字典原子正交化
        [D(:, p(j)), gamma_j, data_indices, unused_sigs, replaced_atoms] = optimize_atom(data, D, p(j), Gamma, unused_sigs, replaced_atoms);
        % 更新稀疏系数
        Gamma(p(j), data_indices) = gamma_j;
        
        if (msgdelta > 0)
            % 估计剩余时间，直接显示
            % 输入参数：tid：计时器，j：目前迭代次数，msgdelta：输出显示的的延迟停顿时间
            timereta(tid, j, msgdelta); 
        end
    end
    if (msgdelta > 0)
        printf('updating atoms: iteration %d/%d', dictsize, dictsize);
    end
    
    %%%%%  计算残差  %%%%%
    if (comperr)
        err(iter) = compute_err(D, Gamma, data);
    end
    if (testgen)
        if (memusage >= MEM_NORMAL)
            G = D' * D;
        end
        GammaG = sparsecode(testdata, D, XtXg, G, thresh);
        % 泛化误差，用验证数据计算
        gerr(iter) = compute_err(D, GammaG, testdata);
    end
    
    %%%%%  清理字典中的原子  %%%%%    
    [D, cleared_atoms] = cleardict(D, Gamma, data, muthresh, unused_sigs, replaced_atoms);
    
    %%%%%  打印信息  %%%%%
    info = sprintf('Iteration %d / %d complete', iter, iternum);
    if (printerr)
        info = sprintf('%s, %s = %.4g', info, errstr, err(iter));
    end
    if (printgerr)
        info = sprintf('%s, test %s = %.4g', info, errstr, gerr(iter));
    end
    if (printreplaced)
        info = sprintf('%s, replaced %d atoms', info, sum(replaced_atoms) + cleared_atoms);
    end
    
    if (printiter)
        disp(info);
        if (msgdelta > 0), disp(' '); end
    end
    
end  % 主程序至此结束

end  % 整个程序至此结束，与第一行相呼应



%%
%%%%% 以下均为核心函数 %%%%%

%% 更新原子 %%
% 输入参数：
%     X：原始数据信号矩阵
%     D：字典
%     j：随机索引值，表示字典中的某列，即表示某个原子
%     Gamma：稀疏系数矩阵
%     unused_sigs：未被表示的原始信号
%     replaced_atoms：被取代的原子
% 输出参数：
%     atom：更新后的原子
%     gamma_j：稀疏系数矩阵中第j行中的非0元素构成的向量
%     data_indices：稀疏系数矩阵中第j行中的非0元素的索引位置构成的向量
%     unused_sigs：未被表示的原始信号
%     replaced_atoms：被取代的原子
%
function [atom, gamma_j, data_indices, unused_sigs, replaced_atoms] = optimize_atom(X, D, j, Gamma, unused_sigs, replaced_atoms)

global exactsvd  % 精确SVD

% data samples which use the atom, and the corresponding nonzero
% coefficients in Gamma
% 返回稀疏系数矩阵中第j行中的非0元素及其索引位置
[gamma_j, data_indices] = sprow(Gamma, j);

% 若稀疏系数矩阵中第j行中的元素全为0
if (length(data_indices) < 1)
    maxsignals = 5000;
    perm = randperm(length(unused_sigs));  % length(unused_sigs)就是还未被表示信号个数，初始时为X的列数，每循环一次减1
    perm = perm(1 : min(maxsignals, end));  % 若perm的维数比5000大，则只取其前5000个
    
    % 计算残差平方
    E = sum((X(:, unused_sigs(perm)) - D * Gamma(:, unused_sigs(perm))) .^ 2);
    [d, i] = max(E);  % d为E的每一列的最大值构成的向量，i为每一列最大值的行号构成的向量
    
    % 从信号数据中取相应的列作为原子，并归一化
    atom = X(:, unused_sigs(perm(i)));
    
    if (j > 1)
        % 将已得到的所有字典原子正交化处理
        [DictionaryUpdate, RR] = qr([D(:, 1:(j - 1)), atom], 0);
        atom = DictionaryUpdate(:, j);  % 更新后的原子取正交矩阵的最后一列
    else
        atom = atom ./ norm(atom);
    end
    
    gamma_j = zeros(size(gamma_j));  % 重新变为0向量
    unused_sigs = unused_sigs([1 : perm(i) - 1, perm(i) + 1 : end]);  % 去掉已经被表示的信号
    replaced_atoms(j) = 1;  % 在被替换原子向量的相应位置，即第j处置1
    return;
end

smallGamma = Gamma(:, data_indices);
Dj = D(:, j);

if (exactsvd)
    % 若需精确SVD计算
    % 此处svds将矩阵分解为mx1,1x1,1xn
    % 改动此处可实现正交字典
    EJ = X(:, data_indices) - D * smallGamma + Dj * gamma_j;  % 去掉当前原子成分所造成的误差矩阵
    [atom, s, gamma_j] = svds(EJ, 1);
    if (j > 1)
        % 将已得到的所有字典原子正交化处理
        [DictionaryUpdate, RR] = qr([D(:, 1:(j - 1)), atom], 0);
        atom = DictionaryUpdate(:, j);  % 更新后的原子取正交矩阵的最后一列
        gamma_j = atom \ EJ;
    else
        gamma_j = s * gamma_j;  % 更新字典原子之后也要更新相应的稀疏系数
    end 
else
    % 若无需精确SVD计算
    % 则将原始信号中的某几个列的线性组合再减去[已用的其他的字典原子所得的向量]，又加上[更新之前的该列的字典原子]作为更新的原子
    atom = collincomb(X, data_indices, gamma_j') - D * (smallGamma * gamma_j') + Dj * (gamma_j * gamma_j');
    
    if (j > 1)
        % 将已得到的所有字典原子正交化处理
        [DictionaryUpdate, RR] = qr([D(:, 1:(j - 1)), atom], 0);  % 正交化过程包含归一化
        atom = DictionaryUpdate(:, j);  % 更新后的原子取正交矩阵的最后一列
    else
        atom = atom ./ norm(atom);  % 归一化
    end

    % 更新稀疏系数
    gamma_j = rowlincomb(atom, X, 1 : size(X, 1), data_indices) - (atom' * D) * smallGamma + (atom' * Dj) * gamma_j;  
end

end


%% 稀疏表示 %% 
% 输入参数
%     data：原始信号矩阵
%     D：字典
%     XtX：残差向量 
%     G：D'*D 
%     thresh：阈值
function Gamma = sparsecode(data, D, XtX, G, thresh)

global CODE_SPARSITY codemode  % 迭代限制性模式：基于稀疏度的限制性模式
global MEM_HIGH memusage  % high级别的内存使用模式
global ompfunc ompparams  % omp函数和相关参数

if (memusage < MEM_HIGH)
    % memusage不是high
    % ompparams{:}：ompparams元胞矩阵转化成的列向量
    % 此处调用的是omp2
    Gamma = ompfunc(D, data, G, thresh, ompparams{:});
    
else  
    % memusage是high
    if (codemode == CODE_SPARSITY)
        % 基于稀疏度的限制性模式，调用的是omp
        Gamma = ompfunc(D' * data, G, thresh, ompparams{:});      
    else
        % 基于残差的限制性模式，此时要给定残差向量XtX，调用的是omp2
        Gamma = ompfunc(D' * data, XtX, G, thresh, ompparams{:});
    end
end

end


%% 计算残差 %%
function err = compute_err(D, Gamma, data)

global CODE_SPARSITY codemode

if (codemode == CODE_SPARSITY)
    % 基于稀疏度限制的残差计算
    err = sqrt(sum(reperror2(data, D, Gamma)) / numel(data));
else
    % 基于误差限制的残差计算
    err = nnz(Gamma) / size(data, 2);
end

end


%% 清理字典中的原子 %%
% 输入参数：
%     D：字典
%     Gamma：稀疏系数矩阵
%     X：原始数据信号矩阵
%     muthresh：互相关阈值
%     unused_sigs：未被表示的原始信号
%     replaced_atoms：被取代的原子
function [D, cleared_atoms] = cleardict(D, Gamma, X, muthresh, unused_sigs, replaced_atoms)

use_thresh = 4;  % at least this number of samples must use the atom to be kept

dictsize = size(D, 2);

% 分块计算残差
err = zeros(1, size(X, 2));
blocks = [1 : 3000 : size(X, 2) size(X, 2) + 1];
for i = 1 : length(blocks) - 1
    err(blocks(i) : blocks(i + 1) - 1) = sum((X(:, blocks(i) : blocks(i + 1) - 1) - D * Gamma(:, blocks(i) : blocks(i + 1) - 1)) .^ 2);
end

cleared_atoms = 0;
usecount = sum(abs(Gamma) > 1e-7, 2);

for j = 1:dictsize
    % compute G(:,j)
    Gj = D' * D(:,j);
    Gj(j) = 0;
    
    % replace atom
    if ( (max(Gj .^ 2) > muthresh^2 || usecount(j) < use_thresh) && ~replaced_atoms(j) )
        [y, i] = max(err(unused_sigs));
        D(:, j) = X(:, unused_sigs(i)) / norm(X(:, unused_sigs(i)));
        unused_sigs = unused_sigs([1 : i - 1, i + 1 : end]);
        cleared_atoms = cleared_atoms + 1;
    end
end

end


%% 分块计算残差的平方和
function err2 = reperror2(X,D,Gamma)

err2 = zeros(1, size(X, 2));
blocksize = 2000;
for i = 1 : blocksize : size(X, 2)
    blockids = i : min(i + blocksize - 1, size(X, 2));
    err2(blockids) = sum((X(:, blockids) - D * Gamma(:, blockids)) .^ 2);
end

end


%% 矩阵按列归一化 %%
% 分块计算节约内存
function Y = colnorms_squared(X)
Y = zeros(1, size(X, 2));
blocksize = 2000;  % 块的大小为2000，若数据列数大于2000，就要分块计算
for i = 1 : blocksize : size(X, 2)
    % 块的编号向量
    blockids = i : min(i + blocksize - 1, size(X, 2));
    Y(blockids) = sum(X(:, blockids) .^ 2);  % 对每一块的当前列进行归一化计算
end

end
