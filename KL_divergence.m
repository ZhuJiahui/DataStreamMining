function distance = KL_divergence(A, B)

if length(A) == length(B)
    for i = 1 : length(A)
        if A(i) < 0.00001
            A(i) = 0.0001;
        end
        if B(i) < 0.00001
            B(i) = 0.0001;
        end
    end
    A = A / sum(A);
    B = B / sum(B);
    distance = sum(A .* log2(A ./ B));
else
    distance = 0;
end

end