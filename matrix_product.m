function solution = matrix_product(A, B, l)
[repeats, ~] = size(A);
solution = zeros(repeats*l, l);

for j = 1:repeats
    solution(((j-1)*l+1):(j*l), :) = A(j, :)' * B(j, :);
end

end