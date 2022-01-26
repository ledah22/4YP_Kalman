function solution = matrix_product(A, B, l)
repeats = length(A)/l;
solution = zeros(repeats*l, l);

for j = 1:repeats
    solution(((j-1)*l+1):j*l, :) = (A(((j-1)*l+1):(j*l), :)') * B(((j-1)*l+1):(j*l), :);
end

end