function solution = matrix_sum(A, l)

repeats = length(A)/l;
solution = zeros(l);

for j = 1:repeats
    solution = solution + A(((j-1)*l+1):(j*l), :);
end

end