function xt_prod = matrix_product(A, B, model.state_dimensions)
repeats = length(A)/length;
solution = zeros(length);

for j = 1:repeats
    solution = solution + A(:, ((j-1)*length+1):(j*length));
end

end