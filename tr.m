function trace = tr(A, l)
trace = zeros(min(l, 1000), 1);

for j=1:l
    trace(j) = A((j-1)*2+1, 1)+A(j*2, 2);
end