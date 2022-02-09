function acert = rt(A, l)
acert = zeros(min(l, 1000), 1);

for j=1:l
    acert(j) = A(j*2, 1)-A((j-1)*2+1, 2);

end