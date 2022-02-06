function [f, a, sigma]  = em_method(model, PttN, xttN, Ptt_1N, xforwN, PforwN)
dimensions = model.state_dimensions;
l = model.param_est_length;
num_rhythms = model.num_rhythms;

P_A = matrix_sum(PttN, dimensions);
P_B = matrix_sum(Ptt_1N, dimensions);
P_C = matrix_sum([PttN((dimensions+1):end, :); PforwN], dimensions);

xtN = [xttN(2:l, :); xforwN];
x_A_prod = matrix_product(xttN, xttN, dimensions);
x_A = matrix_sum(x_A_prod, dimensions);
x_B_prod = matrix_product(xtN, xttN, dimensions);
x_B = matrix_sum(x_B_prod, dimensions);
x_C_prod = matrix_product(xtN, xtN, dimensions);
x_C = matrix_sum(x_C_prod, dimensions);

A_smooth = P_A+x_A;
B_smooth = P_B+x_B;
C_smooth = P_C+x_C;

f = atan(rt(B_smooth, num_rhythms)./tr(B_smooth, num_rhythms));
a = sqrt(rt(B_smooth, num_rhythms).^2+tr(B_smooth, num_rhythms).^2)./tr(A_smooth, num_rhythms);
sigma = (tr(C_smooth, num_rhythms)-a.^2.*tr(A_smooth, num_rhythms))/(2*l);
end