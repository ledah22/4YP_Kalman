function [f, a, sigma]  = debug_parameter_estimation(Ptt, Ptt_1, xtt, model)
if length(xtt) == (model.param_est_length+2)*(2*model.num_rhythms)
    %%Backward Smoothing
    [PttN, xttN, Ptt_1N, xforwN, PforwN] = backward_smoothing(model, Ptt, Ptt_1, xtt);
    
    %%EM parameter estimation
    [f, a, sigma]  = em_method(model, PttN, xttN, Ptt_1N, xforwN, PforwN);
end

end

function [PttN, xttN, Ptt_1N, xforwN, PforwN] = backward_smoothing(model, Ptt, Ptt_1, xtt)
dimensions = model.state_dimensions;
length = model.param_est_length;

PttN = zeros(dimensions * length, dimensions); % For the smoothing step
xttN = zeros(length, dimensions);
Ptt_1N = zeros(dimensions * length, dimensions); % Covariance

coder.varsize('xforwN');
coder.varsize('PforwN');
xforwN = xtt(length+1, :);
PforwN = Ptt((dimensions * length+1):(dimensions * (length+1)), :);


for j= (length+1):(-1):2 % The paper is 0 indexed, MATLAB is 1 indexed
    % Extract matrices for simpler representation
    % Index (j-1) becomes:
    % ((j-2)*dimensions+1):((j-1)*dimensions)
    % Index j becomes:
    % ((j-1)*dimensions+1):(j*dimensions)
    ind_prev = ((j-2)*dimensions+1):((j-1)*dimensions);
    ind_curr = ((j-1)*dimensions+1):(j*dimensions);
    Pt_1_prev = Ptt_1(ind_prev, :);
    Pt_1_curr = Ptt_1(ind_curr, :);
    Pt_prev = Ptt(ind_prev, :);
    xt_prev = xtt(j-1, :);
    
    J = Pt_1_prev\Pt_prev*model.A'; %UNSTABLE INVERSES?
    xttN(j-1, 1:dimensions) = xt_prev+(xforwN-xt_prev*(model.A'))*(J');
    PttN(ind_prev, 1:dimensions) = Pt_prev+J*(PforwN-Pt_1_curr)*J';
    %covariance is symmetric, so Ptt_1N = Pt_1tN = J*PttN hence:
    Ptt_1N(ind_prev, 1:dimensions) = J*PforwN;
    
    xforwN = xttN(j-1, :);
    PforwN = PttN(ind_prev, :);
end

end

function [f, a, sigma, A_smooth, B_smooth, C_smooth]  = em_method(model, PttN, xttN, Ptt_1N, xforwN, PforwN)
dimensions = model.state_dimensions;
l = model.param_est_length;
num_rhythms = model.num_rhythms;

P_A = matrix_sum(PttN, dimensions);
P_B = matrix_sum(Ptt_1N, dimensions);
P_C = matrix_sum([PttN((l+1):((l+1)*dimensions), :); PforwN], dimensions);

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