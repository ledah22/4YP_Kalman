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