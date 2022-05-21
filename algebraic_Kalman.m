function [Ptt, Ptt_1, xtt, xtt_1] = algebraic_Kalman(data, phase_model)
Ptt = zeros(phase_model.phase_model_dimensions*(phase_model.param_est_window+1), phase_model.phase_model_dimensions);
Ptt_1 = zeros(phase_model.phase_model_dimensions*(phase_model.param_est_window+1), phase_model.phase_model_dimensions);
xtt = zeros((phase_model.param_est_window+1), phase_model.phase_model_dimensions);
xtt_1 = zeros((phase_model.param_est_window+1), phase_model.phase_model_dimensions);

x_prev = phase_model.x_0;
P_prev = phase_model.P0;

for t = 1:(phase_model.param_est_window+1)
    xtt_1(t, :) = (phase_model.A*x_prev)';
    Ptt_1((t-1)*(phase_model.phase_model_dimensions)+1:t*(phase_model.phase_model_dimensions), :) = phase_model.A*P_prev*phase_model.A'+ phase_model.Q;
    
    K = Ptt_1((t-1)*(phase_model.phase_model_dimensions)+1:t*(phase_model.phase_model_dimensions), :)*(phase_model.C')*pinv(phase_model.C*Ptt_1((t-1)*(phase_model.phase_model_dimensions)+1:t*(phase_model.phase_model_dimensions), :)*phase_model.C' + phase_model.R);
    epsilon = data(t) - phase_model.C*xtt_1(t, :)';
    
    xtt(t, :) = xtt_1(t, :) + (K*epsilon)';
    Ptt((t-1)*(phase_model.phase_model_dimensions)+1:t*(phase_model.phase_model_dimensions), :) = Ptt_1((t-1)*(phase_model.phase_model_dimensions)+1:t*(phase_model.phase_model_dimensions), :) - K*phase_model.C*Ptt_1((t-1)*(phase_model.phase_model_dimensions)+1:t*(phase_model.phase_model_dimensions), :);
    
    x_prev = xtt(t, :)';
    P_prev = Ptt((t-1)*(phase_model.phase_model_dimensions)+1:t*(phase_model.phase_model_dimensions), :);
end
end