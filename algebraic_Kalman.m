function [Ptt, Ptt_1, xtt, xtt_1] = algebraic_Kalman(data, phase_ns_model)
Ptt = zeros(phase_ns_model.phase_model_dimensions*(phase_ns_model.param_est_window+1), phase_ns_model.phase_model_dimensions);
Ptt_1 = zeros(phase_ns_model.phase_model_dimensions*(phase_ns_model.param_est_window+1), phase_ns_model.phase_model_dimensions);
xtt = zeros((phase_ns_model.param_est_window+1), phase_ns_model.phase_model_dimensions);
xtt_1 = zeros((phase_ns_model.param_est_window+1), phase_ns_model.phase_model_dimensions);

x_prev = phase_ns_model.x_0;
P_prev = phase_ns_model.P0;

for t = 1:(phase_ns_model.param_est_window+1)
    xtt_1(t, :) = (phase_ns_model.A*x_prev)';
    Ptt_1((t-1)*(phase_ns_model.phase_model_dimensions)+1:t*(phase_ns_model.phase_model_dimensions), :) = phase_ns_model.A*P_prev*phase_ns_model.A'+ phase_ns_model.Q2;
    
    K = Ptt_1((t-1)*(phase_ns_model.phase_model_dimensions)+1:t*(phase_ns_model.phase_model_dimensions), :)*(phase_ns_model.C')*pinv(phase_ns_model.C*Ptt_1((t-1)*(phase_ns_model.phase_model_dimensions)+1:t*(phase_ns_model.phase_model_dimensions), :)*phase_ns_model.C' + phase_ns_model.R2);
    epsilon = data(t) - phase_ns_model.C*xtt_1(t, :)';
    
    xtt(t, :) = xtt_1(t, :) + (K*epsilon)';
    Ptt((t-1)*(phase_ns_model.phase_model_dimensions)+1:t*(phase_ns_model.phase_model_dimensions), :) = Ptt_1((t-1)*(phase_ns_model.phase_model_dimensions)+1:t*(phase_ns_model.phase_model_dimensions), :) - K*phase_ns_model.C*Ptt_1((t-1)*(phase_ns_model.phase_model_dimensions)+1:t*(phase_ns_model.phase_model_dimensions), :);
    
    x_prev = xtt(t, :)';
    P_prev = Ptt((t-1)*(phase_ns_model.phase_model_dimensions)+1:t*(phase_ns_model.phase_model_dimensions), :);
end
end