function [PHI, a, sigmaQ, sigmaR]  = modular_parameter_estimation(data, phase_model)
%% Kalman Filtering
[Ptt, Ptt_1, xtt, xtt_1] = algebraic_Kalman(data, phase_model);
%% Backward Smoothing
[PttN, xttN, Ptt_1N, xforwN, PforwN] = backward_smoothing(phase_model, Ptt, Ptt_1, xtt, xtt_1);
%% EM parameter estimation
[PHI, a, sigmaQ, sigmaR]  = em_method(phase_model, data, PttN, xttN, Ptt_1N, xforwN, PforwN);
end