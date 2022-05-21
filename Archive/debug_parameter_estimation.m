function [f, a, sigmaR, dummy]  = debug_parameter_estimation(Ptt, Ptt_1, xtt, xtt_1, model, f_curr, a_curr, sigmaR_curr)

% Test with
% [a1, a2, a3, adummy] = debug_parameter_estimation(out.Ptt(1:(model.param_est_length+1)*2, :), out.Ptt_1(1:2*(model.param_est_length+1), :), out.xtt(1:(model.param_est_length+1), :), model, 10.5, [0.99 0; 0 0.99], 1)

if length(xtt) == model.param_est_length+1
    %%Backward Smoothing
    [PttN, xttN, Ptt_1N, xforwN, PforwN] = backward_smoothing(model, Ptt, Ptt_1, xtt, xtt_1);
    
    %%EM parameter estimation
    [f, a, sigmaR, dummy]  = em_method(model, PttN, xttN, Ptt_1N, xforwN, PforwN);

else
    f = f_curr;
    a = a_curr;
    sigmaR = sigmaR_curr;
    dummy = [1 2; 3 4];
end

    function [PttN, xttN, Ptt_1N, xNfinal, PNfinal] = backward_smoothing(model, Ptt, Ptt_1, xtt, xtt_1)
        length_param = model.param_est_length; %Limit to a max of 5000, BUT product limit is 10000
        
        PttN = zeros(min(2 * length_param, 10000), 2); % For the smoothing step
        xttN = zeros(min(length_param, 5000), 2);
        Ptt_1N = zeros(min(2 * length_param, 10000), 2); % Covariance
        
        xforwN = xtt(length_param+1, :);
        PforwN = Ptt((2 * length_param+1):(2 * (length_param+1)), :);
        
        xNfinal = xforwN;
        PNfinal = PforwN;
        
        for j= (length_param+1):(-1):2 % The paper is 0 indexed, MATLAB is 1 indexed
            % Extract matrices for simpler representation
            % Index (j-1) becomes:
            % ((j-2)*dimensions_param+1):((j-1)*dimensions_param)
            % Index j becomes:
            % ((j-1)*dimensions_param+1):(j*dimensions_param)
            
            ind_prev = max((j-2)*2+1, 1):min((j-1)*2, 50000);
            ind_curr = max((j-1)*2+1, 1):min(j*2, 50000);
            
            Pt_1_prev = Ptt_1(ind_prev, :);
            Pt_1_curr = Ptt_1(ind_curr, :);
            %Pt_prev = Ptt(ind_prev, :);
            xt_prev = xtt_1(j-1, 1:2);
            
            J = Pt_1_prev*model.A'*pinv(Pt_1_curr);
            xttN(j-1, min(1:2, 10:11)) = xt_prev+(xforwN-xt_prev*(model.A'))*(J');
            PttN(ind_prev, min(1:2, 10:11)) = Pt_1_prev+J*(PforwN-Pt_1_curr)*J';
            %covariance is symmetric, so Ptt_1N = Pt_1tN = J*PttN hence:
            Ptt_1N(ind_prev, min(1:2, 10:11)) = J*PforwN;
            
            xforwN = xttN(j-1, :);
            PforwN = PttN(ind_prev, :);
        end
        
    end

    function [f, a, sigmaR, dummy]  = em_method(model, PttN, xttN, Ptt_1N, xNfinal, PNfinal)
        length_param = model.param_est_length;
        num_rhythms = model.num_rhythms;
        
        P_A = matrix_sum(PttN, 2);
        P_B = matrix_sum(Ptt_1N, 2);
        P_C = matrix_sum([PttN((2+1):end, :); PNfinal], 2);
        
        xtN = [xttN(2:min(length_param, 5000), :); xNfinal];
        x_A_prod = matrix_product(xttN, xttN, 2);
        x_A = matrix_sum(x_A_prod, 2);
        x_B_prod = matrix_product(xtN, xttN, 2);
        x_B = matrix_sum(x_B_prod, 2);
        x_C_prod = matrix_product(xtN, xtN, 2);
        x_C = matrix_sum(x_C_prod, 2);
        
        A_smooth = P_A+x_A;
        B_smooth = P_B+x_B;
        C_smooth = P_C+x_C;
        
        denum_temp = tr(B_smooth, num_rhythms);
        denum = denum_temp(1);
        f_temp = 0.5*pi;
        if denum ~= 0
            f_temp = atan(double(rt(B_smooth, num_rhythms))./double(tr(B_smooth, num_rhythms)));
        end
        
        a_temp = sqrt(rt(B_smooth, num_rhythms).^2+tr(B_smooth, num_rhythms).^2)./tr(A_smooth, num_rhythms);
        sigmaR_temp = double(tr(C_smooth, num_rhythms)-(a_temp.^2).*tr(A_smooth, num_rhythms))./double(2*length_param);
        
        %f_temp2 = mod(double(f_temp(1)/pi+50000), 2)*pi;
        f = f_temp(1);
        %f = length(f_temp);
        a = a_temp(1)*eye(2);
        sigmaR = sigmaR_temp(1);
        dummy = P_A;
    end

end