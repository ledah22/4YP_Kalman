function [PHI, a, sigmaQ, sigmaR]  = em_method(phase_ns_model,data, PttN, xttN, Ptt_1N, xNfinal, PNfinal)
        length_param = phase_ns_model.param_est_window;
        num_rhythms = phase_ns_model.num_rhythms;
        
        P_A = matrix_sum(PttN, 2);
        P_B = matrix_sum(Ptt_1N, 2);
        P_C = matrix_sum([PttN((2+1):end, :); PNfinal], 2);
        
        xtN = [xttN(2:length_param, :); xNfinal];
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
        omega_temp = 0.5*pi;
        if denum ~= 0
            omega_temp = atan(double(rt(B_smooth, num_rhythms))./double(tr(B_smooth, num_rhythms)));
        end
        
        a_temp = sqrt(rt(B_smooth, num_rhythms).^2+tr(B_smooth, num_rhythms).^2)./tr(A_smooth, num_rhythms);
        sigmaQ_temp = double(tr(C_smooth, num_rhythms)-(a_temp.^2).*tr(A_smooth, num_rhythms))./double(2*length_param);
        
        omega = omega_temp(1);
        a = a_temp(1)*eye(2);
        
        sigmaR = sum((data(1:(phase_ns_model.param_est_window))'-phase_ns_model.C*xttN').^2+repmat(phase_ns_model.C, [1, phase_ns_model.param_est_window])*PttN*phase_ns_model.C')/double(length_param);
        PHI = omega*a;
        f = omega/(2*pi*phase_ns_model.Ts);
        sigmaQ = sigmaQ_temp(1);
    end