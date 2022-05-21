function [PttN, xttN, Ptt_1N, xNfinal, PNfinal] = backward_smoothing(phase_ns_model, Ptt, Ptt_1, xtt, xtt_1)
        length_param = phase_ns_model.param_est_window; %Limit to a max of 5000, BUT product limit is 10000
        
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
            
            J = Pt_1_prev*phase_ns_model.A'*pinv(Pt_1_curr);
            xttN(j-1, min(1:2, 10:11)) = xt_prev+(xforwN-xt_prev*(phase_ns_model.A'))*(J');
            PttN(ind_prev, min(1:2, 10:11)) = Pt_1_prev+J*(PforwN-Pt_1_curr)*J';
            %covariance is symmetric, so Ptt_1N = Pt_1tN = J*PttN hence:
            Ptt_1N(ind_prev, min(1:2, 10:11)) = J*PforwN;
            
            xforwN = xttN(j-1, :);
            PforwN = PttN(ind_prev, :);
        end
        
    end