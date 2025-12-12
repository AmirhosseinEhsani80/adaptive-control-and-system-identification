classdef ELSClass < handle
   properties
      P_old
      theta_hat_old
      theta_epsilon_old
      phi_epsilon_old
   end
   methods
       function obj = ELSClass(P_init, theta_hat_zero, theta_epsilon_zero)
        if nargin == 3
            obj.P_old = P_init;
            obj.theta_hat_old = [theta_hat_zero;theta_epsilon_zero];
            obj.phi_epsilon_old = zeros([length(theta_epsilon_zero),1]);
        end
       end
       function theta_hat_new = update_ELS(obj, y_real, phi_t_without_noise)
        phi_t = [phi_t_without_noise; obj.phi_epsilon_old];
        P_new = obj.P_old - (obj.P_old * (phi_t * phi_t.') * obj.P_old)/(1+phi_t.'*obj.P_old*phi_t);
        K_t = (obj.P_old * phi_t)/(1+phi_t.'*obj.P_old*phi_t);
        theta_hat_new = obj.theta_hat_old + K_t * (y_real - phi_t.' * obj.theta_hat_old);
        obj.theta_hat_old = theta_hat_new;
        obj.P_old = P_new;
        new_epsilon = y_real - obj.predict_y(phi_t);
        obj.phi_epsilon_old = [new_epsilon; obj.phi_epsilon_old(1:end-1)];
        
       end
       function y_predicted=predict_y(obj, phi_t)
           y_predicted = phi_t.' * obj.theta_hat_old;
       end
   end
end
