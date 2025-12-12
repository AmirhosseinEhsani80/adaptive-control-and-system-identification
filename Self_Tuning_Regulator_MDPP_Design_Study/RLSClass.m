classdef RLSClass < handle
   properties
      P_old
      theta_hat_old
   end
   methods
       function obj = RLSClass(P_init, theta_hat_zero)
        if nargin == 2
            obj.P_old = P_init;
            obj.theta_hat_old = theta_hat_zero;
        end
       end
       function theta_hat_new = update_RLS(obj, y_real, phi_t)

        P_new = obj.P_old - (obj.P_old * (phi_t * phi_t.') * obj.P_old)/(1+phi_t.'*obj.P_old*phi_t);
        K_t = (obj.P_old * phi_t)/(1+phi_t.'*obj.P_old*phi_t);
        theta_hat_new = obj.theta_hat_old + K_t * (y_real - phi_t.' * obj.theta_hat_old);
        obj.theta_hat_old = theta_hat_new;
        obj.P_old = P_new;
       end
       function y_predicted=predict_y(obj, phi_t)
           y_predicted = phi_t.' * obj.theta_hat_old;
       end
   end
end
