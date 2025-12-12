function [SSE_y, SSE_theta, mse_y, mse_theta, theta_hat] = RLS(theta_old, phi)
    theta_hat = zeros(length(theta_old));
    phi = phi;
    [num, den]=tfdata(sysd,'v');
    
    theta_real = [den(2:end), num(2:end)].';
    noise_variance = noise_variance; % system noise variance
    noise = sqrt(noise_variance) * randn(1, num_samples);
    
    P_old = P_zero;
    theta_hat_zero = theta_hat_zero;
    
    for i=1:num_samples
        % generate output and phi_t
        switch i
        case 1 % at time t=1
            phi_t = [zeros(1,4), zeros(1,4)].';
        case 2
            phi_t = [-y_real(i-1),zeros(1,3), u(i-1),zeros(1,3)].';
        case 3
            phi_t = [-y_real(i-1:-1:i-2).',zeros(1,2), u(i-1:-1:i-2),zeros(1,2)].';
        case 4
            phi_t = [-y_real(i-1:-1:i-3).',zeros(1,1), u(i-1:-1:i-3),zeros(1,1)].';
        otherwise
            phi_t = [-y_real(i-1:-1:i-4).' u(i-1:-1:i-4)].';
        end
    
        phi(i, :) = phi_t.';
        y_real(i) = (phi_t.' * [den(2:end), num(2:end)].') + noise(i);
    
        P_new = P_old - (P_old * (phi_t * phi_t.') * P_old)/(1+phi_t.'*P_old*phi_t);
        K_t = (P_old * phi_t)/(1+phi_t.'*P_old*phi_t);
    
        switch i
            case 1 % at time t=1
                theta_hat(:, i) = theta_hat_zero + K_t * (y_real(i) - phi_t.' * theta_hat_zero);
            otherwise
                theta_hat(:, i) = theta_hat(:, i-1) + K_t * (y_real(i) - phi_t.' * theta_hat(:, i-1));
        end
        y_predicted(i) = phi_t.' * theta_hat(:, i);
        
        P_old = P_new;
    end
    
    error_y = y_predicted-y_real;
    SSE_y = norm(error_y,2)^2;
    mse_y = immse(y_real, y_predicted);

    error_theta = theta_real * ones(1,num_samples) - theta_hat;
    SSE_theta = zeros(8,1);
    mse_theta = zeros(8,1);
    for k = 1:length(theta_real)
        SSE_theta(k) = norm(error_theta(k,:),2)^2;
        mse_theta(k) = mean(error_theta(k,:).^2);
    end
    

    % figures
    f1 = figure();
    f1.Position = [-1000 0 1000 500];
    
    plot(1:num_samples, y_real, "--b", 'DisplayName','real output');
    xlabel("sample time");
    hold on;
    plot(1:num_samples, y_predicted, 'DisplayName','predicted output');
    legend('Location','best');
    title("Output")
    saveas(gcf,'images/q2/' + main_title+ "_summary" + '.jpeg')
    close all

    
    
    f2 = figure();
    f2.Position = [-1500 -200 1000 1300];
    
    for i = 1:length(theta_hat(:,1))
        title_text = "θ_%d";
        subplot(ceil(length(theta_hat(:,1)))/2,2,i);
        plot(1:num_samples, theta_hat(i,:), 'DisplayName','Predicted')
        hold on;
        plot(1:num_samples, ones(size(1:num_samples)) * theta_real(i) , 'DisplayName','Real')
        title(sprintf(title_text, i));
        legend('Location','best');
        xlabel("sample number")
    end
    
    saveas(gcf,'images/q2/' + main_title + '.jpeg')
    close all

    f3 = figure();
    f3.Position = [-1000 0 1000 500];
    sysd_predicted = tf([0, theta_hat(5:8, num_samples).'], [1 , theta_hat(1:4, num_samples).'], sysd.Ts); 
    bode(sysd_predicted);
    hold on; 
    h = bodeplot(sysd);
    setoptions(h,'FreqUnits','Hz','PhaseVisible','off');
    legend("Predicted system", "Real system");
    saveas(gcf,'images/q2/' + main_title+ "_bode" + '.jpeg')
    close all

    
end
