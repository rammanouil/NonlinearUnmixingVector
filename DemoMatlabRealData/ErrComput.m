function [RMSE, std] = ErrComput(a,a_est)
    [R,N]=size(a);
    err_var = sum((a-a_est).^2)/R;
    e2mean = mean(err_var);
    e2var = var(err_var);
    RMSE = sqrt(e2mean);
    std = sqrt(e2var);
end