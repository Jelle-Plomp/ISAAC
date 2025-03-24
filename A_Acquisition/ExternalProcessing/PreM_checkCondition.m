function [continue_iters,Nit_stable_out] = PreM_checkCondition(apod_cur, apod_prev, Nit_stable_in, PreM)
    % PreM_checkCondition checks if the current apodization satisfies some condition.
    
    % NOTE: no actual method implemented here yet, but this function
    % provides the infrastructure for doing so in the future.

    % Alternatively to checking apodization stability, one could set
    % criteria on the image profile.
    
    % PreM.condition = 
    % PreM.condition_arg = 
    % PreM.Nit_es = Number of stable iterations after which early stop is
    % forced regardless of whether condition is met or not.
    % PreM.maxchange = 0.1; % 0.1 = 10% --> less than 10% change counts as stable

    % apod_cur : apodization for current iteration
    % apod_prev: apodization of previous iteration
    % Nit_stable_out: number of iterations for which the smooth profile has not
    
    continue_iters = true;
    N_it_stable_out = N_it_stable_in + 1;
    warning('A stability check has not yet been implemented but could be implemented here.')

end
