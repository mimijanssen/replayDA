function params = get_theta_params()
    %% Theta rhythm
    % theta band
    % params.theta.band = [7 12];  
    % epimemage
    params.theta.band = [6 9]; 
    params.theta.sm_twin_s = 0.5; % in s, size of the gaussian smoothing 
    % applied to the theta power before getting the z-score
    params.theta.merge_thr = 0.3; % in s, merge events closer than this
    params.theta.minlen = 0.25; % min len of a theta bout
    % theta phase that will be considered the start of a theta cycle
    params.theta.start_phase = -pi;
    % z-score of theta power above which signal will be considered theta
    % here the through is 0 / 360 as per Wang et al., 2020
    params.theta.sec_zscore_th = -0.2;

    % delta band (used for theta power ratio)
    params.delta.band = [0.5 4]; % Olafsdottir et al 2017

    % params_theta = params.theta;

    % https://www.sciencedirect.com/science/article/pii/S0306452202006693
    % theta 6-9 Hz
    % gamma 40-100 Hz
    % ripples 140-200 Hz
    % some say 250 Hz:
    % https://www.jneurosci.org/content/42/11/2268.abstract 

    %delta/theta

    % delta

end