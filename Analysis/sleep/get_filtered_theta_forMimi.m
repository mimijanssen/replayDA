function [theta_data, delta_data] = get_filtered_theta(params, lfp_tsd_raw)
% filter in theta band using params   
    cfg = [];
    cfg.f = params.theta.band ;
    
    theta_lfp = FilterLFP(cfg, lfp_tsd_raw);
    theta_data.lfp_tsd = theta_lfp;
    % keyboard
    cfg_zsc = [];
    cfg_zsc.output = 'power';
    
    theta_data.power = LFPpower(cfg_zsc, theta_data.lfp_tsd);
    
    % smooth it
    
    % sm_twin_s = 0.3; %this doesn't seem to do anything
    sm_twin_s = params.theta.sm_twin_s;
    % Need to compute how many bins this corresponds to
    t_int = median(diff(theta_data.power.tvec));
    sm_twin_bin =sm_twin_s/t_int;
    theta_data.power_sm = theta_data.power;
    theta_data.power_sm.data = ...
    smoothdata(theta_data.power.data, 'gaussian', sm_twin_bin);
    
    % Z-score theta power
    theta_data.z_power = theta_data.power ;
    theta_data.z_power.data = zscore(theta_data.power_sm.data);                
    
    %% Select data chuncks with "high" theta (or, in this case, some theta)
    cfg_dtc = [];
    cfg_dtc.method = 'raw';
    cfg_dtc.threshold = params.theta.sec_zscore_th;
    cfg_dtc.operation =  '>'; % return intervals where threshold is exceeded
    cfg_dtc.merge_thr = params.theta.merge_thr; % merge events closer than this
    cfg_dtc.minlen = params.theta.minlen; % minimum interval length %Wang et al 2020
    fprintf('detecting high z-score theta\n')
    theta_data.theta_iv = TSDtoIV(cfg_dtc, theta_data.z_power);
    
    % Get theta phase
    phase = mod(angle(hilbert(theta_data.lfp_tsd.data)) + ...
        params.theta.start_phase, 2*pi); % get theta phase
    % Note, the 'start phase' value allows to align the start
    % of the phase to the trough of the oscillation
    
    theta_data.phase = theta_data.power;
    theta_data.phase.data = phase;
    
    % TODO add intervals for each theta cycle?
    
    % Also computing delta filtered lfp, but it is not used for
    % now
    
    cfg = [];
    cfg.f = params.delta.band ;
    cfg.display_filter = 0; 
    delta_lfp = FilterLFP(cfg, lfp_tsd_raw);
    delta_data.lfp_tsd = delta_lfp;
    
    cfg_zsc = [];
    cfg_zsc.output = 'power';
    delta_data.pow_tsd = LFPpower(cfg_zsc, delta_data.lfp_tsd);

    
end