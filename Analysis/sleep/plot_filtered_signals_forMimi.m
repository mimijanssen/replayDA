function [all_axes] = plot_filtered_signals(ses_name, restrict_between, pov_tsd, ...
    lfp_tsd_raw, ...
    ripple_data, theta_data, delta_data, rip_plot_opt, ...
    save_ripple_detect_plot, save_fig_path, params)


    % %% options
    % zscore_th = 4; % limit of raw lfp z-score above or below which we 
    % % consider the signal to be noise

    if isempty(restrict_between)
        st = lfp_tsd_raw.tvec(1);
        et = lfp_tsd_raw.tvec(end);
    else
        st = restrict_between.st;
        et = restrict_between.et;
    end
    
    %% restrict the relevant data
    
    % Speed
    if ~isempty(pov_tsd)
        pov_tr = restrict(pov_tsd, st, et);
    else
        pov_tr = [];
    end
    
    % Raw LFP
    lfp_tsd_tr = restrict(lfp_tsd_raw, st, et);
    lfp_unit = lfp_tsd_tr.units;   

    % % New:z-score of raw LFP
    % lfp_raw_zscore_data = zscore(lfp_tsd_tr.data);
    % lfp_raw_zscore = lfp_tsd_tr;
    % lfp_raw_zscore.data = lfp_raw_zscore_data;

    % Ripple band
    lfp_rip_tr = restrict(ripple_data.lfp_tsd, st, et);
    
    
    % Could add ripple power?
    % ripple envelope
    rip_env_tr = restrict(ripple_data.envelope_tsd, st, et);
    % zscore of envelope
    rip_zenv_tr = restrict(ripple_data.zscored_env_tsd, st, et);
    % candidate ripples
    rip_evts_tr = restrict(ripple_data.evts_iv, st, et);
    
    % Theta band
    lfp_theta_tr = restrict(theta_data.lfp_tsd, st, et);
    % theta power
    theta_pow_tr = restrict(theta_data.power, st, et);
    % Theta phase
    theta_ph_tr = restrict(theta_data.phase, st, et);

    keep_ends = 1;
    theta_ivs = restrict(theta_data.theta_iv, st, et, keep_ends);


    % Delta frequency & power
    lfp_delta_tr = restrict(delta_data.lfp_tsd, st, et);
    delta_pow_tr = restrict(delta_data.pow_tsd, st, et);
    
    % Theta over delta power??
    th_ov_de_tr = theta_pow_tr;
    th_ov_de_tr.data = theta_pow_tr.data./delta_pow_tr.data;
    th_ov_de_tr.units = 'none';
    th_ov_de_tr.label = 'theta_over_delta';
    
    %% Create the figure
    fig_title = [ses_name ' ripple detection'];
    this_pos = get(0, 'Screensize');
    this_fig = figure('Name', fig_title, 'Position', this_pos);
    
    % num_rows =9;
    
    layout_spacing_type = 'none'; % Could also be tight
    tiledlayout('vertical', 'TileSpacing', layout_spacing_type); %Tiledlayout is great!
    all_axes = {};
    ind_plot = 1;
    
    % 1. Speed
    if ~isempty(pov_tsd)
        all_axes{ind_plot} = nexttile;
        plot(pov_tr, 'Color', rip_plot_opt.rip_env_col)
        hold on
        % Show vertical line at low-speed threshold
        yline(thisSdata.ls_th, 'k:')
        ylim([0, 80]);
        ylabel({'Speed', ['(' pov_tsd.units ')']})
        set(gca, 'YColor', rip_plot_opt.rip_env_col);
    
        if ~strcmp(layout_spacing_type, 'none')
            title('Running speed');
        end
        ind_plot = ind_plot +1;
    end
    
    % 2. Raw LFP
    all_axes{ind_plot} = nexttile();
    plot(lfp_tsd_tr, 'Color', rip_plot_opt.lfp_col);
    set(gca, 'YColor', rip_plot_opt.lfp_col);
    ylabel({'raw LFP',['(' lfp_unit ')']}); 
    if ~strcmp(layout_spacing_type, 'none')
        title('Raw LFP')
    end

    % % 2.5 z-scored LFP
    % % lfp_raw_zscore
    % all_axes{ind_plot} = nexttile();
    % plot(lfp_raw_zscore, 'Color', rip_plot_opt.lfp_col);
    % set(gca, 'YColor', rip_plot_opt.lfp_col);
    % ylabel({'z-scored LFP'}); 
    % 
    % % at a horizontal line for 'significant' zscore
    % yline(zscore_th, 'r')
    % yline(-zscore_th, 'r')
    % if ~strcmp(layout_spacing_type, 'none')
    %     title('z-scored LFP')
    % end    
    
    % 3. Ripple-band LFP
    ind_plot = ind_plot +1; all_axes{ind_plot} = nexttile();
    % trying to plot envelope first
    
    
    % Show ripple-band -filtered lfp and detected ripples
    dec_rip_cfg = [];
    dec_rip_cfg.display = 'tsd';
    dec_rip_cfg.bgcol = rip_plot_opt.lfp_col; % tsd color outside iv's
    dec_rip_cfg.fgcol = [0.8500 0.3250 0.0980]; % tsd color within iv's
    
    PlotTSDfromIV(dec_rip_cfg, rip_evts_tr, lfp_rip_tr);
    box off
    all_axes{ind_plot}.XAxis.Visible = 'off';
    ylabel({'ripple-filtered LFP',['(' lfp_unit ')']}); 
    set(gca, 'YColor', rip_plot_opt.lfp_col);
    
    % scaling 
    ind_plot = ind_plot +1; all_axes{ind_plot} = nexttile();

    % hold on
    % Add zscore on 2nd axis
    % yyaxis('right');
    plot(rip_zenv_tr, 'Color', rip_plot_opt.rip_env_col , 'LineStyle', '-', ...
        'LineWidth', 0.3);
    % Also change color of the y axis accordingly
    ylabel({'z-scored envelope','(sd)'}); 
    set(gca, 'YColor', rip_plot_opt.rip_env_col);
    
    % Show z-score limit that is used to detect ripples
    yline(params.ripple.sec_zscore_th, '--', 'Color', rip_plot_opt.rip_zscore_line_col);
    % center on zero
    ylims = ylim;
    ylim([-ylims(2), ylims(2)])                
    
    if ~strcmp(layout_spacing_type, 'none')
        title({'Ripple band', ['(' num2str(params.ripple.band(1)), '-' num2str(params.ripple.band(2)),'Hz)']})
    end
    ind_plot = ind_plot +1;
    
    % could add limit of ripple power or ripple / theta power for
    % detection as well
    
    % 4. Add theta band LFP
    all_axes{ind_plot} = nexttile();
    cline(lfp_theta_tr.tvec, lfp_theta_tr.data, theta_ph_tr.data);

    % Add detected theta periods
     % new: add theta-detected periods (z-score of power higher
    % than threshold)
    % keyboard
    % theta_ivs
    this_num_theta = size(theta_ivs.tstart,1);
    hold on
    y_val = -0.1;
    % keyboard
    % loop through each theta period to plot it
    for theta_i = 1:this_num_theta
        plot([theta_ivs.tstart(theta_i), theta_ivs.tend(theta_i)],...
            [1,1]*y_val,'r-', 'LineWidth',4);
    end

    ylabel({'theta-filtered LFP', 'and phase', ['(' lfp_unit ')']}); 
    set(gca, 'YColor', rip_plot_opt.lfp_col);
    set(gca,'fontsize', rip_plot_opt.label_fs)
    
    % 8 Add z-scored theta power
    ind_plot = ind_plot +1; all_axes{ind_plot} = nexttile();

    % % add this on a new axis on the right
    % yyaxis('right');
    
    theta_pow_z = theta_pow_tr;
    theta_pow_z.data = zscore(theta_pow_tr.data);
    plot(theta_pow_z, 'Color', rip_plot_opt.lfp_col);
    
    % yline(params.ripple.sec_zscore_th, '--', 'Color', rip_zscore_line_col);
    yline(params.theta.sec_zscore_th, '--', 'Color', rip_plot_opt.rip_zscore_line_col);
    
    
    ylabel({'theta power' , 'z-scored'}); 
    set(gca, 'YColor', rip_plot_opt.lfp_col);
    ind_plot = ind_plot+1;
    
    % %% 7 Add delta just to see
    % all_axes{ind_plot} = nexttile();
    % plot(lfp_delta_tr, 'Color', lfp_col);
    % ylabel({'delta-filtered LFP', ['(' lfp_unit ')'] }); 
    % set(gca, 'YColor', lfp_col);
    % 
    % ind_plot = ind_plot+1;
    
    % %% 9 Add Delta power just to see
    % all_axes{ind_plot} = nexttile();
    % plot(delta_pow_tr, 'Color', lfp_col);
    % ylabel({'delta power'}); 
    % 
    % set(gca, 'YColor', lfp_col);
    % ind_plot = ind_plot+1;
    % keyboard
    %             %% 10 Add theta/delta power
    %             all_axes{ind_plot} = nexttile();
    %             plot(th_ov_de_tr, 'Color', lfp_col);
    %             ylabel({'theta/delta power'}); 
    % 
    %             set(gca, 'YColor', lfp_col);
    % %             ind_plot = ind_plot+1;
    
    
    %% Formatting
    for ax_i = 1:length(all_axes)
        set(all_axes{ax_i},'fontsize', rip_plot_opt.label_fs)
        box(all_axes{ax_i}, 'off');
        if ax_i == length(all_axes) % Last axis
            xlabel('Time (s)')
        else
            if ~rip_plot_opt.show_x_axis
                all_axes{ax_i}.XAxis.Visible = 'off'; % remove x-axis
            end
        end
    end
    linkaxes([all_axes{:}], 'x')
     % keyboard
    %             keyboard
    
    %% Saving
    
    if save_ripple_detect_plot
        % save_fig_path = paths.lfp_plots_path;
        if ~isfolder(save_fig_path)
            mkdir(save_fig_path);
            disp(['Created folder: ' save_fig_path])
        end
    
        this_fig.InvertHardcopy = 'off';
    %                 this_fig.PaperPositionMode = 'auto';
        this_fig.Color = 'w';
    
        disp(['Saving figure to ' save_fig_path fig_title ])
    
        print (this_fig,'-dpng','-r300',[save_fig_path fig_title])              
	    % Close the fig to save memory
        close(this_fig)
    end
end