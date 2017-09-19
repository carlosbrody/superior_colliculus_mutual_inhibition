# DON'T MODIFY THIS FILE -- the source is in file ProAnti.ipynb


"""
    plot_PA(t, U, V; fignum=1, clearfig=true, rule_and_delay_period=1, target_period=1, post_target_period=1,
        other_unused_params...)

Helper function for plotting ProAnti results
"""
function plot_PA(t, U, V; fignum=1, clearfig=true, rule_and_delay_period=1, target_period=1, post_target_period=1,
    other_unused_params...)
    figure(fignum)
    if clearfig; clf(); end
    
    ax1 = subplot(3,1,1)
    h = plot(t, V'); 
    setp(h[1], color=[0, 0, 1])
    setp(h[2], color=[1, 0, 0])
    setp(h[3], color=[1, 0.5, 0.5])
    setp(h[4], color=[0, 1, 1])
    ylabel("V")

    ax = gca()
    oldlims = [ylim()[1]+0.1, ylim()[2]-0.1]
    ylim(minimum([V[:];oldlims[1]])-0.1, maximum([V[:];oldlims[2]])+0.1)
    yl = [ylim()[1], ylim()[2]]
    vlines([rule_and_delay_period, 
            rule_and_delay_period+target_period,
            rule_and_delay_period+target_period+post_target_period], 
            -0.05, 1.05, linewidth=2)
    if yl[1]<0.02
        yl[1] = -0.02
    end
    if yl[2]>0.98
        yl[2] = 1.02
    end
    ylim(yl)
    grid(true)
    remove_xtick_labels(ax1)
        
    ax2 = subplot(3,1,2)
    hu = plot(t, U')
    oldlims = [ylim()[1]+0.1, ylim()[2]-0.1]
    ylim(minimum([U[:];oldlims[1]])-0.1, maximum([U[:];oldlims[2]])+0.1)
    setp(hu[1], color=[0, 0, 1])
    setp(hu[2], color=[1, 0, 0])
    setp(hu[3], color=[1, 0.5, 0.5])
    setp(hu[4], color=[0, 1, 1])
    ylabel("U"); 
    vlines([rule_and_delay_period, 
            rule_and_delay_period+target_period,
            rule_and_delay_period+target_period+post_target_period], 
            ylim()[1], ylim()[2], linewidth=2)
    remove_xtick_labels(ax2)

    grid(true)
    
    subplot(3,1,3)
    delta = V[1,:] - V[4,:]
    hr = plot(t, delta)
    oldlims = [ylim()[1]+0.1, ylim()[2]-0.1]
    ylim(minimum([delta[:];oldlims[1]])-0.1, maximum([delta[:];oldlims[2]])+0.1)
    vlines([rule_and_delay_period, 
            rule_and_delay_period+target_period,
            rule_and_delay_period+target_period+post_target_period], 
            ylim()[1], ylim()[2], linewidth=2)
    xlabel("t"); ylabel("Pro R - Pro L")
    grid(true)
        
end




