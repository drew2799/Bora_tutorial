using Printf
using LaTeXStrings

function plot_measurements(s, flatdata, cov; ells=(0, 2, 4), show=true, save_fn=nothing)

    Ns = length(s)
    
    p1 = scatter(
        s, (s.^2).*flat_data[1:Ns];
        yerr = (s.^2).*sqrt.(diag(cov)[1:Ns]),
        seriestype = :scatter,
        marker = (:circle, 4),
        color = :dodgerblue,
        markerstrokecolor = :dodgerblue,
        capsize = 5,
        label = "",#"ℓ = 0",
        xlabel = L"r [Mpc/h]",
        ylabel = L"\xi_0(r)",
    )

    if 2 in ells
    
        p2 = scatter(
            s, (s.^2).*flat_data[Ns+1:2*Ns];
            yerr = (s.^2).*sqrt.(diag(cov)[Ns+1:2*Ns]),
            seriestype = :scatter,
            marker = (:circle, 4),
            color = :coral,
            markerstrokecolor = :coral,
            capsize = 5,
            label = "",#"ℓ = 2",
            xlabel = L"r [Mpc/h]",
            ylabel = L"\xi_2(r)",
        )
        

        if 4 in ells
        
            p3 = scatter(
                s, (s.^2).*flat_data[2*Ns+1:3*Ns];
                yerr = (s.^2).*sqrt.(diag(cov)[2*Ns+1:3*Ns]),
                seriestype = :scatter,
                marker = (:circle, 4),
                color = :forestgreen,
                markerstrokecolor = :forestgreen,
                capsize = 5,
                label = "",#"ℓ = 4",
                xlabel = L"r [Mpc/h]",
                ylabel = L"\xi_4(r)",
            )

            plt = plot(p1, p2, p3; layout = (1, 3), size = (1200, 400), margin = 5Plots.mm, xformatter=:latex, yformatter=:latex)
        else
            plt = plot(p1, p2; layout = (1, 2), size = (800, 400), margin = 5Plots.mm, xformatter=:latex, yformatter=:latex)
        end
    else
        plt = plot(p1; layout = (1, 1), size = (400, 400), margin = 5Plots.mm, xformatter=:latex, yformatter=:latex)
    end

    if save_fn != nothing 
        savefig(save_fn, plt)
    end
    if show
        display(plt)
    end
end

function plot_bestfit(s, flatdata, cov, bestfit; ells=(0, 2, 4), show=true, save_fn=nothing)

    Ns = length(s)

    if 2 in ells
        if 4 in ells
            layout = (1,3)
            figsize = (1200,400)
        else
            layout = (1,2)
            figsize = (800,400)
        end
    else
        layout = (1,1)
        figsize = (400,400)
    end
    
    plt = plot(layout=layout, size=figsize, margin = 5Plots.mm, xformatter=:latex, yformatter=:latex)

    plot!(s, (s.^2).* bestfit[1,:], color=:navy, subplot=1, label="Best Fit")
    scatter!(plt,
        s, (s.^2).*flat_data[1:Ns];
        yerr = (s.^2).*sqrt.(diag(cov)[1:Ns]),
        seriestype = :scatter,
        marker = (:circle, 4),
        color = :dodgerblue,
        alpha=0.65,
        markerstrokecolor = :dodgerblue,
        capsize = 5,
        label = "",#"ℓ = 0",
        xlabel = L"r [Mpc/h]",
        ylabel = L"\xi_0(r)",
        subplot=1
    )

    if 2 in ells
    
        plot!(s, (s.^2).* bestfit[2,:], color=:firebrick, subplot=2, label="Best Fit")
        scatter!(plt,
            s, (s.^2).*flat_data[Ns+1:2*Ns];
            yerr = (s.^2).*sqrt.(diag(cov)[Ns+1:2*Ns]),
            seriestype = :scatter,
            marker = (:circle, 4),
            color = :coral,
            alpha=0.65,
            markerstrokecolor = :coral,
            capsize = 5,
            label = "",#"ℓ = 2",
            xlabel = L"r [Mpc/h]",
            ylabel = L"\xi_2(r)",
            subplot=2
        )
    end

    if 4 in ells
        plot!(s, (s.^2).* bestfit[3,:], color=:darkgreen, subplot=3, label="Best Fit")
        scatter!(plt,
            s, (s.^2).*flat_data[2*Ns+1:3*Ns];
            yerr = (s.^2).*sqrt.(diag(cov)[2*Ns+1:3*Ns]),
            seriestype = :scatter,
            marker = (:circle, 4),
            color = :forestgreen,
            alpha=0.65,
            markerstrokecolor = :forestgreen,
            capsize = 5,
            label = "",#"ℓ = 4",
            xlabel = L"r [Mpc/h]",
            ylabel = L"\xi_4(r)",
            subplot=3
        )
    end

    if save_fn != nothing 
        savefig(save_fn, plt)
    end
    if show
        display(plt)
    end
end