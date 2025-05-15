using Plots
using LaTeXStrings
using RecipesBase

struct CoherentMCData
    t_range::Vector{Float64}
    dσdt::Vector{Float64}
    dσdt_err::Vector{Float64}
    tcent_hera::Vector{Float64}
    dσcoh_hera::Vector{Float64}
    Δtot_hera::Vector{Float64}
end


@recipe function f(data::CoherentMCData; custom_label="GBW")
    fontfamily --> "Computer Modern"
    framestyle --> :box
    legendfontsize --> 10
    labelfontsize --> 14
    tickfontsize --> 11
    size --> (500, 400)
    foreground_color_legend --> nothing
    background_color_legend --> nothing
    dpi --> 600

    xlabel --> L"|t|\;[\mathrm{GeV}^2]"
    ylabel --> L"\mathrm{d}\sigma/\mathrm{d}|t|\;[\mathrm{nb}/\mathrm{GeV}^2]" 
    xlims --> (0.0, 1.0)
    yticks --> :auto
    yaxis --> :log10

    @series begin
        seriestype := :path
        label := custom_label
        linewidth := 1
        ribbon := data.dσdt_err
        data.t_range, data.dσdt
    end

    @series begin
        seriestype := :scatter
        yerror := data.Δtot_hera
        label := "Coherent H1"
        color := :blue
        marker := :utriangle
        markersize := 3.5
        data.tcent_hera, data.dσcoh_hera
    end
end

struct CohIncohMCData
    t_range::Vector{Float64}
    dσdt_coh::Vector{Float64}
    dσdt_coh_err::Vector{Float64}
    dσdt_incoh::Vector{Float64}
    tcent_coh_hera::Vector{Float64}
    dσdt_coh_hera::Vector{Float64}
    Δtot_coh_hera::Vector{Float64}
    tcent_incoh_hera::Vector{Float64}
    dσdt_incoh_hera::Vector{Float64}
    Δtot_incoh_hera::Vector{Float64}
end

@recipe function f(data::CohIncohMCData)
    fontfamily --> "Computer Modern"
    framestyle --> :box
    legendfontsize --> 10
    labelfontsize --> 14
    tickfontsize --> 11
    size --> (500, 400)
    foreground_color_legend --> nothing
    background_color_legend --> nothing
    dpi --> 600

    xlabel --> L"|t|\;[\mathrm{GeV}^2]"
    ylabel --> L"\mathrm{d}\sigma/\mathrm{d}t\;[\mathrm{nb}/\mathrm{GeV}^2]"
    xlims --> (0.0, 2.5)
    ylims --> (10^(-1), 10^3)
    yticks --> :auto
    yaxis --> :log10

    @series begin
        seriestype := :path
        color := :blue
        label := "Coherent CQ"
        linewidth := 1
        # ribbon := data.dσdt_coh_err
        data.t_range, data.dσdt_coh
    end

    @series begin
        seriestype := :path
        color := :red
        linewidth := 1
        label := "Incoherent CQ"
        data.t_range, data.dσdt_incoh
    end

    @series begin
        seriestype := :scatter
        yerror := data.Δtot_coh_hera
        label := "Coherent H1"
        color := :blue
        marker := :utriangle
        markersize := 3.5
        data.tcent_coh_hera, data.dσdt_coh_hera
    end

    @series begin
        seriestype := :scatter
        yerror := data.Δtot_incoh_hera
        label := "Incoherent H1"
        color := :red
        marker := :utriangle
        markersize := 3.5
        data.tcent_incoh_hera, data.dσdt_incoh_hera
    end
end

"""
Plot multipe configurations for the circular membrane
"""

Base.@kwdef struct SingleConfiguration
    coeff_dict::Dict{Tuple{Int,Int}, Float64}
    L::Float64
    Nx::Int
    Ny::Int
    func_type::String
    alpha::Float64 = 0.125
    env_func::Function = gaussenv
    a::Float64
end

@recipe function f(cfg::SingleConfiguration)
    x_vals = range(-cfg.L, cfg.L, length=cfg.Nx)
    y_vals = range(-cfg.L, cfg.L, length=cfg.Ny)
    X = repeat(x_vals, 1, cfg.Ny)
    Y = repeat(y_vals', cfg.Nx, 1)

    # dens = circmemb_2D(X, Y, cfg.coeff_dict; a=cfg.a)
    if cfg.func_type=="circmemb"
        dens = circmemb_2D(X, Y, cfg.coeff_dict; a=cfg.a)
    elseif cfg.func_type=="Tp"
        dens = Tp_2D(X, Y, cfg.coeff_dict; alpha=cfg.alpha, envfunc=cfg.env_func, a=cfg.a)
    else
        error("Unknown function type: $func_type")
    end

    seriestype := :heatmap
    xlabel := L"x"
    ylabel := L"y"
    color := :inferno
    cbar_title := L"T(\theta, r)"
    aspect_ratio := :equal
    xlims := (-cfg.L, cfg.L)
    ylims := (-cfg.L, cfg.L)
    size := (440, 400)

    x_vals, y_vals, dens
end

Base.@kwdef struct MultipleConfigurations
    mmax::Int
    nmax::Int
    L::Float64
    Nx::Int
    Ny::Int
    func_type::String
    alpha::Float64 = 0.125
    env_func::Function = gaussenv
    a::Float64
end

@recipe function f(mc::MultipleConfigurations)
    fontfamily --> "Computer Modern"
    x_vals = range(-mc.L, mc.L, length=mc.Nx)
    y_vals = range(-mc.L, mc.L, length=mc.Ny)
    X = repeat(x_vals, 1, mc.Ny)
    Y = repeat(y_vals', mc.Nx, 1)

    layout := (mc.mmax+1, mc.nmax)
    size := (900, 900)
    dpi := 900
    # colorbar := true
    # right_margin := 10Plots.mm
    rowgap := 0
    colgap := 0
    labelfontsize := 10
    tickfontsize := 8
    titlefontsize := 11

    # Flatten all heatmaps into a tuple of series
    for m in 0:mc.mmax
        for n in 1:mc.nmax
            coeff_dict = Dict((m, n) => 1.0)
            if mc.func_type=="circmemb"
                dens = circmemb_2D(X, Y, coeff_dict; a=mc.a)
            elseif mc.func_type=="Tp"
                dens = Tp_2D(X, Y, coeff_dict; alpha=mc.alpha, envfunc=gaussenv, a=mc.a)
            else
                error("Unknown function type: $func_type")
            end
            @series begin
                seriestype := :heatmap
                xlabel := L"x"
                ylabel := L"y"
                color := :inferno
                title := L"(%$m,%$n)"
                aspect_ratio := :equal
                colorbar := false
                xlims := (-mc.L, mc.L)
                ylims := (-mc.L, mc.L)
                x_vals, y_vals, dens
            end
        end
    end
end