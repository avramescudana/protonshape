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


@recipe function f(data::CoherentMCData; custom_label="GBW", hera=true)
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

    if hera
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

    @series begin
        seriestype := :path
        label := custom_label
        linewidth := 1
        ribbon := data.dσdt_err
        data.t_range, data.dσdt
    end
end

struct CohIncohMCData
    t_range::Vector{Float64}
    dσdt_coh::Vector{Float64}
    dσdt_coh_err::Vector{Float64}
    dσdt_incoh::Vector{Float64}
    dσdt_incoh_err::Vector{Float64}
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
        # label := "Coherent CQ"
        label := "Coherent"
        linewidth := 1.5
        # ribbon := data.dσdt_coh_err
        ribbon := min.(data.dσdt_coh_err, data.dσdt_coh .- 1e-10)
        # alpha := 0.25 
        # yerror := data.dσdt_coh_err
        data.t_range, data.dσdt_coh
    end

    @series begin
        seriestype := :path
        color := :red
        linewidth := 1.5
        # label := "Incoherent CQ"
        label := "Incoherent"
        # ribbon := data.dσdt_incoh_err
        ribbon := min.(data.dσdt_incoh_err, data.dσdt_incoh .- 1e-10)
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
    α::Float64 = 1.0
    env_func::Function = gaussenv
    a::Float64
end

@recipe function f(cfg::SingleConfiguration)
    fontfamily --> "Computer Modern"
    x_vals = range(-cfg.L, cfg.L, length=cfg.Nx)
    y_vals = range(-cfg.L, cfg.L, length=cfg.Ny)
    X = repeat(x_vals, 1, cfg.Ny)
    Y = repeat(y_vals', cfg.Nx, 1)

    # dens = circmemb_2D(X, Y, cfg.coeff_dict; a=cfg.a)
    if cfg.func_type=="circmemb"
        dens = circmemb_2D(X, Y, cfg.coeff_dict; a=cfg.a)
    elseif cfg.func_type=="Tp"
        dens = Tp_2D(X, Y, cfg.coeff_dict; α=cfg.α, envfunc=cfg.env_func, a=cfg.a)
    else
        error("Unknown function type: $func_type")
    end

    seriestype := :heatmap
    xlabel := L"x\;\mathrm{[GeV^{-1}]}"
    ylabel := L"y\;\mathrm{[GeV^{-1}]}"
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
    α::Float64 = 1.0
    env_func::Function = gaussenv
    a::Float64
    amp::Float16 = 1.0
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
            coeff_dict = Dict((m, n) => mc.amp)
            if mc.func_type=="circmemb"
                dens = circmemb_2D(X, Y, coeff_dict; a=mc.a)
            elseif mc.func_type=="Tp"
                dens = Tp_2D(X, Y, coeff_dict; α=mc.α, envfunc=gaussenv, a=mc.a)
            else
                error("Unknown function type: $func_type")
            end
            @series begin
                seriestype := :heatmap
                xlabel := L"x\;\mathrm{[GeV^{-1}]}"
                ylabel := L"y\;\mathrm{[GeV^{-1}]}"
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

Base.@kwdef struct SuperposedModesPlot
    mmax::Int
    nmax_list::Vector{Int}   # List of nmax values for columns
    L::Float64
    Nx::Int
    Ny::Int
    func_type::String
    α::Float64
    env_func::Function
    a::Float64
    amp::Float64
end

@recipe function f(sp::SuperposedModesPlot)
    fontfamily --> "Computer Modern"
    x_vals = range(-sp.L, sp.L, length=sp.Nx)
    y_vals = range(-sp.L, sp.L, length=sp.Ny)
    X = repeat(x_vals, 1, sp.Ny)
    Y = repeat(y_vals', sp.Nx, 1)

    # Now: rows = m, columns = nmax
    layout := (sp.mmax+1, length(sp.nmax_list))
    size := (900, 900)
    dpi := 900
    rowgap := 0
    colgap := 0
    labelfontsize := 10
    tickfontsize := 8
    titlefontsize := 11

    for m in 0:sp.mmax
        for (col, nmax) in enumerate(sp.nmax_list)
            coeff_dict = Dict((m, n) => sp.amp for n in 1:nmax)
            if sp.func_type == "circmemb"
                dens = circmemb_2D(X, Y, coeff_dict; a=sp.a)
            elseif sp.func_type == "Tp"
                dens = Tp_2D(X, Y, coeff_dict; α=sp.α, envfunc=sp.env_func, a=sp.a)
            else
                error("Unknown function type: $(sp.func_type)")
            end
            @series begin
                seriestype := :heatmap
                xlabel := L"x\;\mathrm{[GeV^{-1}]}"
                ylabel := L"y\;\mathrm{[GeV^{-1}]}"
                color := :inferno
                # title := L"(%$m,%$nmax)"
                # title := L"" * join(["(%$m,%$n)" for n in 1:nmax], " + ")
                labelstr = join(["(" * string(m) * "," * string(n) * ")" for n in 1:nmax], " + ")
                title := LaTeXString(labelstr)
                aspect_ratio := :equal
                colorbar := false
                xlims := (-sp.L, sp.L)
                ylims := (-sp.L, sp.L)
                x_vals, y_vals, dens
            end
        end
    end
end

Base.@kwdef struct SampledSuperposedModesPlot
    m_list::Vector{Int}
    nmax_list::Vector{Int}
    L::Float64
    Nx::Int
    Ny::Int
    func_type::String
    α::Float64
    env_func::Function
    a::Float64
    σ::Float64
    params_shape_base::NamedTuple
end

@recipe function f(sp::SampledSuperposedModesPlot)
    fontfamily --> "Computer Modern"
    x_vals = range(-sp.L, sp.L, length=sp.Nx)
    y_vals = range(-sp.L, sp.L, length=sp.Ny)
    X = repeat(x_vals, 1, sp.Ny)
    Y = repeat(y_vals', sp.Nx, 1)

    layout := (length(sp.m_list), length(sp.nmax_list))
    size := (900, 900)
    dpi := 900
    rowgap := 0
    colgap := 0
    labelfontsize := 10
    tickfontsize := 8
    titlefontsize := 11

    for (row, m) in enumerate(sp.m_list)
        for (col, nmax) in enumerate(sp.nmax_list)
            # Set up params_shape for this panel, including σ
            params_shape = merge(sp.params_shape_base, (
                mn = (m, 1),
                nvals = nmax,
                σ = sp.σ,
            ))
            amps_vec = sample_amp_dict_samem_multin(params_shape)
            full_dict = amps_vec[1]
            coeff_dict_panel = Dict((m, n) => full_dict[(m, n)] for n in 1:nmax)
            if sp.func_type == "circmemb"
                dens = circmemb_2D(X, Y, coeff_dict_panel; a=sp.a)
            elseif sp.func_type == "Tp"
                dens = Tp_2D(X, Y, coeff_dict_panel; α=sp.α, envfunc=sp.env_func, a=sp.a)
            else
                error("Unknown function type: $(sp.func_type)")
            end
            @series begin
                seriestype := :heatmap
                xlabel := L"x\;\mathrm{[GeV^{-1}]}"
                ylabel := L"y\;\mathrm{[GeV^{-1}]}"
                color := :inferno
                labelstr = join(["(" * string(m) * "," * string(n) * ")" for n in 1:nmax], " + ")
                title := LaTeXString(labelstr)
                aspect_ratio := :equal
                colorbar := false
                xlims := (-sp.L, sp.L)
                ylims := (-sp.L, sp.L)
                x_vals, y_vals, dens
            end
        end
    end
end