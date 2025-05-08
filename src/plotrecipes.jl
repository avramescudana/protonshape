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