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
    legendfontsize --> 11
    labelfontsize --> 14
    tickfontsize --> 12
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
        ribbon := data.dσdt_err
        data.t_range, data.dσdt
    end

    @series begin
        seriestype := :scatter
        yerror := data.Δtot_hera
        label := "Coherent H1"
        color := :blue
        marker := :utriangle
        markersize := 5
        data.tcent_hera, data.dσcoh_hera
    end

    # @series begin
    #     seriestype := :scatter
    #     yerror := [0.2]
    #     label := "Coherent H1"
    #     color := :blue
    #     marker := :utriangle
    #     markersize := 5
    #     [NaN], [NaN]  
    # end
end
