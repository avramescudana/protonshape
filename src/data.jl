function read_coherent_data(file_path::String)
    tcent_hera = Float64[]
    dσcoh_hera = Float64[]
    Δtot_hera = Float64[]

    for line in eachline(file_path)
        if startswith(strip(line), "#") || isempty(strip(line))
            continue
        end

        columns = split(line)
        push!(tcent_hera, parse(Float64, columns[4]))
        push!(dσcoh_hera, parse(Float64, columns[5]))
        push!(Δtot_hera, parse(Float64, columns[6]))
    end

    return tcent_hera, dσcoh_hera, Δtot_hera
end