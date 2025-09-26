
export Sample

struct Sample 
    ids::Vector{String}
end

Base.length(s::Sample) = length(s.ids)


function Sample(s::StationaryPopulation, n::Int)
    if n < 2    
        throw(ArgumentError("Number of samples must be at least 2."))
    end
    if n > s.size
        throw(ArgumentError("Number of samples cannot exceed population size ($(s.size))."))
    end
    ids = [string("pop") for i in 1:n]
    return Sample(ids)
end


function Sample(d::Demography, n::Int; population::Union{Nothing,String}=nothing)
    if n < 2    
        throw(ArgumentError("Number of samples must be at least 2."))
    end

    if isnothing(population) 
        if length(d.populations) > 1
            throw(ArgumentError("Population must be specified when demography has multiple populations. Use `population` keyword argument."))
        else
            population = d.populations[1].id
        end
    end

    if isnothing(get_population_index_by_id(d, population))
        throw(ArgumentError("Population name '$population' not found in demography populations."))
    end

    ids = [population for i in 1:n]
    return Sample(ids)
end

