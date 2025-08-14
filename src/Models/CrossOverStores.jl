module CrossoverStores


export AbstractCrossoverStore, MemoryCrossoverStore, VectorCrossoverStore, 
    newid!, getparentat, gettime, genome_length

abstract type AbstractCrossoverStore end


Base.show(io::IO, c::AbstractCrossoverStore) = print(io, "$(typeof(c))")


# -------------------------------------------------


struct MemoryCrossoverStore <: AbstractCrossoverStore
    data::Vector{Int64}
    genome_length::Int64
end


"""
    MemoryCrossoverStore(genome_length::Int64)

Create a new `MemoryCrossoverStore` with an empty data vector and a genome length of `genome_length`.
"""
MemoryCrossoverStore(genome_length) = MemoryCrossoverStore(Vector{Int64}(), genome_length)



"""
    newid!(c::MemoryCrossoverStore, id1, id2, time, vpos)

Push a new cross-over event to the store. The cross-over event is between
`id1` and `id2` at time `time` and the cross-over positions are stored in `vpos`.

Return the id of the new event.
"""
function newid!(c::MemoryCrossoverStore, id1::Int64, id2::Int64, time::Float64, vpos::Vector{Int64})
    myid = nextid(c)
    push!(c.data, id1)
    push!(c.data, id2)
    push!(c.data, reinterpret(Int64, time))
    push!(c.data, length(vpos))
    @assert issorted(vpos)
    Base.append!(c.data, vpos)
    myid
end

"""
    newid!(c::MemoryCrossoverStore, time)

Push a new unrelated initial event to the store and return the id.
"""
newid!(c::MemoryCrossoverStore, time::Float64) = newid!(c, 0, 0, time, Int64[])
    

"""
    getparentat(c::MemoryCrossoverStore, i::Int64, pos::Int64)

Return the parent id and the next break position after `pos` in the genome.
"""
function getparentat(c::MemoryCrossoverStore, i::Int64, pos::Int64)
    npos = getnpos(c, i)
    j = searchsortedfirst(c.data, pos, i + 4, i + 3 + npos, Base.Order.Forward)
    pid = iseven(j) == iseven(i) ? getid1(c, i) : getid2(c, i)
    nextbreak = j > i + 3 + npos ? c.genome_length : c.data[j]
    pid::Int64, nextbreak::Int64
end



nextid(c::MemoryCrossoverStore) = length(c.data) + 1
getid1(c::MemoryCrossoverStore, i) = @inbounds c.data[i]
getid2(c::MemoryCrossoverStore, i) = @inbounds c.data[i+1]
gettime(c::MemoryCrossoverStore, i) = @inbounds reinterpret(Float64, c.data[i+2])
getnpos(c::MemoryCrossoverStore, i) = @inbounds c.data[i+3]
getvpos(c::MemoryCrossoverStore, i) = @inbounds c.data[i+4:i+3+getnpos(c, i)]


# -------------------------------------------------

struct CrossoverEvent
    id1::Int64
    id2::Int64
    time::Float64
    vpos::Vector{Int64}
end


struct VectorCrossoverStore <: AbstractCrossoverStore
    data::Vector{CrossoverEvent}
    genome_length::Int64
end


"""
    VectorCrossoverStore(genome_length::Int64)

Create a new `VectorCrossoverStore` with an empty data vector and a genome length of `genome_length`.
"""
VectorCrossoverStore(genome_length) = VectorCrossoverStore(Vector{CrossoverEvent}(), genome_length)



"""
    newid!(c::VectorCrossoverStore, id1, id2, time, vpos)

Push a new cross-over event to the store. The cross-over event is between
`id1` and `id2` at time `time` and the cross-over positions are stored in `vpos`.

Return the id of the new event.
"""
function newid!(c::VectorCrossoverStore, id1::Int64, id2::Int64, time::Float64, vpos::Vector{Int64})
    myid = nextid(c)
    push!(c.data, CrossoverEvent(id1, id2, time, vpos))
    myid
end

"""
    newid!(c::VectorCrossoverStore, time)

Push a new unrelated initial event to the store and return the id.
"""
newid!(c::VectorCrossoverStore, time::Float64) = newid!(c, 0, 0, time, Int64[])
    

"""
    getparentat(c::VectorCrossoverStore, i::Int64, pos::Int64)

Return the parent id and the next break position after `pos` in the genome.
"""
function getparentat(c::VectorCrossoverStore, i::Int64, pos::Int64)
    vpos = getvpos(c, i)
    npos = length(vpos)
    j = searchsortedfirst(vpos, pos)
    pid = isodd(j) ? getid1(c, i) : getid2(c, i)
    nextbreak = j > npos ? c.genome_length : vpos[j]
    pid::Int64, nextbreak::Int64
end



nextid(c::VectorCrossoverStore) = length(c.data) + 1
getid1(c::VectorCrossoverStore, i) = @inbounds c.data[i].id1
getid2(c::VectorCrossoverStore, i) = @inbounds c.data[i].id2
gettime(c::VectorCrossoverStore, i) = @inbounds c.data[i].time
getnpos(c::VectorCrossoverStore, i) = length(getvpos(c, i))
getvpos(c::VectorCrossoverStore, i) = @inbounds c.data[i].vpos




end # module

