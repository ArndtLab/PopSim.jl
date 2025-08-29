using TestItems
using TestItemRunner
@run_package_tests verbose=true

  




@testitem "StationaryPopulation" begin
    sd = StationaryPopulation(2, 100)
    @test sd.ploidy == 2
    @test sd.size == 100
    @test string(sd) == "StationaryPopulation(ploidy=2, size=100)"

    @test_throws ArgumentError StationaryPopulation(3, 100)  # Invalid ploidy
    @test_throws ArgumentError StationaryPopulation(0)    # Invalid size


    @test APop.ploidy(sd) == 2
    @test size(sd) == 100


    @test TNvector(sd, 4400) == [4400, 100]
end

@testitem "TN for varying pop" begin
    d = Demography()
    add_population!(d, Population(id = "pop1", description = "Population 1", size = 100, growth_rate = 0.0, time_offset = 0))
    set_end_time!(d, 5000)

    @test TNvector(d, 2222) == [2222, 100]

    add_event!(d, PopulationSizeEvent(4000, "pop1", 10))
    add_event!(d, PopulationSizeEvent(4020, "pop1", 10000))

    println(APop.summary(d))

    @test TNvector(d, 2222) == [2222, 100, 20, 10, 980, 10000]
end

@testitem "Demography" begin
    d = Demography()
    @test d.start_time == 0
    @test d.end_time == 0
    @test length(d.populations) == 0
    @test length(d.events) == 0

    add_population!(d, Population(id = "pop1", description = "Population 1", size = 100, growth_rate = 0.0, time_offset = 0))
    set_end_time!(d, 1000)
    @test length(d.population_sizes) == 1
    @test d.population_sizes[1][0] == 100
    @test d.population_sizes[1][1000] == 100

    add_event!(d, PopulationSizeEvent(500, "pop1", 150))
    @test d.population_sizes[1][500-1] == 100
    @test d.population_sizes[1][500] == 150
    @test d.population_sizes[1][1000] == 150

    add_population!(d, Population(id = "pop2", description = "Population 2", size = 400, growth_rate = 0.0, time_offset = 0))
    set_migration!(d, "pop1", "pop2", 0.1)
    @test d.migration[2, 1] == 0.1
    @test d.migration[2, 2] == 0.9

    add_population!(d, Population(id = "pop3", description = "Population 3", size = 400, growth_rate = 0.0, time_offset = 0))
    @test d.migration[2, 1] == 0.1
    @test d.migration[2, 2] == 0.9
    @test d.migration[3, 3] == 1.0

    @test_throws ArgumentError TNvector(d, 1000)  
    @test_throws ArgumentError APop.get_population_index_by_id(d, "pop4")
end



@testitem "Genome" begin
    ur = UniformRate(0.1)
    g = Genome(ur, ur, 1000)
    @test length(g) == 1000
    @test APop.rate(recombination(g)) == 0.1
    @test APop.rate(mutation(g)) == 0.1
    @test string(g) == "Genome(recombination=Uniform(rate=0.1), mutation=Uniform(rate=0.1), length=1000)"

    g = Genome(recombination_rate=0.1, mutation_rate=0.1, length=1000)
    @test length(g) == 1000
    @test string(g) == "Genome(recombination=Uniform(rate=0.1), mutation=Uniform(rate=0.1), length=1000)"

    g = Genome(recombination_rate=0.1, length=1000)
    @test length(g) == 1000
    @test string(g) == "Genome(recombination=Uniform(rate=0.1), mutation=Uniform(rate=0.0), length=1000)"

end








@testitem "Segment" begin
    s = Segment(1, 10)
    @test first(s) == 1
    @test last(s) == 10
    @test length(s) == 10  # Inclusive length
    @test string(s) == "Segment(first=1, last=10)"

    s = Segment(5, 5)
    @test first(s) == 5
    @test last(s) == 5
    @test length(s) == 1

end



@testitem "Rate Distributions" begin
    # using StatsBase
    rate = 0.1
    u = UniformRate(rate)
    @test average_rate(u) == rate

    S = 100000
    for L in [1,10,100,1000], dt in [0.1, 1.0, 10]
        x = mapreduce(+, 1:S) do i
            length(APop.sample(u, dt, 1, L))
        end
        @test abs(x / (L * S * rate * dt) - 1) < 0.1
    end

    nu = NonUniformRate([0.1, 0.2, 0.3, 0.4])
    @test average_rate(nu) == 0.25

    L = 10000
    rates = [i % 10 == 1 ? 0.1 : 0.01 for i in 1:L]
    nu = NonUniformRate(rates)
    @test average_rate(nu) ≈ 0.019

    s = APop.sample(nu, 1.0, 1, L)
    @test all(x -> x in 1:L, s)
    ss = map(1:10) do m
        sum(i -> i % 10 == m % 10, s)
    end
    @test ss[1] > 4 * ss[2]

    
    s = APop.sample(nu, 100000.0, 2, 10)
    @test all(x -> x in 2:10, s)
    ss = map(1:10) do m
        sum(i -> i % 10 == m % 10, s)
    end
    @test abs(sum(ss) - 100000.0 * 0.01 * 9) < 0.1 * 100000.0 * 0.01 * 9
end

@testitem "UniformRateDistribution & JC correction" begin
    for i in 1:100
        rate = 0.1
        L = 100
        u = UniformRate(rate)

        r = APop.sample(u, 10*L/rate, 1, L; multiple_hits = :ignore)
        @test all(x -> x in 1:L, r)
        @test length(r) > 2*L

        r = APop.sample(u, 100*L/rate, 1, L; multiple_hits = :as_one)
        @test all(x -> x in 1:L, r)
        @test allunique(r)
        @test length(r) == L

        L = 10000
        r = APop.sample(u, 100*L/rate, 1, L; multiple_hits = :JCcorrect)
        @test all(x -> x in 1:L, r)
        @test allunique(r)
        @test 0.65 * L < length(r) < 0.85 * L
    end
end





@testitem "ARG" begin
    tree = CoalescentTreeTwoLineages(1, 100.0)
    @test tree.root_id == 1
    @test tree.timespan == 100.0

    # Test ARG
    segment = Segment(1, 10)
    arg = ARGsegment(segment, tree)
    @test arg.segment == segment

    @test first(arg) == 1
    @test last(arg) == 10
    @test length(arg) == 10
    @test timespan(arg) == 100.0

end

@testitem "IBS kwargs" begin
    tau = 100.0
    mut = UniformRate(0.1)
    tree = CoalescentTreeTwoLineages(1, tau)

    ibds = [ARGsegment(Segment(1, 1000), tree), ARGsegment(Segment(1001, 2000), tree)]
    ibs = collect(APop.IBSIterator(ibds,mut))
    @test length(ibs) > 3000
    @test sum(length, ibs) === 2000

    ibs = collect(APop.IBSIterator(ibds, mut, multiple_hits = :as_one))
    @test length(ibs) in [1999,2000,2001]
    @test sum(length, ibs) === 2000

    ibs = collect(APop.IBSIterator(ibds, mut, multiple_hits = :JCcorrect))
    @test length(ibs) < 1750
    @test sum(length, ibs) === 2000

end



@testitem "MemoryCrossoverStore" begin
   
    using APop.CrossoverStores

    c = MemoryCrossoverStore(100)

    @test CrossoverStores.nextid(c) == 1

    CrossoverStores.newid!(c, Int64(0), Int64(0), 0.0, Int64[])
    @test CrossoverStores.nextid(c) > 1
    @test CrossoverStores.getid1(c, 1) == 0
    @test CrossoverStores.getid2(c, 1) == 0
    @test CrossoverStores.gettime(c, 1) == 0.0
    @test CrossoverStores.getnpos(c, 1) == 0
    @test CrossoverStores.getvpos(c, 1) == Int64[]


    i = CrossoverStores.nextid(c)
    j = CrossoverStores.newid!(c, Int64(12), Int64(17), 4.5e5, [2,5,6])
    @test i == j


    @test CrossoverStores.getid1(c, i) == 12
    @test CrossoverStores.getid2(c, i) == 17
    @test CrossoverStores.gettime(c, i) == 4.5e5
    @test CrossoverStores.getnpos(c, i) == 3
    @test CrossoverStores.getvpos(c, i) == [2,5,6]

    @test CrossoverStores.getparentat(c, i, 1) == (12, 2)
    @test CrossoverStores.getparentat(c, i, 2) == (12, 2)
    @test CrossoverStores.getparentat(c, i, 3) == (17, 5)
    @test CrossoverStores.getparentat(c, i, 4) == (17, 5)
    @test CrossoverStores.getparentat(c, i, 5) == (17, 5)
    @test CrossoverStores.getparentat(c, i, 6) == (12, 6)
    @test CrossoverStores.getparentat(c, i, 7) == (17, c.genome_length)


    i = CrossoverStores.nextid(c)
    CrossoverStores.newid!(c, Int64(12), Int64(17), 4.5e5, [2,5,6])

    @test CrossoverStores.getid1(c, i) == 12
    @test CrossoverStores.getid2(c, i) == 17
    @test CrossoverStores.gettime(c, i) == 4.5e5
    @test CrossoverStores.getnpos(c, i) == 3
    @test CrossoverStores.getvpos(c, i) == [2,5,6]

    @test CrossoverStores.getparentat(c, i, 1) == (12, 2)
    @test CrossoverStores.getparentat(c, i, 2) == (12, 2)
    @test CrossoverStores.getparentat(c, i, 3) == (17, 5)
    @test CrossoverStores.getparentat(c, i, 4) == (17, 5)
    @test CrossoverStores.getparentat(c, i, 5) == (17, 5)
    @test CrossoverStores.getparentat(c, i, 6) == (12, 6)
    @test CrossoverStores.getparentat(c, i, 7) == (17, c.genome_length)


end



@testitem  "VectorCrossoverStore" begin

    using APop.CrossoverStores

    c = VectorCrossoverStore(100)

    @test CrossoverStores.nextid(c) == 1

    CrossoverStores.newid!(c, 0, 0, 0.0, Int64[])
    @test CrossoverStores.nextid(c) > 1
    @test CrossoverStores.getid1(c, 1) == 0
    @test CrossoverStores.getid2(c, 1) == 0
    @test CrossoverStores.gettime(c, 1) == 0.0
    @test CrossoverStores.getnpos(c, 1) == 0
    @test CrossoverStores.getvpos(c, 1) == Int64[]


    i = CrossoverStores.nextid(c)
    j = CrossoverStores.newid!(c, 12, 17, 4.5e5, [2,5,6])
    @test i == j


    @test CrossoverStores.getid1(c, i) == 12
    @test CrossoverStores.getid2(c, i) == 17
    @test CrossoverStores.gettime(c, i) == 4.5e5
    @test CrossoverStores.getnpos(c, i) == 3
    @test CrossoverStores.getvpos(c, i) == [2,5,6]

    @test CrossoverStores.getparentat(c, i, 1) == (12, 2)
    @test CrossoverStores.getparentat(c, i, 2) == (12, 2)
    @test CrossoverStores.getparentat(c, i, 3) == (17, 5)
    @test CrossoverStores.getparentat(c, i, 4) == (17, 5)
    @test CrossoverStores.getparentat(c, i, 5) == (17, 5)
    @test CrossoverStores.getparentat(c, i, 6) == (12, 6)
    @test CrossoverStores.getparentat(c, i, 7) == (17, c.genome_length)


    i = CrossoverStores.nextid(c)
    CrossoverStores.newid!(c, 12, 17, 4.5e5, [2,5,6])

    @test CrossoverStores.getid1(c, i) == 12
    @test CrossoverStores.getid2(c, i) == 17
    @test CrossoverStores.gettime(c, i) == 4.5e5
    @test CrossoverStores.getnpos(c, i) == 3
    @test CrossoverStores.getvpos(c, i) == [2,5,6]

    @test CrossoverStores.getparentat(c, i, 1) == (12, 2)
    @test CrossoverStores.getparentat(c, i, 2) == (12, 2)
    @test CrossoverStores.getparentat(c, i, 3) == (17, 5)
    @test CrossoverStores.getparentat(c, i, 4) == (17, 5)
    @test CrossoverStores.getparentat(c, i, 5) == (17, 5)
    @test CrossoverStores.getparentat(c, i, 6) == (12, 6)
    @test CrossoverStores.getparentat(c, i, 7) == (17, c.genome_length)
end



@testitem "Individual" begin
    i = Individual(1, 2)
    @test i.alleles[1] == 1
    @test i.alleles[2] == 2 
    @test ploidy(i) == 2
    @test i[1] == 1
    @test i[2] == 2
    @test length(i) == 2
    @test randallele(i) in i.alleles
end

@testitem "WrightFisher - ARG etc" begin
    
    population_size = 100
    mutation_rate = 2e-8
    recombination_rate = 1e-8
    L = 1_000_000_000

    
    d = Demography()
    add_population!(d, Population(id = "pop1", description = "Population 1", size = population_size, growth_rate = 0.0, time_offset = 0))
    set_end_time!(d, 4000)
    
    g = Genome(UniformRate(recombination_rate), UniformRate(mutation_rate),  L)

    model = WrightFisher()


    anc = sim_ancestry(model, d, g)
    @test length(anc.alives) == length(d.populations)




    indv = anc.alives[1][1]
    @test indv[1] > 2 * population_size
    @test indv[2] > 2 * population_size

    ARG = collect(get_ARGsegments(anc, [indv[1], indv[2] ]))
    @test length(ARG) > 0
    @test all(seg -> length(seg) > 0, ARG)
    @test sum(length, ARG) == L
    @test all(seg -> iscoalescent(seg), ARG)
    @test all(seg -> 0.0 <= timespan(seg) < Inf, ARG)

    ibs = collect(IBSIterator(ARG, anc.genome.mutation))
    @test length(ibs) > 0
    @test all(seg -> length(seg) >= 0, ibs)
    @test sum(length, ibs) == L

    ibs = collect(IBSIterator(ARG, anc.genome.mutation, multiple_hits = :as_one))
    @test length(ibs) > 0
    @test all(seg -> length(seg) > 0, ibs)
    @test sum(length, ibs) == L

    ibs = collect(IBSIterator(ARG, anc.genome.mutation, multiple_hits = :JCcorrect))
    @test length(ibs) > 0
    @test all(seg -> length(seg) > 0, ibs)
    @test sum(length, ibs) == L


    indv2 = anc.alives[1][2]
    ARGmulti = collect(get_ARGsegments(anc, [indv[1], indv[2], indv2[1]]))
    @test length(ARGmulti) > 0
    @test all(seg -> length(seg) >= 0, ARGmulti)
    @test sum(length, ARGmulti) == L
    @test all(seg -> iscoalescent(seg), ARGmulti)
    @test all(seg -> 0.0 <= timespan(seg) < Inf, ARGmulti)

    mutARGmulti = collect(IBMIterator(ARGmulti, anc.genome.mutation))
    @test length(mutARGmulti) >= 0
    @test all(seg -> length(seg) >= 0, mutARGmulti)
    @test sum(length, mutARGmulti) == L
    @test all(seg -> iscoalescent(seg), mutARGmulti)
    @test all(seg -> 0.0 <= timespan(seg) < Inf, mutARGmulti)
    @test typeof(mutARGmulti[1]) == APop.ARGsegment{Int64, APop.CoalescentTree{Vector{APop.MutatedBranch}, Float64}}

    IBSmutARGmulti = collect(IBSIterator(mutARGmulti, 1,2))
    @test length(IBSmutARGmulti) > 0
    @test all(seg -> length(seg) >= 0, IBSmutARGmulti)
    @test sum(length, IBSmutARGmulti) == L
end




@testitem "WrightFisher - Events" begin
    
    population_size = 100
    mutation_rate = 2e-8
    recombination_rate = 1e-8
    L = 1_000_000_000

    
    d = Demography()
    add_population!(d, Population(id = "pop1", description = "Population 1", size = population_size))
    add_population!(d, Population(id = "pop2a", description = "Population 2A", size = 200))
    add_population!(d, Population(id = "pop2b", description = "Population 2B", size = population_size))
    add_population!(d, Population(id = "pop3", description = "Population 3", size = population_size))
    
    add_event!(d, PopulationSizeEvent(2, "pop1",  150))
    add_event!(d, PopulationSplitEvent(4, "pop1", "pop2a", "pop2b"))
    add_event!(d, PopulationMergeEvent(8, "pop2a", "pop2b", "pop3"))
    
    set_start_time!(d, 0)
    set_end_time!(d, 12)
    

    @test startswith(APop.summary(d), "Demography with 4 populations, 3 events, start time 0, end time 12")

    println(APop.summary(d))


    g = Genome(UniformRate(recombination_rate), UniformRate(mutation_rate),  L)

    model = WrightFisher()


    anc = sim_ancestry(model, d, g)
    @test length(anc.alives) == length(d.populations)


end



@testitem "WrightFisher - Migration Sampler" begin
    
    population_size = 100
    mutation_rate = 2e-8
    recombination_rate = 1e-8
    L = 1_000_000_000

    
    d = Demography()
    add_population!(d, Population(id = "pop1", description = "Population 1", size = population_size))
    add_population!(d, Population(id = "pop2a", description = "Population 2A", size = 200))
    add_population!(d, Population(id = "pop2b", description = "Population 2B", size = population_size))
    add_population!(d, Population(id = "pop3a", description = "Population 3A", size = population_size))
    add_population!(d, Population(id = "pop3b", description = "Population 3B", size = population_size))

    set_migration!(d, "pop2a", "pop2b", 0.1)
    set_migration!(d, "pop3a", "pop3b", 0.5)
    set_migration!(d, "pop3b", "pop3a", 0.5)

    set_start_time!(d, 0)
    set_end_time!(d, 12)
    
    s = APop.get_migration_parent_pop_sampler(d)
    @test length(s) == length(d.populations)

    i = 1
    pp = map(k -> s[i](), 1:10000) 
    @test all(==(1), pp)

    i = 2
    pp = map(k -> s[i](), 1:10000) 
    @test all(==(2), pp)

    i = 3
    pp = map(k -> s[i](), 1:10000) 
    @test all(in([2,3]), pp)
    @test length(pp) * 0.15 > sum(==(2), pp) > 0.05 * length(pp)
    @test sum(==(3), pp) > 0.8 * length(pp)

    i = 4
    pp = map(k -> s[i](), 1:10000) 
    @test all(in([4,5]), pp)
    @test length(pp) * 0.6 > sum(==(4), pp) > 0.4 * length(pp)
    @test length(pp) * 0.6 > sum(==(5), pp) > 0.4 * length(pp)

    i = 5
    pp = map(k -> s[i](), 1:10000) 
    @test all(in([4,5]), pp)
    @test length(pp) * 0.6 > sum(==(4), pp) > 0.4 * length(pp)
    @test length(pp) * 0.6 > sum(==(5), pp) > 0.4 * length(pp)
end

@testitem "Hudson StatefulWithDefaultIterator" begin
    using APop.HudsonModel

    si = HudsonModel.StatefulWithDefaultIterator([1, 2, 3], 0)
    @test HudsonModel.nextitem(si) == 1
    @test HudsonModel.nextitem(si) == 2
    @test HudsonModel.nextitem(si) == 3
    @test HudsonModel.nextitem(si) == 0
    @test HudsonModel.nextitem(si) == 0
    @test HudsonModel.nextitem(si) == 0
    @test HudsonModel.nextitem(si) == 0

    v = [Segment(1,10), Segment(20,30)]
    si = HudsonModel.StatefulWithDefaultIterator(v, 400)
    @test HudsonModel.nextitem(si) == (1, 10)
    @test HudsonModel.nextitem(si) == (20, 30)
    @test HudsonModel.nextitem(si) == (400, 400)
    @test HudsonModel.nextitem(si) == (400, 400)
    @test HudsonModel.nextitem(si) == (400, 400)
    @test HudsonModel.nextitem(si) == (400, 400)
end




@testitem "Hudson distribute" begin
    using APop.HudsonModel
    bps = [10, 20, 30, 100, 110, 120, 130, 140, 150, 160]

    vi = [Segment(2, 4)]
    v1, v2 = HudsonModel.distribute(vi, bps)
    @test length(v1) == 1 && length(v2) == 0

    vi = [Segment(1, 1)]
    v1, v2 = HudsonModel.distribute(vi, bps)
    @test length(v1) == 1 && length(v2) == 0


    vi = [Segment(172, 175)]
    v1, v2 = HudsonModel.distribute(vi, bps)
    @test length(v1) + length(v2) == 1


    vi = [Segment(12, 14)]
    v1, v2 = HudsonModel.distribute(vi, bps)
    @test length(v1) == 0 && length(v2) == 1

    vi = [Segment(2, 4), Segment(12, 14)]
    v1, v2 = HudsonModel.distribute(vi, bps)
    @test length(v1) == 1 && length(v2) == 1


    vi = [Segment(2, 10)]
    v1, v2 = HudsonModel.distribute(vi, bps)
    @test length(v1) == 1 && length(v2) == 0

    vi = [Segment(11, 14)]
    v1, v2 = HudsonModel.distribute(vi, bps)
    @test length(v1) == 0 && length(v2) == 1

    vi = [Segment(2, 10), Segment(11, 14)]
    v1, v2 = HudsonModel.distribute(vi, bps)
    @test length(v1) == 1 && length(v2) == 1

    vi = [Segment(2, 14)]
    v1, v2 = HudsonModel.distribute(vi, bps)
    @test length(v1) == 1 && length(v2) == 1


    vi = [Segment(2, 24)]
    v1, v2 = HudsonModel.distribute(vi, bps)
    @test length(v1) == 2 && length(v2) == 1

    vi = [Segment(2, 34)]
    v1, v2 = HudsonModel.distribute(vi, bps)
    @test length(v1) == 2 && length(v2) == 2

    vi = Segment{Int}[]
    v1, v2 = HudsonModel.distribute(vi, bps)
    @test length(v1) == 0 && length(v2) == 0

    vi = [Segment(2, 34)]
    v1, v2 = HudsonModel.distribute(vi, Int[])
    @test length(v1) == 1 && length(v2) == 0

end


@testitem "Hudson distribute rate" begin
    using APop.HudsonModel

    vi = [Segment(1, 1000000)]
    v1, v2 = HudsonModel.distribute(vi, UniformRate(1e-3))

    @test length(v1) > 0
    @test length(v2) > 0
    @test sum(length, v1) + sum(length, v2) == 1000000

    vi = [Segment(1, 1000000), Segment(2000001, 3000000)]
    v1, v2 = HudsonModel.distribute(vi, UniformRate(1e-3))

    @test length(v1) > 0
    @test length(v2) > 0
    @test sum(length, v1) + sum(length, v2) == 2000000
end


@testitem "Hudson coalesce" begin
    using APop.HudsonModel

    n = 2
    t = -1
    tmax = 0
    vc = Vector{ARGsegment{Int, CoalescentTreeTwoLineages}}()
    v1 = [Segment(1, 2), Segment(3, 4)]
    v2 = [Segment(5, 6), Segment(7, 8)]

    v = HudsonModel.coalesce(v1, v2, vc, t, tmax, n)
    @test length(v) == 4 && length(vc) == 0

    v1 = [Segment(1, 2), Segment(3, 4)]
    v2 = [Segment(15, 16), Segment(17, 18)]

    v = HudsonModel.coalesce(v1, v2, empty!(vc), t, tmax, n)
    @test length(v) == 4 && length(vc) == 0
    @test vcat(v1, v2) == v
    
    v1 = [Segment(1, 20)]
    v2 = [Segment(4,6)] 
    v = HudsonModel.coalesce(v1, v2, empty!(vc), t, tmax, n)
    @test length(v) == 2 && length(vc) == 1
    @test first(vc[1]) == 4
    @test last(vc[1]) == 6
    @test timespan(vc[1]) ==  tmax - t

    v = HudsonModel.coalesce(v2, v1, empty!(vc), t, tmax, n)
    @test length(v) == 2 && length(vc) == 1


    v1 = [Segment(1, 20)]
    v2 = [Segment(4,4)] 
    v = HudsonModel.coalesce(v1, v2, empty!(vc), t, tmax, n)
    @test length(v) == 2 && length(vc) == 1

    v = HudsonModel.coalesce(v2, v1, empty!(vc), t, tmax, n)
    @test length(v) == 2 && length(vc) == 1

    v1 = [Segment(1, 2), Segment(3, 4)]
    v2 = [Segment(2, 3), Segment(4, 5)]
    v = HudsonModel.coalesce(v1, v2, empty!(vc), t, tmax, n)
    @test length(v) == 2 && length(vc) == 3

    v1 = [Segment(1, 2), Segment(3, 4)]
    v = HudsonModel.coalesce(v1, v1, empty!(vc), t, tmax, n)
    @test length(v) == 0 && length(vc) == 2

    v1 = [Segment(1,1), Segment(3, 3)]
    v = HudsonModel.coalesce(v1, v1, empty!(vc), t, tmax, n)
    @test length(v) == 0 && length(vc) == 2
end



@testitem "coalesce loop" begin
    using APop.HudsonModel
    t = 1
    tmax = 2
    n = 2
    
    a1 = 100
    for l1 in 1:5
        for l2 in 1:5
            for a2 in a1-l2-1 : a1 + l1 + 1
                v1 = [Segment(a1, a1 + l1)]
                v2 = [Segment(a2, a2 + l2)]
                vc = Vector{ARGsegment{Int, CoalescentTreeTwoLineages}}()
                v = HudsonModel.coalesce(v1, v2, vc, t, tmax, n)
                @test length(vc) <= 1
                @test length(v) <= 2
                @test length(v) + length(vc) >= 1
            end
        end
    end
end



@testitem "Hudson distribute & coalesce" begin
    using APop.HudsonModel

    v1 = [Segment(1, 100)]
    v2 = [Segment(1, 100)]
    vc = Vector{ARGsegment{Int, CoalescentTreeTwoLineages}}()
    n = 2
    tmax = 3

    bps = [30,50,70,90]
    v11, v12 = HudsonModel.distribute(v1, bps)

    bps = [10,50,60]
    v21, v22 = HudsonModel.distribute(v2, bps)

    t = 1
    v3 = HudsonModel.coalesce(v11, v21, vc, t, tmax, n)
    v4 = HudsonModel.coalesce(v12, v22, vc, t, tmax, n)

    t = 2
    v5 = HudsonModel.coalesce(v3, v4, vc, t, tmax, n)

    sort!(vc, by = first)
    @test length(v5) == 0
    @test sum(length, vc) == 100

end


@testitem "Hudson StationaryPopulation" begin
    using APop.HudsonModel


    for genome_length in [1000, 10000, 1000000],
            population_size in [100, 1000],
            recombination_rate in [1.0e-8, 1.0e-9],
            mutation_rate in [1.0e-9, 1.0e-10]
            

    
        d = Demography()
        add_population!(d, Population(id = "pop1", size = population_size))
        set_end_time!(d, 4000)
    
        g = Genome(UniformRate(recombination_rate), UniformRate(mutation_rate),  genome_length)

        model = Hudson()

        anc = sim_ancestry(model, d, g, 2)
        ibds = APop.HudsonModel.get_ARGsegments(anc) 

        @test sum(length, ibds) == genome_length
        @test sum(length, IBSIterator(ibds, mutation(g))) == genome_length

    end
end


@testitem "Hudson StationaryPopulation Multi" begin
    using APop.HudsonModel


    for genome_length in [1000, 10000, 1000000],
            population_size in [100, 1000],
            recombination_rate in [1.0e-8, 1.0e-9],
            mutation_rate in [1.0e-9, 1.0e-10],
            n in [3,4]

    
        d = Demography()
        add_population!(d, Population(id = "pop1", size = population_size))
        set_end_time!(d, 4000)
    
        g = Genome(UniformRate(recombination_rate), UniformRate(mutation_rate),  genome_length)

        model = Hudson()

        anc = sim_ancestry(model, d, g, n)
        ibds = APop.HudsonModel.get_ARGsegments(anc) 

        @test sum(length, ibds) == genome_length
        ibm = collect(IBMIterator(ibds, g.mutation))
        @test sum(length, ibm) == genome_length
        ibs = IBSIterator(ibm, 1,2)
        @test sum(length, ibs) == genome_length

    end
end

