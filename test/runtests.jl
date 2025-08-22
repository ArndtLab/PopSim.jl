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

end



@testitem "Genome" begin
    ur = UniformRate(0.1)
    g = Genome(ur, ur, 1000)
    @test length(g) == 1000
    @test rate(recombination(g)) == 0.1
    @test rate(mutation(g)) == 0.1
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
 
    @test sum(length, ARG) == L


    ibs = collect(APop.IBSIteratorNonMutated(ARG, anc.genome.mutation))
    @test length(ibs) > 0
    @test all(seg -> length(seg) > 0, ibs)
    @test sum(length, ibs) == L



    ibs = collect(APop.IBSIteratorNonMutated(ARG, anc.genome.mutation, multiple_hits = :as_one))
    @test length(ibs) > 0
    @test all(seg -> length(seg) > 0, ibs)
    @test sum(length, ibs) == L


    ibs = collect(APop.IBSIteratorNonMutated(ARG, anc.genome.mutation, multiple_hits = :JCcorrect))
    @test length(ibs) > 0
    @test all(seg -> length(seg) > 0, ibs)
    @test sum(length, ibs) == L
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



    g = Genome(UniformRate(recombination_rate), UniformRate(mutation_rate),  L)

    model = WrightFisher()


    anc = sim_ancestry(model, d, g)
    @test length(anc.alives) == length(d.populations)


end