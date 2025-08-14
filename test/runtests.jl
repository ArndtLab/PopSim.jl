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



