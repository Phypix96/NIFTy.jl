#using NIFTy
using Test

@testset "NIFTy.jl" begin
    # Write your tests here.

    @testset "Domain" begin
        @testset "Functions" begin
            for s in ((3,), (3,3), (3,3,3))
                dom = RGDomain(s...)
                @test lengths(dom) == map(x -> 1., size(dom))
                @test distances(dom) == 1 ./ size(dom)

                for harmonic in (true, false)
                    l = ntuple(x -> 1/x, length(s))
                    dom = RGDomain(s...; lengths = l, harmonic = harmonic)
                    @test isharmonic(dom) == harmonic
                    @test size(dom) == s
                    @test shape(dom) == s
                    @test dims(dom) == length(s)
                    @test lengths(dom) == l
                    @test distances(dom) == lengths(dom) ./ size(dom)

                    @test isharmonic(getcodomain(dom)) == !harmonic
                    @test getcodomain(getcodomain(dom)) == dom
                end
            end
        end

        @testset "Domain Tuple" begin
            nelem = 3
            doms = []
            for dim = 1:3
                s = ntuple(x -> nelem, dim)
                push!(doms, RGDomain(s..., harmonic = (dim%2 == 0)))
            end
            for i in 1:length(doms)
                dt = DomainTuple(doms[1:i]...)
                @test domains(dt) == tuple(doms[1:i]...)
                @test shape(dt) == map(size, tuple(doms[1:i]...))
                @test size(dt) == ntuple(x -> nelem, div(i * (i + 1), 2))
                @test length(dt) == nelem^(div(i * (i + 1), 2))
            end
        end
    end

    @testset "Field" begin
        dom = 
        f = 
    end
end
