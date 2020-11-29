#using NIFTy
using Test

@testset "NIFTy.jl" begin
    # Write your tests here.

    @testset "Domain" begin
        #TODO put dimension as parameter in testset
        @testset "Basics" begin
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

    @testset "Field $dim-dim of $T" for dim = 1:3, T = [Int, Float64, Float32, ComplexF64]
        s = ntuple(x -> 3, dim)
        dom = RGDomain(s...)

        f = randn(T, dom)
        g = similar(f)
        @test eltype(f) == T
        for func in [size, shape, length]
            @test @eval $func($f) == $func($g) == $func($dom)
        end

        val = fill(1, s.+1)
        @test_throws AssertionError Field(dom, val)

        #TODO make val imutable if assigned to field!
        v1 = T(4)
        v2 = T(2)
        val = fill(v1, s)
        f = Field(dom, val)
        val = fill(v2, s)
        g = Field(dom, val)
        for op in [:+, :-, :*, :/]
            h = @eval broadcast($op, $f, $g)
            @test getdomain(h) == dom
            @test h.val == fill((@eval $T($op($v1, $v2))), s)
        end
        !(T <: Integer) && for op in [:exp, :log, :sin, :cos, :tan]
            h = @eval broadcast($op, $g)
            @test getdomain(h) == dom
            @test h.val == fill((@eval $T($op($v2))), s)
        end
    end
end
