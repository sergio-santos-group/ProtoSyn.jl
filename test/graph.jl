
function f(r,p)
    @test hasparent(r) == false
    @test hasparent(p) == false
    @test isparent(p,r) == false
    @test isparent(r,p) == false
    @test isparent(nothing,r) == false
    @test haschildren(r) == false
    @test haschildren(p) == false
    @test r ∉ p.children
end

@testset "parenthood" begin
    r = Residue("A", 1)
    p = Residue("B", 2)
    
    @test r isa ProtoSyn.AbstractDigraph
    @test p isa ProtoSyn.AbstractDigraph
    
    f(r,p)
    
    x = setparent!(r, p)
    @test hasparent(r) == true
    @test hasparent(p) == false
    @test haschildren(p) == true
    @test haschildren(r) == false
    @test isparent(p, r)
    @test isparent(nothing, r) == false
    @test r ∈ p.children
    @test x===r
    
    x = popparent!(r)
    @test x===r
    f(r,p)

    setparent!(r, p)
    x = popchild!(p)
    @test x===r
    f(r,p)

    @test popchild!(p)===nothing

end
