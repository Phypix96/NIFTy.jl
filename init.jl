using Revise, OhMyREPL, BenchmarkTools
using Pkg
Pkg.activate(".")
using NIFTy

dom = RGDomain(10,10)
f = randn(size(dom))
