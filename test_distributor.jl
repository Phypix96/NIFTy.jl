using NIFTy

dom = RGDomain(500,500, harmonic = true)
pdom = PowerDomain(dom)

f1 = randn(dom)
f2 = randn(dom)
kernel1= Field(pdom, map(x -> 1/(1+x)^1.4, pdom._kvec))
kernel2 = Field(pdom, map(x -> 1/(1+x)^3, pdom._kvec))
full_field = MultiField(MultiField(f1, kernel1), MultiField(f2, kernel2))

distributor = DistributionOperator((dom, pdom), dom, pdom._pindices)
op1 = exp(hartley(distributor))
op2 = scale(target(op1), 80)(op1)
full_op = op1 + op2

using Plots
heatmap(full_op(full_field))

