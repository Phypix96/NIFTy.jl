import Base.eachindex

abstract type AbstractField end

struct Field <: AbstractField
    domain
    val
end

struct MultiField

end

function in(f::AbstractField)
    error()
end

eachindex(f::AbstractField) = eachindex(f.val)
