abstract type Frame end

abstract type BodyFrame <: Frame end
struct Body <: BodyFrame end

abstract type WindFrame <: Frame end
struct Wind <: WindFrame end

abstract type GlobalFrame <: Frame end
struct Global <: GlobalFrame end

# function rotation_matrix(::Body, ::Global, model, X, freestream)

# end

@inline function R_from_wind(::Body, α, θ)
    R = ST.@SArray [-cos(α) sin(α); sin(α) cos(α)]
end

@inline function R_from_wind(::Wind, α, θ)
    R = ST.@SArray [1 0; 0 1]
end

@inline function R_from_wind(::Global, α, θ)
    δ = θ - α
    R = ST.@SArray [-cos(δ) -sin(δ); -sin(δ) cos(δ)]
end

@inline function R_from_body(::Body, α, θ)
    R = ST.@SArray [1 0; 0 1]
end

@inline function R_from_body(::Wind, α, θ)
    R = ST.@SArray [-cos(α) sin(α); sin(α) cos(α)]
end

@inline function R_from_body(::Global, α, θ)
    R = ST.@SArray [cos(θ) -sin(θ); sin(θ) cos(θ)]
end

@inline function R_from_global(::Body, α, θ)
    R = ST.@SArray [cos(θ) sin(θ); -sin(θ) cos(θ)]
end

@inline function R_from_global(::Wind, α, θ)
    δ = θ - α
    R = ST.@SArray [-cos(δ) -sin(δ); -sin(δ) cos(δ)]
end

@inline function R_from_global(::Global, α, θ)
    R = ST.@SArray [1 0; 0 1]
end
