module ClassicalMap

using LinearAlgebra, StaticArrays

export map, tangent_map

@inline map(R::SVector{3,Float64}, params::NamedTuple) =
    map(R, params.p, params.k, params.kprime)

@inline function map(R::SVector{3,Float64}, p::Real, k::Real, kprime::Real)
    X, Y, Z = R

    # Free rotation about the y-axis
    cp, sp = cos(p), sin(p)
    X3 = X * cp + Z * sp
    Y3 = Y
    Z3 = Z * cp - X * sp

    # First nonlinear kick: rotation about z
    ckz, skz = cos(k * Z3), sin(k * Z3)
    X2 = X3 * ckz - Y3 * skz
    Y2 = Y3 * ckz + X3 * skz
    Z2 = Z3

    # Second nonlinear kick: rotation about x
    ckp, skp = cos(kprime * X2), sin(kprime * X2)
    X1 = X2
    Y1 = Y2 * ckp - Z2 * skp
    Z1 = Z2 * ckp + Y2 * skp

    return SVector{3,Float64}(X1, Y1, Z1)
end

@inline tangent_map(R::SVector{3,Float64}, params::NamedTuple) =
    tangent_map(R, params.p, params.k, params.kprime)

@inline function tangent_map(
    R::SVector{3,Float64},
    p::Real,
    k::Real,
    kprime::Real
)
    X, Y, Z = R

    # Free rotation about y
    cp, sp = cos(p), sin(p)
    X3 = X * cp + Z * sp
    Y3 = Y
    Z3 = Z * cp - X * sp

    Jp = @SMatrix [
        cp   0   sp
        0    1   0
       -sp   0   cp
    ]

    # First nonlinear kick: z rotation
    ckz, skz = cos(k * Z3), sin(k * Z3)

    X2 = X3 * ckz - Y3 * skz
    Y2 = Y3 * ckz + X3 * skz
    Z2 = Z3

    Jz = @SMatrix [
        ckz   -skz   -k * (X3 * skz + Y3 * ckz)
        skz    ckz    k * (X3 * ckz - Y3 * skz)
        0      0      1
    ]

    # Second nonlinear kick: x rotation
    ckp, skp = cos(kprime * X2), sin(kprime * X2)

    Jx = @SMatrix [
        1    0     0
       -kprime * (Y2 * skp + Z2 * ckp)   ckp   -skp
        kprime * (Y2 * ckp - Z2 * skp)   skp    ckp
    ]

    return Jx * Jz * Jp
end

end
