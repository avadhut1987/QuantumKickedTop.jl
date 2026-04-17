module FloquetSystem
    using QuantumOptics
    export floquet

    function floquet(s, p, k, kp)
        # p = 2π
        b = SpinBasis(s)
        Jx = dense(0.5 * sigmax(b))
        Jy = dense(0.5 * sigmay(b))
        Jz = dense(0.5 * sigmaz(b))

        return exp(-1im * kp * Jx * Jx / (2s)) *
           exp(-1im * k * Jz * Jz / (2s)) *
           exp(-1im * p * Jy)
    end
end