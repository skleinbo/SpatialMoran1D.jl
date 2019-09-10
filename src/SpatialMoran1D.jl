module SpatialMoran1D

    export  ETKomaP,
            ETKomaErfP,
            ETotalKoma,
            ETMRCAKoma,
            ESegSitesKoma,
            EPiKoma,
            PSegSitesMoran,
            ESegSitesMoran,
            TMoran,
            ETMoran,
            PiMoran

    using Distributions

    """
        ETKomaP(k, n)

    Expected value for the kth coalescent time under the spatial model with PBC.
    # Arguments
    * `k`: Lineages present
    * `n`: Population size
    """
    ETKomaP(k, n)   = n^3/pi * (1/(k-1)^2-1/k^2)*(1-1/k) + (1-2/pi)
    """
        ETKomaP(k, n, m)

    Expected value for the kth coalescent time under the spatial model with PBC.
    # Arguments
    * `k`: Lineages present
    * `n`: Population size
    * `m`: Tweak parameter
    """
    ETKomaP(k, n, m) = n^3/pi *(1/(k-m)^2-1/(k+1-m)^2)*((k-m)/(k-m+1))+(1-2/pi)
    """
        ETKomaErfP(k, n, C0::Real=1.)

    Expected value for the kth coalescent time under the spatial model with PBC.
    Exact solution via inverse of error function.
    # Arguments
    * `k`: Lineages present
    * `n`: Population size
    * `C0=1.0`: Initial concentration; 1.0 is fully occupied.
    """
    ETKomaErfP(k, n, C0::Real=1.) = 1/4/C0^2*(1/erfinv((k-1)/n)^2-1/erfinv(k/n)^2)*(1-1/k)*n
    """
        ETotalKoma(n)

    Expected value for total length of a genealogy under the spatial model.
    # Arguments
    * `n`: Population size
    """
    ETotalKoma(n) = sum( k->ETKomaP(k,n)*k, 2:n )
    """
        ETMRCAKoma(n)

    Expected value for time to the MRCA under the spatial model.
    # Arguments
    * `n`: Population size
    """
    ETMRCAKoma(n) = sum( k->ETKomaP(k,n), 2:n )
    """
        ESegSitesKoma(n, μ)

    Expected value for the number of segregating sites under the spatial model.
    # Arguments
    * `n`: Population size
    * `μ`: Mutation rate
    """
    ESegSitesKoma(n, μ) = μ*ETotalKoma(n)
    """
        EPiKoma(n, μ)

    Expected value for the number of pairwise differences under the spatial model.
    # Arguments
    * `n`: Population size
    * `μ`: Mutation rate
    """
    EPiKoma(n, μ) = μ/12*n^2*(n+1)

    """
        PSegSitesMoran(k, n, θ)

    Probability Mass Function of the number of segregating sites under Moran.
    # Arguments
    * `n`: Population size
    * `θ`: Mutation rate
    """
    PSegSitesMoran(k, n, θ) = sum( i->(-1)^i*binomial(n-1,i-1)*(i-1)/(θ+i-1)*(θ/(θ+i-1))^k, 2:n)

    """
        ESegSitesMoran(n, μ)

    Expected number of segregrating sites under Moran.
    # Arguments
    * `n`: Population size
    * `θ`: Mutation rate
    """
    ESegSitesMoran(n, μ) = n*μ*sum( i->1/i, 1:n-1)

    """
        TMoran(t, k, n)

    Distribution of coalescence times according to Kingman. Moran-timescale `2/n^2`
    # Arguments
    * `k`: Lineages present
    * `n`: Population size
    """
    TMoran(t, k, n) = pdf( Exponential(binomial(k,2)), t/n^2*2)

    """
        ETMoran(k, n) = n^2/2 / binomial(k,2)

    Expected value for the kth coalescence time under Moran.
    # Arguments
    * `k`: Lineages present
    * `n`: Population size
    """
    ETMoran(k, n) = n^2/2 / binomial(k,2)

    """
        PiMoran(n, θ)

    Expected number of pairwise differences under Moran.
    # Arguments
    * `θ`: Mutation rate
    * `n`: Population size
    """
    PiMoran(n, θ) = n*θ

end # module
