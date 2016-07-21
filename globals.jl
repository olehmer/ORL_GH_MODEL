"""
Constants and structs used throughout the code
"""

module Globals

type BAND
    start_wn::Float64
    end_wn::Float64
    n::Float64 #line-dependent exponent
    hvl::Float64 #hvl - the energy of the lower energy state in the transition
    p0::Float64 #reference pressure for band measurements
    T0::Float64 #reference temperature for band measurements
    gamma0::Float64 #the reference line width at T0
    S0::Float64 #the average linestrength at T0 
    R2::Float64 #the fitted Malkmus value R^2
    S::Float64 #the fitted Malkmus value S

    BAND() = new(0,0,0,0,0,0,0,0,0,0)
end


type LAYER
    T::Float64 #temperature of layer
    p1::Float64 #pressure at bottom of layer
    p2::Float64 #pressure at top of layer

    flux_up::Vector{Float64}
    flux_down::Vector{Float64}
    abs::Vector{Float64} #the absorption in this layer

    bands::Vector{BAND} #vector of bands

    LAYER() = new(0,0,0,[],[],[],[])
    LAYER(num_bands) = new(0,0,0,zeros(num_bands),zeros(num_bands),
                                 zeros(num_bands),[])
end


KB = 1.38E-23 #Boltzmann constant in [m2 kg s-2 K-1]
H = 6.626E-34 #Planck constant in [J s]
C = 2.99792458E8 #Speed of light [m s-1]                                          
SIGMA = 5.67E-8 #Stefan-Boltzmann constant

UP = 0
DOWN = 1

function cc_relation(;p=-1, t=-1, gas=1)
    """
    Returns the Clausius-Clapeyron temp or pressure given the other

    p - pressure in [Pa]
    t - temperature [K]
    gas - the gas in question, defaults to H2O
    """

    R = 0
    Lv = 0
    p0 = 0
    t0 = 0

    if 1 == gas #water
        R = 461.89 #gas constant for water [J K-1 kg-1]
        Lv = 2.425E6 #latent heat of water vapor [J kg-1]
        t0 = 273.15 #reference temp [K]
        p0 = 611.0 #reference pressure [Pa]
    end

    result = 0
    if p >= 0
        #we're given pressure, solve for t
        result = t0/(1-log(p/p0)*R*t0/Lv)
    elseif t >= 0
        result = p0*exp(Lv/R/t0*(1-t0/t))
    end
    return result
end

function planck_function(T, v1, v2)
    """
    Return an array of the flux for each wavenumber v between v1 and v2 at 
    temperature T. NOTE: v1 < v1

    If v1 and v2 are not integers the floor(v1) and floor(v2) will be used

    The returned flux is scaled by pi to account for the flux from a hemisphere
    of isotropic radiation.

    Parameters:
    T - the temperature of the body
    v1 - the first wavenumber to consider
    v2 - the final wavenumber to consider

    Returns:
    Bv - an array with containing v2-v1 entries (one per wavenumber), each entry
    is the flux in [W m-2] at the given wavenumber
    """

    v1 = floor(v1)
    v2 = floor(v2)
    if v2-v1 < 1
        return []
    end

    vs = linspace(v1,v2,(v2-v1)+1) 
    Bv = zeros(length(vs))

    c1 = 2*H/C^2*(100*C)^3*C*100 #[W m-2 sr-1 cm4]
    c2 = H/KB*(100*C) #[K cm]
    for i=1:length(Bv)
        Bv[i] = c1*vs[i]^3/(exp(c2*vs[i]/T)-1)
    end

    return pi*Bv #scale by pi to account for flux over hemisphere
end

export BAND, LAYER, KB, H, C, SIGMA, UP, DOWN, cc_relation, planck_function
end #end module
