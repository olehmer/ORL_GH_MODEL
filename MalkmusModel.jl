#=
Owen Lehmer - 7/10/2016

This model was originally developed to calculate the greenhouse effect felt
on early Io, to see if it could enter a runaway state from the heat of Jupiter's
accretion. It is based on the description of the Malkmus model from:

    Principles of Planetary Climates, by R.T. Pierrehumbert (2010)
    The 3rd edition was used for this model (2014) and all page/equation number
    references pertain to that edition. Throughout the code the book will be
    referred to as PPC. Most of the equations are from chapter 4.4
=#

module Malkmus
using Globals

function layer_absorption(layer::LAYER, g)
    """
    This is the top level function that will generate the absorption for all
    bands in the layer
    """

    for (i,band) in enumerate(layer.bands)
        ls = weighted_path_strong(g, layer.p1, layer.p2, layer.T, band)
        p = (layer.p1 + layer.p2)/2
        band_sum = malkmus_equivalent_width_sum(band.R2, band.S, p, band.p0, ls)
        delta = (band.end_wn - band.start_wn)*100 #convert from cm-1 to m-1
        trans = goody_random_overlap_transmittance(delta, band_sum)
        
        layer.abs[i] = 1 - trans #store the absorption
    end


end


function malkmus_equivalent_width_sum(R2, S, p2, p0, ls)
    """
    This will compute the two-parameter fit of the equivalent widths over the
    band in question. This is equation 4.71 of PPC, page 233.

    Parameters:
        R2 - R is a fitted parameter, R2 is R^2
        S - another fitted parameter
        p2 - pressure at top of the layer
        p0 - refence pressure
        ls - weighted path in the strong case
    """

    band_sum = 2*(R2/S)*(p2/p0)*(sqrt(1+S^2/R2*(p0/p2)^2*ls)-1)
    return band_sum
end


function weighted_path_strong(g,p1,p2,T,band::BAND)
    """
    The weighted path length for the strong absoprtion limit. This is equation
    4.67 on page 229 of PPC.

    For this equation we've made the assumption that our layer is isothermal
    and the mixing ratio is constant within the layer as well. It is assumed 
    that the light is directly overhead so cos(theta) from equation 4.67 is 1.

    Parameters:
        g - gravity [m/s]
        p1 - pressure at bottom of layer
        p2 - pressure at top of layer
        T - temperature of layer
        band - the band to compute on

    Returns:
        ls - the weighted path of the strong case
    """

    S_T = linestrength(T, band)
    S_T0 = linestrength(band.T0, band)

    ls = 1/g*S_T/S_T0*(band.T0/T)^band.n/band.p0*0.5*(p1^2-p2^2)

    return ls
end


function linestrength(T, band::BAND)
    """
    Calculates the line strength of the given band at temperature T
    """

    S = band.S0*(T/band.T0)^band.n*exp(H*C/KB*band.hvl*(1/T-1/band.T0))
end


function goody_random_overlap_transmittance(delta, weight_sum)
    """
    Calculates the transmittance of the band in question assuming the lines in
    the band are not correlated.

    See equation 4.70 on page 231 of PPC

    Parameters:
        delta - the width of the band
        weight_sum - the Malkmus weighted sum

    Returns:
        t - the transmittance of the band
    """


    D = 1.66 #diffusivity factor, to account for our 1D model - see p.76 of 
             #David's book for a description
    t = exp(-1/delta*weight_sum*D)
    return t
end

export layer_absorption
end #end module



