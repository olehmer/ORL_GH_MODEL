#=
Owen Lehmer 7/12/16

This file will run the greenhouse model calculations


=#

include("globals.jl")
include("MalkmusModel.jl") 

using PyPlot
using Malkmus
using Globals

function run_model()
    bands = read_hitran_data_h2o() #get all the absorption data into a layer

    num_layers = 100
    pressure_top = 1000 #pressure at the top of model [Pa]
    initial_temp = 300 #initial surface temp [K]
    Jupiter_radius = 69911513 #[m]
    Jupiter_temp = 1000 #[k]
    orbital_dist = 4.22E8 #[m]
    iteration_limit = 200 #hopefully it'll converge before it hits this limit

    pressure_bot = cc_relation(t=initial_temp)

    layers::Vector{LAYER} = []

    #create a log spacing of the pressures so we'll use smaller layers in high
    #pressure levels, +1 pressure to account for top of last layer
    pressures = get_pressure_array(pressure_bot, pressure_top, num_layers)

    for i in 1:num_layers
        l = LAYER(length(bands))
        l.T = initial_temp
        l.p1 = pressures[i]
        l.p2 = pressures[i+1]
        l.bands = bands #just a reference to the bands array
        push!(layers, l)
    end

    get_band_toa_flux(orbital_dist,Jupiter_radius,Jupiter_temp,
                      bands[1].start_wn, bands[end].end_wn, bands)

    #begin radiative loop!
    count = 0
    stop = false #set to true if temp has converged
    while count < iteration_limit && stop == false
        count += 1



end


function get_band_toa_flux(dist, rad, temp, v1, v2, bands)
    """
    Get the flux in [W m-2] for each band for inputs:

    dist - the distance in [m] from source to planet
    rad - the radius [m] of the source
    temp - the temperature of the source
    v1 - the starting wavenumber for Planck's function
    v2 - the ending wavenumber for Planck's function
    bands - the band array to store the flux at top of atmosphere in

    NOTE: v1 < v2
    """

    Bv = planck_function(temp,v1,v2)

    area = 4*pi*rad^2
    for band in bands
        flux_sum = sum(Bv[floor(band.start_wn):floor(band.end_wn)]) #get the band flux
        toa_flux = flux_sum*area/(4*pi*dist^2)
        band.toa_flux = toa_flux
    end
end







function get_pressure_array(p_bot, p_top, num_layers)
    """
    Get the log spaced pressure array from pressure p_bot to pressure p_top
    """
    ps = p_bot + p_top - logspace(log10(p_top), log10(p_bot),num_layers+1)
    return ps
end



function read_hitran_data_h2o(wn_per_band=100)
    """
    This function will read the data from the HITRAN database for a pure H20
    atmosphere for wavenumbers between 1 and 7000. It will read values for:
            1. wavenumber
            2. line strength
            3. line width
            4. temperature scaling parameter n
            5. lower state energy

    Returns:
        A layer object with all the bands initialized
    """

    hitran_data = readdlm("HITRAN/HITRAN_H20_WN1_to_WN7000.txt",',', skipstart=1)
    num_bands = Int(floor(size(hitran_data)[1]/wn_per_band))

    hitran_wns = hitran_data[:,2] #cm-1
    hitran_S = hitran_data[:,3]/(18*1.66E-27)/100 #convert from cm/molecule to m/kg
    hitran_n = hitran_data[:,4] #dimensionless
    hitran_gamma = hitran_data[:,5]*(9.869E-4) #convert from [cm-1 atm-1] to [m-1 Pa-1]
    hitran_vl = hitran_data[:,6]


    bands::Vector{BAND} = []

    for bn in 1:num_bands
        band = BAND()
        start_wn = (bn-1)*wn_per_band + 1
        end_wn = bn*wn_per_band
        band.start_wn = hitran_wns[start_wn]
        band.end_wn = hitran_wns[end_wn] 
        band.p0 = 101325 #1 atm in [Pa]
        band.T0 = 296 #reference temp in [K]

        gamma0_sum = 0
        S_sum = 0
        R2_sum = 0
        hvl_sum = 0

        for wn in start_wn:end_wn
            S_sum += hitran_S[wn]
            gamma0_sum += hitran_gamma[wn]
            R2_sum += hitran_S[wn]*hitran_gamma[wn]
            hvl_sum += hitran_vl[wn]
        end

        band.hvl = hvl_sum/wn_per_band
        band.S0 = S_sum/wn_per_band
        band.S = S_sum
        band.gamma0 = gamma0_sum/wn_per_band
        band.R2 = R2_sum

        band.flux_up = 0
        band.flux_down = 0
        band.trans = 0
        band.q = 1

        push!(bands, band)
            
    end


    #all the HITRAN data is loaded into the bands now

    return bands
end



run_model()
