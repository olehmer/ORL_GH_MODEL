#=
Owen Lehmer 7/12/16

This file will run the greenhouse model calculations


=#

include("globals.jl")
include("MalkmusModel.jl") 

using PyPlot
using Malkmus
using Globals

function run_model(;
        num_layers = 100,
        pressure_top = 100, #pressure at the top of model [Pa]
        surface_T = -1, #initial surface temp [K]
        iteration_limit = 300, #hopefully it'll converge before it hits this limit
        Jupiter_temp = 1000, #[k]
        albedo = 0.1,
        ts = 10000)  #timestep for model [s]

    bands = read_hitran_data_h2o() #get all the absorption data into a layer
    bands_toa_flux = zeros(length(bands)) #flux at TOA in each band

    Jupiter_radius = 69911513 #[m]
    orbital_dist = 4.22E8 #[m] orbital distance of Io
    gravity = 1.796 #[m/s2] Io gravity
    cp = 2.02E3 #specific heat of water gas [J kg-1 K-1]

    if surface_T < 0
        #the user hasn't specified a temp, just set it to be in balance
        #with incoming flux
        flux = SIGMA*Jupiter_temp^4*(4*pi*Jupiter_radius^2)/(4*pi*orbital_dist^2)
        surface_T = (flux*(1-albedo)/(SIGMA*4))^0.25
    end
    initial_temp = surface_T


    pressure_bot = cc_relation(t=surface_T)

    if pressure_bot < pressure_top
        println("The pressure at the bottom of the atmosphere is lower than the
        upper limit")
        return
    end

    layers::Vector{LAYER} = []

    #create a log spacing of the pressures so we'll use smaller layers in high
    #pressure levels, +1 pressure to account for top of last layer
    pressures = get_pressure_array(pressure_bot, pressure_top, num_layers)

    for i in 1:num_layers
        l = LAYER(length(bands))
        l.T = surface_T 
        l.bands = bands #just a reference to the bands array
        push!(layers, l)
    end

    update_pressure(layers,pressures)

    get_band_toa_flux(orbital_dist,Jupiter_radius,Jupiter_temp,
                      bands[1].start_wn, bands[end].end_wn, 
                      bands_toa_flux,bands,albedo)

    layers[end].flux_down = bands_toa_flux #set incoming flux on top layer

    #begin radiative loop!
    count = 0
    converged = false

    while count < iteration_limit && converged == false
        count += 1

        for i=num_layers:-1:1 #go from top to bottom
            surface_T = update_layer_at_index2(layers,i,gravity,ts,cp,surface_T)
        end

        #now that we have the surface temperature, update the atmospheric
        #pressure from the CC relation
        pressure_bot = cc_relation(t=surface_T)
        if pressure_bot < pressure_top
            #=
            In this case we are in a very low pressure atmosphere, rather than
            set the surface pressure accordingly (typically this happens when
            the pressure is <<1 at the surface) just assume there is no 
            greenhouse effect and set the surface temperature to be in balance
            with the incoming flux, aka the start temp
            =#

            pressures = get_pressure_array(cc_relation(t=initial_temp), 
                                           pressure_top,
                                           num_layers)
        else
            pressures = get_pressure_array(cc_relation(t=surface_T), 
                                           pressure_top,
                                           num_layers)
        end

        #update each layer
        update_pressure(layers,pressures)


        println("COUNT=",count,", Surface Temp = ", Int(round(surface_T)), ", TOA Temp = ",Int(round(layers[end].T)))
            
    end

    #temporary code
    temps = zeros(length(layers))
    ps = zeros(length(layers))
    for i=1:length(layers)
        temps[i] = layers[i].T
        ps[i] = (layers[i].p1 + layers[i].p2)/2

        @printf("LAYER %3d  T=%3.4f, P=%0.3e\n",i,temps[i],ps[i])
    end


    plot(temps,ps)
    yscale("log")
    ylim(minimum(ps),maximum(ps))
    ylim(ylim()[end:-1:1])
    xlabel("Temperature [K]")
    ylabel("Pressure [Pa]")
end

function update_pressure(layers::Vector{LAYER}, pressures)
    for i in 1:length(layers)
        layers[i].p1 = pressures[i]
        layers[i].p2 = pressures[i+1]
    end
end


function update_layer_at_index2(layers, i, gravity, ts, cp, surface_T)

    layer_absorption(layers[i], gravity)

    absorbed_flux = layers[i].abs.*(layers[i].flux_up+layers[i].flux_down)

    planck = planck_function(layers[i].T, layers[i].bands[1].start_wn,
                             layers[i].bands[end].end_wn) #emitted flux
    emitted_flux = convert_planck_to_bands(planck, layers[i].bands)
    emitted_flux = (emitted_flux*2).*layers[i].abs

    layers[i].T -= gravity/cp*(sum(emitted_flux-absorbed_flux))/
                        (layers[i].p1-layers[i].p2)*ts

    if i < length(layers)
        layers[i+1].flux_up = 0.5*emitted_flux + 
                              (1-layers[i].abs).*layers[i].flux_up
    else
        println("TOA flux diff: ",sum(layers[i].flux_up-layers[i].flux_down))
        println("TOA temp     : ",layers[i].T)
    end

    if i > 1
        #not the bottom layer
        layers[i-1].flux_down = 0.5*emitted_flux + 
                                (1-layers[i].abs).*layers[i].flux_down
    else
        surface_emiss = planck_function(surface_T, layers[i].bands[1].start_wn,
                                     layers[i].bands[end].end_wn) #emitted flux 
        #convert the emission array to a band averaged array
        surface_emiss = convert_planck_to_bands(surface_emiss, layers[i].bands)

        surface_abs = 0.5*emitted_flux + (1-layers[i].abs).*layers[i].flux_down

        capacity = 200000 #convergence does not depend on this value
        #TODO make convergence run dependent to speed up convergence
        surface_T -= sum(surface_emiss-surface_abs)*ts/capacity

        println("surface flux diff: ",sum(surface_emiss-surface_abs))

        layers[i].flux_up = surface_emiss
    end

    return surface_T
end

function update_layer_at_index(layers, i, gravity, ts, cp, surface_T)
    layer_absorption(layers[i], gravity)

    #get the absorbed flux
    a_flux = layers[i].abs.*(layers[i].flux_up+layers[i].flux_down)
    e_flux_temp = planck_function(layers[i].T, layers[i].bands[1].start_wn,
                             layers[i].bands[end].end_wn) #emitted flux

    e_flux = convert_planck_to_bands(e_flux_temp, layers[i].bands)
    e_flux = layers[i].abs.*e_flux #Kirchoff's law

    T_new = layers[i].T - gravity/cp*((sum(e_flux)-sum(a_flux))/
            (layers[i].p1-layers[i].p2))*ts

    layers[i].T = T_new

    trans = 1 - layers[i].abs

    #set the flux on the surrounding layers
    if i < length(layers) #this isn't the top, so set the layer above
        layers[i+1].flux_up = 0.5*e_flux + trans.*layers[i].flux_up
    end

    if i > 1 #not the bottom layer
        layers[i-1].flux_down = 0.5*e_flux + trans.*layers[i].flux_down
    else
        #this is the bottom layer
        cp_liquid = 4200 #specific heat of water, [J K-1 kg-1]
        io_area = 4*pi*(1.821E6)^2
        ocean_mass = io_area*0.5*1000 #approx mass of 0.5m ocean

        s_e_flux_temp = planck_function(surface_T, layers[i].bands[1].start_wn,
                      layers[i].bands[end].end_wn) #emitted flux
        s_e_flux = convert_planck_to_bands(s_e_flux_temp, layers[i].bands)

        layer_flux_down = 0.5*e_flux + trans.*layers[i].flux_down


        println("net flux up from surface: ",sum(s_e_flux-layer_flux_down))

        delta_T = sum(layer_flux_down-s_e_flux)*io_area/ocean_mass/cp_liquid*ts

        #we assume the surface absorbs everything in the IR
        surface_T = surface_T + delta_T

        surface_T = surface_T<0?0:surface_T

        layers[i].flux_up = s_e_flux
    end
    return surface_T
end

function get_band_toa_flux(dist, rad, temp, v1, v2, bands_toa_flux, 
                           bands::Vector{BAND}, albedo)
    """
    Get the flux in [W m-2] for each band for inputs:

    dist - the distance in [m] from source to planet
    rad - the radius [m] of the source
    temp - the temperature of the source
    v1 - the starting wavenumber for Planck's function
    v2 - the ending wavenumber for Planck's function
    bands_toa_flux - the array to store the flux at top of atmosphere in
    bands - the array of bands

    NOTE: v1 < v2

    Output is stored in bands_toa_flux
    """

    Bv = planck_function(temp,v1,v2)

    area = 4*pi*rad^2
    for i in 1:length(bands_toa_flux)
        band = bands[i]
        #get the band flux sum over the hemisphere 
        flux_sum = sum(Bv[Int(floor(band.start_wn)):Int(floor(band.end_wn))]) 

        bands_toa_flux[i] = (1-albedo)*flux_sum*area/(4*pi*dist^2)/4 #div by 4 since 
                                       #spherical and rotating body
    end
end


function convert_planck_to_bands(Bv, bands)
    """
    Convert an array of planck function values to bands. Bv is assumed to contain
    the planck function for each wavelength in bands
    """
    band_flux = zeros(length(bands))
    for (i,band) in enumerate(bands)
        flux_sum = sum(Bv[Int(floor(band.start_wn)):Int(floor(band.end_wn))])
        band_flux[i] = flux_sum
    end
    return band_flux
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


        push!(bands, band)
            
    end


    #all the HITRAN data is loaded into the bands now

    return bands
end



run_model()
