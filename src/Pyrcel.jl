module Pyrcel

using PyCall

struct Aerosol
    N::Float64      # Number concentration [cm-3]
    μ::Float64      # Geometric mean diameter [μm]
    σ::Float64      # Geometric standard deviation [-]
    κ::Float64      # Hygroscopicity parameter [-]
    bins::Int       # Number of bins [-]
    label::String   # Name of the mode
end

struct Model
    T::Float64      # Temperature [K]
    p::Float64      # Pressure [Pa]
    S::Float64      # Initial mixing ratio [-]
    α::Float64      # Accomodation coefficient [-]
    w::Float64      # Updraft velocity [m/s]
end

struct SizeDistribution
    A::Any                        # Input parameters [[N1,Dg1,σg1], ...] or DMA
    De::Vector{<:AbstractFloat}   # bin edges
    Dp::Vector{<:AbstractFloat}   # bin midpoints
    ΔlnD::Vector{<:AbstractFloat} # ΔlnD of the grid
    S::Vector{<:AbstractFloat}    # spectral density
    N::Vector{<:AbstractFloat}    # number concentration per bin
    form::Symbol                  # form of the size distribution [:lognormal, ....]
end

function __init__()
    py"""
    import pyrcel as pm
    """    
end

function translate_composition(c)
    distribution = py"pm.Lognorm"(mu = 0.5 * c.μ, sigma = c.σ, N = c.N)
    return py"pm.AerosolSpecies"(c.label, distribution, kappa = c.κ, bins = c.bins)
end

function run(compositions, initial)
    initial_aerosols = map(translate_composition, compositions)
    instance = py"pm.ParcelModel"(
        initial_aerosols,
        initial.w,
        initial.T,
        initial.S,
        initial.p,
        accom = initial.α,
        console = false,
    )
    
    parcel_trace, aer_out = py"model.run"(
        t_end = 300.0 / initial.w,
        output_dt = 1.0 / initial.w,
        solver = "cvode",
        output = "dataframes",
        terminate = false,
    )

    z = parcel_trace.z.values
    T = parcel_trace.T.values[3, :]
    wc = parcel_trace.wc.values .* 1000.0
    s = parcel_trace.S.values .* 100.0

    Tf = T[end]
    smax = maximum(s)
    traces = map(c -> aer_out[c.label], compositions)

    Ni = map(i -> i.Nis .* 1e-6, initial_aerosols)
    De = map(i -> 2.0 * i.rs, initial_aerosols)
    Dp = map(D -> sqrt.(D[1:end-1] .* D[2:end]), De)
    ΔlnD = map(D -> log.(D[1:end-1] ./ D[2:end]), De)
    e = 1:length(De)
    S = map(i -> Ni[i] ./ ΔlnD[i], e)
    l = map(i -> Symbol(compositions[i].label), e)
    𝕟 = map(i -> SizeDistribution([[]], De[i], Dp[i], ΔlnD[i], S[i], Ni[i], l[i]), e)

    Nds = map(1:length(initial_aerosols)) do i
        eq, kn, alpha, phi = py"pm.binned_activation"(
            smax / 100,
            Tf,
            traces[i].values[end, :],
            initial_aerosols[i],
        )
        eq * initial_aerosols[i].total_N
    end

    Nt = mapfoldl(x -> x.total_N, +, initial_aerosols)
    af = sum(Nds) ./ Nt
    CDNC = sum(Nds)
    return (
        z = z,
        T = T,
        wc = wc,
        s = s,
        smax = smax,
        Nt = Nt,
        af = af,
        CDNC = CDNC,
        Nds = Nds,
        𝕟 = 𝕟,
        traces = traces,
        parcel = parcel_trace,
    )
end

end 