module Compose

using Dates
using HDF5
using Statistics: mean

#

# write(*,*) i2,' temperature T                                      [MeV]        '
# write(*,*) i2,' baryon number density n_b                          [fm^-3]      '
# write(*,*) i2,' hadronic charge fraction Y_q                       []           '
# write(*,*) i2,' magnetic field strength B                          [G]          '

const thermo_qtys = Dict(1 => "pressure p                                         [MeV fm^-3]  ",
                         2 => "entropy per baryon S                               []           ",
                         3 => "shifted baryon chemical potential mu_b-m_n         [MeV]        ",
                         4 => "charge chemical potential mu_q                     [MeV]        ",
                         5 => "lepton chemical potential mu_l                     [MeV]        ",
                         6 => "scaled free energy per baryon F/m_n-1              []           ",
                         7 => "scaled internal energy per baryon E/m_n-1          []           ",
                         8 => "scaled enthalpy energy per baryon H/m_n-1          []           ",
                         9 => "scaled free enthalpy per baryon G/m_n-1            []           ",
                         10 => "derivative dp/dn_b|E                               [MeV]        ",
                         11 => "derivative p/dE|n_b                                [fm^-3]      ",
                         12 => "square of speed of sound (c_s)^2                   []           ",
                         13 => "specific heat capacity at constant volume c_V      []           ",
                         14 => "specific heat capacity at constant pressure c_p    []           ",
                         15 => "adiabatic index Gamma                              []           ",
                         16 => "expansion coefficient at constant pressure alpha_p [MeV^-1]     ",
                         17 => "tension coefficient at constant volume beta_V      [fm^-3]      ",
                         18 => "isothermal compressibility kappa_T                 [MeV^-1 fm^3]",
                         19 => "isentropic compressibility kappa_S                 [MeV^-1 fm^3]",
                         20 => "free energy per baryon F                           [MeV]        ",
                         21 => "internal energy per baryon E                       [MeV]        ",
                         22 => "enthalpy per baryon H                              [MeV]        ",
                         23 => "free enthalpy per baryon G                         [MeV]        ",
                         24 => "energy density epsilon                             [MeV/fm^3]   ")

# Additional thermodynamic quantities

# Derivative quantities

# Pairs

# Quadruples
# write(*,*) i2,' total number fraction Y of particle set with index     ',idx_q(i1)
# write(*,*) i2,' average mass number A_av of particle set with index    ',idx_q(i1)
# write(*,*) i2,' average proton number Z_av of particle set with index  ',idx_q(i1)
# write(*,*) i2,' average neutron number N_av of particle set with index ',idx_q(i1)

# Microscopic quantities

# Error quantities

#

function clip(range::UnitRange, mask::UnitRange)
    @assert all(==(1), step(range))
    imin = max(first(range), first(mask))
    imax = min(last(range), last(mask))
    return imin:imax
end

function linterp(x0::T, y0::T, x1::T, y1::T, x::T,
                 extrap::Bool=false) where {T<:AbstractFloat}
    if !extrap
        if x0 < x1
            @assert x0 <= x && x <= x1
        elseif x0 > x1
            @assert x1 <= x && x <= x0
        else
            # We expect the input to have a non-empty range
            @assert false
            @assert x == x0
        end
    end
    y = (x1 - x) / (x1 - x0) * y0 + (x - x0) / (x1 - x0) * y1
    if !extrap
        if y0 < y1
            @assert y0 <= y && y <= y1
        elseif y0 > y1
            @assert y1 <= y && y <= y0
        else
            @assert y == y0
        end
    end
    return y
end

function linterp(x0::Number, y0::Number, x1::Number, y1::Number, x::Number,
                 extrap::Bool=false)
    return linterp(promote(x0, y0, x1, y1, x)..., extrap)
end

function linterp0(arr::Array{T,3}, i::Int, j::Int,
                  k::Int)::T where {T<:AbstractFloat}
    return arr[i, j, k]
end
function linterp1(arr::Array{T,3}, i::Int, j::Int, k::Int,
                  x::T)::T where {T<:AbstractFloat}
    return (1 - x) * linterp0(arr, i, j, k) + x * linterp0(arr, i + 1, j, k)
end
function linterp2(arr::Array{T,3}, i::Int, j::Int, k::Int, x::T,
                  y::T)::T where {T<:AbstractFloat}
    return (1 - y) * linterp1(arr, i, j, k, x) +
           y * linterp1(arr, i, j + 1, k, x)
end
function linterp3(arr::Array{T,3}, i::Int, j::Int, k::Int, x::T, y::T,
                  z::T)::T where {T<:AbstractFloat}
    @assert 0 ≤ x ≤ 1
    @assert 0 ≤ y ≤ 1
    @assert 0 ≤ z ≤ 1
    return (1 - z) * linterp2(arr, i, j, k, x, y) +
           z * linterp2(arr, i, j, k + 1, x, y)
end

#

# de/dt > 0
function find_outliers_dedt(table_e::AbstractArray{<:Real,3})
    outliers = falses(size(table_e))
    nnb, nt, nyq = size(table_e)

    for iyq in 1:nyq, inb in 1:nnb

        # Find outliers
        prev_it = 1
        it = prev_it + 1
        while it ≤ nt
            if table_e[inb, it, iyq] ≤ table_e[inb, prev_it, iyq]
                # We found an inconsistency. Either this point it or
                # the previous point prev_it is an outlier that needs
                # to be removed.
                inbrange = clip((inb - 1):(inb + 1), 1:nnb)
                itrange = clip((prev_it - 1):(it + 1), 1:nt)
                iyqrange = clip((iyq - 1):(iyq + 1), 1:nyq)
                avg_e = mean(table_e[inbrange, itrange, iyqrange])
                # Assume that the point further from the average is
                # the outlier
                if abs(table_e[inb, prev_it, iyq] - avg_e) >
                   abs(table_e[inb, it, iyq] - avg_e)
                    # Disable the previous point
                    outliers[inb, prev_it, iyq] = true
                    # Next, compare the current point to the one
                    # before the previous point
                    prev_it -= 1
                    while 1 ≤ prev_it && outliers[inb, prev_it, iyq]
                        prev_it -= 1
                    end
                    if 1 ≤ prev_it
                        # keep it
                    else
                        # All points before this are disabled: move on to
                        # the next point
                        prev_it = it
                        it += 1
                    end
                else
                    # Disable the current point
                    outliers[inb, it, iyq] = true
                    # keep prev_it
                    it += 1
                end
            else
                # There was no inconsistency. Move on to the next
                # point pair
                prev_it = it
                it += 1
            end
        end

        # Check consistency for all non-outliers

        # Find the first non-disabled point
        prev_it = 1
        while prev_it ≤ nt && outliers[inb, prev_it, iyq]
            prev_it += 1
        end
        # Check consistency for all following points
        for it in (prev_it + 1):nt
            if !outliers[inb, it, iyq]
                @assert table_e[inb, it, iyq] > table_e[inb, prev_it, iyq]
                prev_it = it
            end
        end
    end

    return outliers
end

function correct_outliers_dedt(table_e::AbstractArray{<:Real,3},
                               outliers::AbstractArray{Bool,3})
    @assert size(table_e) == size(outliers)
    table_e_new = copy(table_e)
    nnb, nt, nyq = size(table_e)

    for iyq in 1:nyq, it in 1:nt, inb in 1:nnb
        if outliers[inb, it, iyq]
            # Find two valid neighbouring points
            extrap = false
            prev_it = it - 1
            while 1 ≤ prev_it && outliers[inb, prev_it, iyq]
                prev_it -= 1
            end
            next_it = it + 1
            while next_it ≤ nt && outliers[inb, next_it, iyq]
                next_it += 1
            end
            if prev_it < 1
                extrap = true
                prev_it = next_it
                next_it += 1
                while next_it ≤ nt && outliers[inb, next_it, iyq]
                    next_it += 1
                end
            end
            if next_it > nt
                extrap = true
                next_it = prev_it
                prev_it -= 1
                while 1 ≤ prev_it && outliers[inb, prev_it, iyq]
                    prev_it -= 1
                end
            end
            @assert 1 ≤ prev_it ≤ nt
            @assert 1 ≤ next_it ≤ nt
            # Linear interpolation
            table_e_new[inb, it, iyq] = linterp(prev_it,
                                                table_e[inb, prev_it, iyq],
                                                next_it,
                                                table_e[inb, next_it, iyq], it,
                                                extrap)
        end
    end

    return table_e_new
end

#

function main()
    filename = "$(ENV["HOME"])/data/eos31/eos/eoscompose.h5"
    filename2 = "$(ENV["HOME"])/data/eos31/eos/eoscompose_corrected.h5"

    file = h5open(filename, "r")

    # Metadata
    date = Date(read(attributes(file["metadata"])["date"]), dateformat"y/m/d")
    time = Time(read(attributes(file["metadata"])["time"]),
                dateformat"H\hM\mS\ss")

    # Parameters
    pointsnb = read(attributes(file["Parameters"])["pointsnb"])[]
    pointst = read(attributes(file["Parameters"])["pointst"])[]
    pointsyq = read(attributes(file["Parameters"])["pointsyq"])[]
    tabulation_scheme = read(attributes(file["Parameters"])["tabulation_scheme"])[] # 0 explicit, 1 loop

    @assert pointsnb >= 0
    @assert pointst >= 0
    @assert pointsyq >= 0

    nb = read(file["Parameters"]["nb"])
    t = read(file["Parameters"]["t"])
    yq = read(file["Parameters"]["yq"])

    @assert length(nb) == pointsnb
    @assert length(t) == pointst
    @assert length(yq) == pointsyq

    @assert all(nb[i] > nb[i - 1] for i in 2:pointsnb)
    @assert all(t[i] > t[i - 1] for i in 2:pointst)
    @assert all(yq[i] > yq[i - 1] for i in 2:pointsyq)

    # Thermodynamic quantities
    pointsqty = read(attributes(file["Thermo_qty"])["pointsqty"])[]
    index_thermo = read(file["Thermo_qty"]["index_thermo"])
    thermo = read(file["Thermo_qty"]["thermo"])

    @assert pointsqty >= 0
    @assert length(index_thermo) == pointsqty
    @assert length(Set(index_thermo)) == length(index_thermo)
    @assert all(idx ∈ keys(thermo_qtys) for idx in index_thermo)

    # This order is hard-coded into `hdf5compose.f90`. It is different
    # from the order given in `eos.parameters`.
    @assert size(thermo) == (pointsnb, pointst, pointsyq, pointsqty)
    @info "Read EOS table with ($pointsnb, $pointst, $pointsyq) entries for $pointsqty quantities"

    inverse_index_thermo = Dict(idx => n
                                for (n, idx) in enumerate(index_thermo))

    # Ensure dE/dT > 0

    qty_nb(inb, it, iyq) = nb[inb] # baryon number density n_b
    qty_t(inb, it, iyq) = t[it]    # temperature T
    qty_yq(inb, it, iyq) = yq[iyq] # hadronic charge fraction Y_q
    ne = 21                        # internal energy per baryon E
    ie = inverse_index_thermo[ne]
    qty_e(inb, it, iyq) = thermo[inb, it, iyq, ie]

    @info "Looking for outliers..."
    table_e = @view thermo[:, :, :, ie]
    outliers = find_outliers_dedt(table_e)
    npoints = length(outliers)
    noutliers = count(outliers)
    @info "dE/dT ≤ 0 at $noutliers of $npoints points ($(100 * noutliers / npoints)%)"

    @info "Correcting outliers..."
    table_e_new = correct_outliers_dedt(table_e, outliers)
    outliers_new = find_outliers_dedt(table_e_new)
    any(outliers_new) && error("Could not correct outliers")
    @info "Corrected all outliers"

    # Save result
    @info "Writing result..."

    # Copy HDF5 file
    cp(filename, filename2; force=true)
    # Update thermodynamics table
    file2 = h5open(filename2, "r+")
    file2["Thermo_qty"]["thermo"][:, :, :, ie] = table_e
    close(file2)

    @info "Done."
    return nothing
end

end
