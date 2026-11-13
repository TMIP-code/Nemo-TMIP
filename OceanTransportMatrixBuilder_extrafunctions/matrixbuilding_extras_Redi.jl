using OceanTransportMatrixBuilder
import OceanTransportMatrixBuilder: vertical_diffusion_operator_sparse_entries
import OceanTransportMatrixBuilder: preallocate_sparse_entries
import OceanTransportMatrixBuilder: globalverticalfacetriadderivative
import OceanTransportMatrixBuilder: Icoord
import OceanTransportMatrixBuilder: Jcoord
import OceanTransportMatrixBuilder: horizontalindex
import OceanTransportMatrixBuilder: i₊₁
import OceanTransportMatrixBuilder: i₋₁
import OceanTransportMatrixBuilder: j₊₁
import OceanTransportMatrixBuilder: j₋₁
import OceanTransportMatrixBuilder: k₊₁
import OceanTransportMatrixBuilder: k₋₁
import OceanTransportMatrixBuilder: verticalfacearea
import OceanTransportMatrixBuilder: pushTmixingvalues!

using SparseArrays, LinearAlgebra

# ---------------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------------

"""
    redi_wetvector(field3D, indices)

Linearize a 3D field onto the wet-cell (`Lwet`) index space, matching the
ordering used by all `redi_*` sparse operators below.
"""
redi_wetvector(field3D, indices) = [field3D[indices.C[L]] for L in indices.Lwet]

"""
    redi_taper(Smag; Smax = 0.01, Swidth = 0.1Smax)

Danabasoglu & McWilliams (1995) tapering function for isopycnal slopes:
`f(S) = 0.5*(1 + tanh((Smax - S)/Swidth))`. See earlier discussion for
rationale (smooth taper vs. Cox-style hard clipping).
"""
redi_taper(Smag; Smax = 0.01, Swidth = 0.1Smax) = 0.5 * (1 + tanh((Smax - Smag) / Swidth))

"""
    redi_ldf_c2d_aht(rn_Ud, L; knn = 1)

NEMO's `ldf_c2d` scaling (see `ldftra.F90` in NEMO source code):
`aht(i,j) = (rn_Ud/2) * max(e1,e2)^knn`, using `rn_Ud` directly (no
reference-length normalization).
"""
redi_ldf_c2d_aht(rn_Ud, L; knn = 1) = 0.5 * rn_Ud * L^knn

"""
    redi_ldf_c2d_Tpoint_scale(gridmetrics, C𝑖srf)

Approximate T-point `max(e1t,e2t)` as the average of neighbouring-face
distances (NEMO's `ldf_c2d` only defines U-/V-point values for `TRA`).
"""
function redi_ldf_c2d_Tpoint_scale(gridmetrics, C𝑖srf)
    (; distance_to_neighbour_2D) = gridmetrics
    e1t = 0.5 * (distance_to_neighbour_2D[:west][C𝑖srf] + distance_to_neighbour_2D[:east][C𝑖srf])
    e2t = 0.5 * (distance_to_neighbour_2D[:south][C𝑖srf] + distance_to_neighbour_2D[:north][C𝑖srf])
    return max(e1t, e2t)
end

"""
    redi_kappa_local(gridmetrics, indices; κRedi = nothing, ldf_c2d = false, rn_Ud = nothing, knn = 1)

Return an N-vector (wet-cell space) of the local nominal Redi diffusivity:
a constant `κRedi`, or NEMO's `ldf_c2d` per-column scaling if
`ldf_c2d = true`.
"""
function redi_kappa_local(gridmetrics, indices; κRedi = nothing, ldf_c2d = false, rn_Ud = nothing, knn = 1)
    ldf_c2d && isnothing(rn_Ud) && error("ldf_c2d=true requires rn_Ud (NEMO's lateral diffusive velocity, m/s).")
    !ldf_c2d && isnothing(κRedi) && error("κRedi is required when ldf_c2d=false.")
    (; C) = indices
    if ldf_c2d
        return [redi_ldf_c2d_aht(rn_Ud, redi_ldf_c2d_Tpoint_scale(gridmetrics, horizontalindex(C[L])); knn) for L in indices.Lwet]
    else
        return fill(κRedi, length(indices.Lwet))
    end
end

"""
    redi_mixedlayer_mask(gridmetrics, indices, mlotst)

Return `Ω::BitVector` (length N), true for wet cells at or below the
local mixed layer depth `mlotst[i,j]` — the same `zt`-vs-`mlotst`
comparison used elsewhere, inverted to EXCLUDE the mixed layer.
"""
function redi_mixedlayer_mask(gridmetrics, indices, mlotst)
    (; zt) = gridmetrics
    mlotst = mlotst |> Array
    return replace(reshape(zt, 1, 1, length(zt)) .>= mlotst, missing => false)[indices.Lwet]
end

# ---------------------------------------------------------------------------
# Discrete gradient operators (Gx, Gy, Gz)
# ---------------------------------------------------------------------------

"""
    redi_horizontal_gradient_operators(gridmetrics, indices)

Build `Gx, Gy :: SparseMatrixCSC` (N×N), where `(Gx*T)[i]` and `(Gy*T)[i]`
are centered (or one-sided, at the wet-domain edge) finite-difference
estimates of `∂T/∂x`, `∂T/∂y` at cell `i`. Depends only on the grid and
wet mask — reusable across months (does not depend on `rho3D`/`mlotst`).
"""
function redi_horizontal_gradient_operators(gridmetrics, indices)
    (; gridtopology, distance_to_neighbour_2D) = gridmetrics
    (; wet3D, Lwet, Lwet3D, C, N) = indices

    𝑖s_x, 𝑗s_x, vals_x = Int[], Int[], Float64[]
    𝑖s_y, 𝑗s_y, vals_y = Int[], Int[], Float64[]

    for 𝑖 in eachindex(Lwet)
        L𝑖 = Lwet[𝑖]
        C𝑖 = C[L𝑖]
        C𝑖srf = horizontalindex(C𝑖)

        C𝑗W = i₋₁(C𝑖, gridtopology); C𝑗E = i₊₁(C𝑖, gridtopology)
        validW = !isnothing(C𝑗W) && wet3D[C𝑗W]
        validE = !isnothing(C𝑗E) && wet3D[C𝑗E]
        if validW || validE
            n = (validW ? 1 : 0) + (validE ? 1 : 0)
            if validW
                𝑗W = Lwet3D[C𝑗W]; dW = distance_to_neighbour_2D[:west][C𝑖srf]
                push!(𝑖s_x, 𝑖); push!(𝑗s_x, 𝑖);  push!(vals_x,  1/(n*dW))
                push!(𝑖s_x, 𝑖); push!(𝑗s_x, 𝑗W); push!(vals_x, -1/(n*dW))
            end
            if validE
                𝑗E = Lwet3D[C𝑗E]; dE = distance_to_neighbour_2D[:east][C𝑖srf]
                push!(𝑖s_x, 𝑖); push!(𝑗s_x, 𝑖);  push!(vals_x, -1/(n*dE))
                push!(𝑖s_x, 𝑖); push!(𝑗s_x, 𝑗E); push!(vals_x,  1/(n*dE))
            end
        end

        C𝑗S = j₋₁(C𝑖, gridtopology); C𝑗N = j₊₁(C𝑖, gridtopology)
        validS = !isnothing(C𝑗S) && wet3D[C𝑗S]
        validN = !isnothing(C𝑗N) && wet3D[C𝑗N]
        if validS || validN
            n = (validS ? 1 : 0) + (validN ? 1 : 0)
            if validS
                𝑗S = Lwet3D[C𝑗S]; dS = distance_to_neighbour_2D[:south][C𝑖srf]
                push!(𝑖s_y, 𝑖); push!(𝑗s_y, 𝑖);  push!(vals_y,  1/(n*dS))
                push!(𝑖s_y, 𝑖); push!(𝑗s_y, 𝑗S); push!(vals_y, -1/(n*dS))
            end
            if validN
                𝑗N = Lwet3D[C𝑗N]; dN = distance_to_neighbour_2D[:north][C𝑖srf]
                push!(𝑖s_y, 𝑖); push!(𝑗s_y, 𝑖);  push!(vals_y, -1/(n*dN))
                push!(𝑖s_y, 𝑖); push!(𝑗s_y, 𝑗N); push!(vals_y,  1/(n*dN))
            end
        end
    end

    Gx = sparse(𝑖s_x, 𝑗s_x, vals_x, N, N)
    Gy = sparse(𝑖s_y, 𝑗s_y, vals_y, N, N)
    return Gx, Gy
end

"""
    redi_vertical_gradient_operator(gridmetrics, indices, Ω)

Build `Gz :: SparseMatrixCSC` (N×N), where `(Gz*T)[i]` is a centered (or
one-sided) finite-difference estimate of `∂T/∂z` at cell `i`, restricted
to `Ω` (e.g. excluding the mixed layer): a cell outside `Ω` gets an
all-zero row, and a neighbor outside `Ω` is treated as unavailable
(falling back to a one-sided estimate, or zero if neither neighbor
qualifies). Must be rebuilt whenever `Ω` changes (e.g. monthly, if `Ω`
comes from `mlotst`) — unlike `Gx`/`Gy`, it is not reusable across months.
"""
function redi_vertical_gradient_operator(gridmetrics, indices, Ω)
    (; gridtopology, zt) = gridmetrics
    (; Lwet, Lwet3D, C, N) = indices

    𝑖s, 𝑗s, vals = Int[], Int[], Float64[]

    for 𝑖 in eachindex(Lwet)
        Ω[𝑖] || continue
        L𝑖 = Lwet[𝑖]
        C𝑖 = C[L𝑖]
        _, _, k = C𝑖.I

        C𝑖T = k₋₁(C𝑖, gridtopology)
        𝑖T = isnothing(C𝑖T) ? missing : Lwet3D[C𝑖T]
        validT = !ismissing(𝑖T) && Ω[𝑖T]

        C𝑖B = k₊₁(C𝑖, gridtopology)
        𝑖B = isnothing(C𝑖B) ? missing : Lwet3D[C𝑖B]
        validB = !ismissing(𝑖B) && Ω[𝑖B]

        if validT || validB
            n = (validT ? 1 : 0) + (validB ? 1 : 0)
            if validT
                _, _, kT = C𝑖T.I
                dT = abs(zt[k] - zt[kT])
                push!(𝑖s, 𝑖); push!(𝑗s, 𝑖);  push!(vals,  1/(n*dT))
                push!(𝑖s, 𝑖); push!(𝑗s, 𝑖T); push!(vals, -1/(n*dT))
            end
            if validB
                _, _, kB = C𝑖B.I
                dB = abs(zt[kB] - zt[k])
                push!(𝑖s, 𝑖); push!(𝑗s, 𝑖);  push!(vals, -1/(n*dB))
                push!(𝑖s, 𝑖); push!(𝑗s, 𝑖B); push!(vals,  1/(n*dB))
            end
        end
    end

    return sparse(𝑖s, 𝑗s, vals, N, N)
end

"""
    redi_slopes(Gx, Gy, Gz, ρ; mindRhodz = 1e-10, Smax = 0.01, Swidth = 0.1Smax)

Compute cell-centered isopycnal slope components `Sx, Sy`, their squared
magnitude `S2`, and the DM95 taper factor, from the density field `ρ`
(wet-cell vector) and the shared gradient operators. `dρdz` is
regularized to have magnitude at least `mindRhodz` (including the exact-
zero case, fixed from the earlier `sign(0)==0` bug).
"""
function redi_slopes(Gx, Gy, Gz, ρ; mindRhodz = 1e-10, Smax = 0.01, Swidth = 0.1Smax)
    dρdx = Gx * ρ
    dρdy = Gy * ρ
    dρdz = Gz * ρ
    dρdz_reg = map(x -> x == 0.0 ? mindRhodz : sign(x) * max(abs(x), mindRhodz), dρdz)
    Sx = dρdx ./ dρdz_reg
    Sy = dρdy ./ dρdz_reg
    S2 = Sx.^2 .+ Sy.^2
    taper = [redi_taper(s; Smax, Swidth) for s in sqrt.(S2)]
    return (; Sx, Sy, S2, taper)
end

# ---------------------------------------------------------------------------
# TκHredi — horizontal (isotropic) component
# ---------------------------------------------------------------------------

"""
    redi_horizontal_operator_sparse_entries(; gridmetrics, indices, ΩH, κ_local)

Isotropic horizontal diffusion with a per-column diffusivity `κ_local`
(N-vector), face value = average of the two columns' `κ_local`. Same
symmetric-under-volume-weighting construction (shared `a`, `d` per face)
as the base package's own horizontal operator.
"""
function redi_horizontal_operator_sparse_entries(; gridmetrics, indices, ΩH, κ_local)
    (; v3D, edge_length_2D, thkcello, gridtopology, distance_to_neighbour_2D) = gridmetrics
    (; wet3D, Lwet, Lwet3D, C, N) = indices

    𝑖s, 𝑗s, Tvals = preallocate_sparse_entries(8N)
    ny = size(wet3D, 2)
    κface(a, b) = 0.5 * (κ_local[a] + κ_local[b])

    for 𝑖 in eachindex(Lwet)
        ΩH[𝑖] || continue
        L𝑖 = Lwet[𝑖]
        C𝑖 = C[L𝑖]
        C𝑖srf = horizontalindex(C𝑖)
        i, j, k = C𝑖.I
        V = v3D[C𝑖]
        C𝑗W = i₋₁(C𝑖, gridtopology)
        if !isnothing(C𝑗W)
            𝑗W = Lwet3D[C𝑗W]
            if !ismissing(𝑗W) && ΩH[𝑗W]
                iW, jW, _ = C𝑗W.I
                aij = verticalfacearea(edge_length_2D, thkcello, i, j, k, :west)
                aji = verticalfacearea(edge_length_2D, thkcello, iW, jW, k, :east)
                a = min(aij, aji)
                d = distance_to_neighbour_2D[:west][C𝑖srf]
                pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗W, κface(𝑖, 𝑗W), a, d, V)
            end
        end
        C𝑗E = i₊₁(C𝑖, gridtopology)
        if !isnothing(C𝑗E)
            𝑗E = Lwet3D[C𝑗E]
            if !ismissing(𝑗E) && ΩH[𝑗E]
                iE, jE, _ = C𝑗E.I
                aij = verticalfacearea(edge_length_2D, thkcello, i, j, k, :east)
                aji = verticalfacearea(edge_length_2D, thkcello, iE, jE, k, :west)
                a = min(aij, aji)
                d = distance_to_neighbour_2D[:east][C𝑖srf]
                pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗E, κface(𝑖, 𝑗E), a, d, V)
            end
        end
        C𝑗S = j₋₁(C𝑖, gridtopology)
        if !isnothing(C𝑗S)
            𝑗S = Lwet3D[C𝑗S]
            if !ismissing(𝑗S) && ΩH[𝑗S]
                iS, jS, _ = C𝑗S.I
                aij = verticalfacearea(edge_length_2D, thkcello, i, j, k, :south)
                aji = verticalfacearea(edge_length_2D, thkcello, iS, jS, k, :north)
                a = min(aij, aji)
                d = distance_to_neighbour_2D[:south][C𝑖srf]
                pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗S, κface(𝑖, 𝑗S), a, d, V)
            end
        end
        C𝑗N = j₊₁(C𝑖, gridtopology)
        if !isnothing(C𝑗N)
            𝑗N = Lwet3D[C𝑗N]
            if !ismissing(𝑗N) && ΩH[𝑗N]
                iN, jN, _ = C𝑗N.I
                oppdir = (j == ny) ? :north : :south
                aij = verticalfacearea(edge_length_2D, thkcello, i, j, k, :north)
                aji = verticalfacearea(edge_length_2D, thkcello, iN, jN, k, oppdir)
                a = min(aij, aji)
                d = distance_to_neighbour_2D[:north][C𝑖srf]
                pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗N, κface(𝑖, 𝑗N), a, d, V)
            end
        end
    end
    return 𝑖s, 𝑗s, Tvals
end

# ---------------------------------------------------------------------------
# TκVredi — vertical component
# ---------------------------------------------------------------------------

"""
    redi_vertical_operator_sparse_entries(; gridmetrics, indices, ρ, gh2, κ_local, Smax, Swidth, mindRhodz, Ω)

The `κRedi·S²` diagonal vertical term, evaluated per interface (not
cell-centered — this covers every interior interface including those
adjacent to the surface/bottom, unlike the cell-centered `Gz`-based
approach used for the cross term). `ρ` and `gh2 = (Gx*ρ)²+(Gy*ρ)²` are
wet-cell vectors, shared with `redi_slopes`/the cross term for
consistency.
"""
function redi_vertical_operator_sparse_entries(;
        gridmetrics, indices, ρ, gh2, κ_local, Smax = 0.01, Swidth = 0.1Smax, mindRhodz = 1e-10, Ω,
    )
    (; v3D, area2D, zt, gridtopology) = gridmetrics
    (; Lwet, Lwet3D, C, N) = indices

    𝑖s, 𝑗s, Tvals = preallocate_sparse_entries(4 * N)

    for 𝑖 in eachindex(Lwet)
        Ω[𝑖] || continue
        L𝑖 = Lwet[𝑖]
        C𝑖 = C[L𝑖]
        i, j, k = C𝑖.I
        V = v3D[C𝑖]
        a = area2D[i, j]
        ρ𝑖 = ρ[𝑖]
        gh2𝑖 = gh2[𝑖]

        C𝑗B = k₊₁(C𝑖, gridtopology)
        if !isnothing(C𝑗B)
            𝑗B = Lwet3D[C𝑗B]
            if !ismissing(𝑗B) && Ω[𝑗B]
                _, _, k′ = C𝑗B.I
                d = abs(zt[k] - zt[k′])
                dρdz_reg = max(abs((ρ[𝑗B] - ρ𝑖) / d), mindRhodz)
                gh2f = 0.5 * (gh2𝑖 + gh2[𝑗B])
                S2 = gh2f / dρdz_reg^2
                f = redi_taper(sqrt(S2); Smax, Swidth)
                pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗B, κ_local[𝑖] * S2 * f, a, d, V)
            end
        end
        C𝑗T = k₋₁(C𝑖, gridtopology)
        if !isnothing(C𝑗T)
            𝑗T = Lwet3D[C𝑗T]
            if !ismissing(𝑗T) && Ω[𝑗T]
                _, _, k′ = C𝑗T.I
                d = abs(zt[k] - zt[k′])
                dρdz_reg = max(abs((ρ𝑖 - ρ[𝑗T]) / d), mindRhodz)
                gh2f = 0.5 * (gh2𝑖 + gh2[𝑗T])
                S2 = gh2f / dρdz_reg^2
                f = redi_taper(sqrt(S2); Smax, Swidth)
                pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗T, κ_local[𝑖] * S2 * f, a, d, V)
            end
        end
    end
    return 𝑖s, 𝑗s, Tvals
end

# ---------------------------------------------------------------------------
# TκVcross — vertical-horizontal cross term
# ---------------------------------------------------------------------------

"""
    redi_cross_operator(; Gx, Gy, Gz, Sx, Sy, taper, κ_local, Ω, V)

Build the horizontal-vertical cross-term correction as

    TκVcross = -Diagonal(1 ./ V) * (Gx'*Kxz*Gz + Gz'*Kxz*Gx + Gy'*Kyz*Gz + Gz'*Kyz*Gy)

with `Kxz = Diagonal(κ_local.*Sx.*taper)`, `Kyz = Diagonal(κ_local.*Sy.*taper)`.
Each `A'*K*B + B'*K*A` term is symmetric regardless of grid irregularity, so the bracketed
sum is exactly symmetric, and `Diagonal(V)*TκVcross` is therefore exactly
symmetric too.

Cells/interfaces outside `Ω` (or where `Gz`'s row is entirely zero, e.g.
no valid vertical neighbor even after Ω-filtering) contribute nothing —
`Kxz`/`Kyz` are explicitly zeroed for `.!Ω`, not just relied on to taper
away via a large regularized slope.

Note that this is a single cell-centered estimate per cell, not the
4-triad-per-face scheme, so it will not bit-match NEMO's own cross-term magnitude.
"""
function redi_cross_operator(; Gx, Gy, Gz, Sx, Sy, taper, κ_local, Ω, V)
    Kxz = κ_local .* Sx .* taper
    Kyz = κ_local .* Sy .* taper
    Kxz[.!Ω] .= 0.0
    Kyz[.!Ω] .= 0.0

    Dcross = Gx' * sparse(Diagonal(Kxz)) * Gz + Gz' * sparse(Diagonal(Kxz)) * Gx +
             Gy' * sparse(Diagonal(Kyz)) * Gz + Gz' * sparse(Diagonal(Kyz)) * Gy

    return Diagonal(-1 ./ V) * Dcross
end

# ---------------------------------------------------------------------------
# Wrapper
# ---------------------------------------------------------------------------

"""
    redi_build_operators(;
        gridmetrics, indices, rho3D, mlotst = nothing,
        κRedi = nothing, ldf_c2d = false, rn_Ud = nothing, knn = 1,
        Smax = 0.01, Swidth = 0.1Smax, mindRhodz = 1e-10,
    )

Build a self-consistent set of Redi operator pieces, returned as
`(; TκHredi, TκVredi, TκVcross)` .
"""
function redi_build_operators(;
        gridmetrics, indices, rho3D, mlotst = nothing,
        κRedi = nothing, ldf_c2d = false, rn_Ud = nothing, knn = 1,
        Smax = 0.01, Swidth = 0.1Smax, mindRhodz = 1e-10,
    )
    (; v3D) = gridmetrics
    (; N) = indices

    rho3D = rho3D |> Array
    ρ = redi_wetvector(rho3D, indices)
    any(isnan, ρ) && error("rho3D contains NaNs at wet points.")

    V = redi_wetvector(v3D, indices)
    (any(isnan, V) || any(V .<= 0)) && error("Non-positive or NaN cell volume(s) among wet cells — check gridmetrics.v3D.")

    Ω = isnothing(mlotst) ? trues(N) : redi_mixedlayer_mask(gridmetrics, indices, mlotst)
    κ_local = redi_kappa_local(gridmetrics, indices; κRedi, ldf_c2d, rn_Ud, knn)

    Gx, Gy = redi_horizontal_gradient_operators(gridmetrics, indices)
    Gz = redi_vertical_gradient_operator(gridmetrics, indices, Ω)
    (; Sx, Sy, taper) = redi_slopes(Gx, Gy, Gz, ρ; mindRhodz, Smax, Swidth)
    gh2 = (Gx * ρ).^2 .+ (Gy * ρ).^2

    @debug "Building TκHredi"
    ΩH = trues(N)
    𝑖s, 𝑗s, Tvals = redi_horizontal_operator_sparse_entries(; gridmetrics, indices, ΩH, κ_local)
    any(isnan.(Tvals)) && error("TκHredi contains NaNs.")
    TκHredi = sparse(𝑖s, 𝑗s, Tvals, N, N)

    @debug "Building TκVredi"
    𝑖s, 𝑗s, Tvals = redi_vertical_operator_sparse_entries(;
        gridmetrics, indices, ρ, gh2, κ_local, Smax, Swidth, mindRhodz, Ω,
    )
    any(isnan.(Tvals)) && error("TκVredi contains NaNs.")
    TκVredi = sparse(𝑖s, 𝑗s, Tvals, N, N)

    @debug "Building TκVcross"
    TκVcross = redi_cross_operator(; Gx, Gy, Gz, Sx, Sy, taper, κ_local, Ω, V)
    any(isnan, TκVcross.nzval) && error("TκVcross contains NaNs.")

    return (; TκHredi, TκVredi, TκVcross)
end