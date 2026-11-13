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

"""
    vertical_diffusion_operator_sparse_entries_avt(; gridmetrics, indices, avt, Ω)

Like `vertical_diffusion_operator_sparse_entries`, but instead of a single
scalar `κV`, this takes `avt`, a 3D array of diffusivities from the parent model
(same shape/indexing as `gridmetrics.v3D`), where `avt[i,j,k]` gives the
diffusivity to apply at the top interface of cell `(i,j,k)`.

For the flux between a cell at level `k` and its Bottom neighbor at level
`k+1`, that shared interface is the top interface of the `k+1` cell, so
the correct lookup is `avt` at the neighbor's index, not the cell's own
index. For the Top-neighbor direction, the shared interface is the cell's
own top interface, so `avt` is looked up at the cell's own index.
"""
function vertical_diffusion_operator_sparse_entries_avt(; gridmetrics, indices, avt, Ω)

    # Unpack model grid
    (; v3D, area2D, zt, gridtopology) = gridmetrics
    # Unpack indices
    (; wet3D, Lwet, Lwet3D, C, N) = indices

    avt = avt |> Array

    any(isnan, avt[Lwet]) && error("avt contains NaNs at wet points.")

    𝑖s, 𝑗s, Tvals = preallocate_sparse_entries(4 * N) # 2 × 2 directions for vertical diffusion

    for 𝑖 in eachindex(Lwet)
        Ω[𝑖] || continue # only continue if inside Ω
        L𝑖 = Lwet[𝑖]
        C𝑖 = C[L𝑖]
        i, j, k = C𝑖.I
        V = v3D[C𝑖]
        a = area2D[i, j]
        # From Bottom (shared interface is the TOP interface of the Bottom neighbor)
        C𝑗B = k₊₁(C𝑖, gridtopology)
        if !isnothing(C𝑗B)
            𝑗B = Lwet3D[C𝑗B]
            if !ismissing(𝑗B) && Ω[𝑗B] # only continue if inside Ω
                _, _, k′ = C𝑗B.I
                d = abs(zt[k] - zt[k′])
                κV = avt[C𝑗B] # top interface of the neighbor below = interface between i and j
                pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗B, κV, a, d, V)
            end
        end
        # From Top (shared interface is the TOP interface of the cell itself)
        C𝑗T = k₋₁(C𝑖, gridtopology)
        if !isnothing(C𝑗T)
            𝑗T = Lwet3D[C𝑗T]
            if !ismissing(𝑗T) && Ω[𝑗T] # only continue if inside Ω
                _, _, k′ = C𝑗T.I
                d = abs(zt[k] - zt[k′])
                κV = avt[C𝑖] # top interface of the cell itself = interface between i and j
                pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗T, κV, a, d, V)
            end
        end
    end
    return 𝑖s, 𝑗s, Tvals
end

"""
    buildTκVavt(; gridmetrics, indices, avt)

Build the vertical diffusivity operator directly from a full 3D field `avt`
of diffusivities output by the parent model, where `avt[i,j,k]` is the
diffusivity applied at the top interface of cell `(i,j,k)` (same
shape/indexing as `gridmetrics.v3D`).
"""
function buildTκVavt(; gridmetrics, indices, avt)

    N = indices.N

    # Entire ocean, same as buildTκVdeep
    Ω = trues(N)

    @debug "Building TκVavt"
    𝑖s, 𝑗s, Tvals = vertical_diffusion_operator_sparse_entries_avt(; gridmetrics, indices, avt, Ω)

    any(isnan.(Tvals)) && error("TκVavt contains NaNs.")

    TκVavt = sparse(𝑖s, 𝑗s, Tvals, N, N)

    return TκVavt
end
