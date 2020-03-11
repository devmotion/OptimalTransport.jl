"""
    earthmover(a, b, C)

Compute the optimal transport map
```math
argmin_{P ∈ U(a, b)} ⟨C, P⟩
```
and the dual potentials using the network simplex method.

# References

Peyré, G., & Cuturi, M.. (2018). Computational Optimal Transport. [arXiv:1803.00567](https://arxiv.org/abs/1803.00567).
"""
function earthmover(
    a::AbstractVector{<:Real},
    b::AbstractVector{<:Real},
    C::AbstractMatrix{<:Real}
)
    P = SparseArrays.spzeros(Base.promote_eltype(a, b), length(a), length(b))
    earthmover!(P, a, b, C)
end

function earthmover!(P, a, b, C)
    # find a feasible initial transport map
    initialmap!(P, a, b)

    # construct the initial graph
    n = length(a)
    m = length(b)
    graph = LightGraphs.SimpleGraph(n + m)
    for ij in findall(!iszero, P)
        i, j = Tuple(ij)

        # convert to vertices
        u = i
        v = j + n

        # add new edge
        LightGraphs.add_edge!(graph, u, v)
    end

    # compute the dual potentials
    T = Base.promote_eltype(P, C)
    f = Vector{T}(undef, n)
    g = Vector{T}(undef, m)
    dualpotentials!!(f, g, C, graph)

    while true
        # check if the dual potentials are feasible
        res = findfirstinfeasible(f, g, C)
        res === nothing && break

        # obtain vertices of the edge that will be added
        i, j = res
        u = i
        v = j + n

        # check if these vertices are already connected
        path = LightGraphs.a_star(graph, v, u)

        if !isempty(path)
            # determine the largest possible increase and which edge to remove
            e = path[1]
            θmin = P[LightGraphs.dst(e), LightGraphs.src(e) - n]
            kmin = 1
            for k in 3:2:length(path)
                e = path[k]
                θ = P[LightGraphs.dst(e), LightGraphs.src(e) - n]
                if θ < θmin
                    θmin = θ
                    kmin = k
                end
            end

            # update the transport map
            P[i, j] += θmin
            for k in 1:length(path)
                e = path[k]
                if isodd(k)
                    # try to avoid floating point errors
                    if k == kmin
                        P[LightGraphs.dst(e), LightGraphs.src(e) - n] = 0
                    else
                        P[LightGraphs.dst(e), LightGraphs.src(e) - n] -= θmin
                    end
                else
                    P[LightGraphs.src(e), LightGraphs.dst(e) - n] += θmin
                end
            end

            # remove the replaced edge
            LightGraphs.rem_edge!(graph, path[kmin])
        end

        # add the new edge
        LightGraphs.add_edge!(graph, u, v)

        # recompute the dual potentials
        dualpotentials!!(f, g, C, graph)
    end

    P, f, g
end

function findfirstinfeasible(f, g, C)
    n = length(f)
    m = length(g)

    @inbounds for j in 1:m, i in 1:n
        f[i] + g[j] > C[i, j] && return (i, j)
    end

    nothing
end

function dualpotentials!!(f, g, C, graph)
    n = length(f)
    m = length(g)
    
    visited = falses(n+m)
    S = Int[]
    @inbounds for v in LightGraphs.vertices(graph)
        # do not revisit processed vertices
        visited[v] && continue

        # handle first vertex in a new tree
        if v > n
            g[v - n] = 0
        else
            f[v] = 0
        end
        push!(S, v)
        visited[v] = true

        # handle all other vertices in this tree by DFS
        while !isempty(S)
            s = S[end]
            u = 0
            for t in LightGraphs.outneighbors(graph, s)
                if !visited[t]
                    u = t
                    if t > n
                        g[t - n] = C[s, t - n] - f[s]
                    else
                        f[t] = C[t, s - n] - g[s - n]
                    end
                    break
                end
            end
            if u == 0
                pop!(S)
            else
                visited[u] = true
                push!(S, u)
            end
        end
    end

    nothing
end

function initialmap!(P, a, b)
    n = length(a)
    m = length(b)
    size(P) == (n, m) ||
        throw(DimensionMismatch("dimensions of histograms and transport map do not match"))

    # initialize row and column sums
    T = eltype(P)
    r = T(a[1])
    c = T(b[1])

    i = j = 1
    while i ≤ n && j ≤ m
        # update initial map
        t = min(r, c)
        P[i, j] = t

        # update r and c
        r -= t
        c -= t

        if iszero(r)
            i += 1
            if i ≤ n
                r = T(a[i])
            end
        end

        if iszero(c)
            j += 1
            if j ≤ m
                c = T(b[j])
            end
        end
    end

    P
end