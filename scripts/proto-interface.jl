### A Pluto.jl notebook ###
# v1.0.4

using Markdown
using InteractiveUtils

# ╔═╡ 821010c2-b1be-11f1-b1ad-0541a2e225f7
begin
	using Pkg
	Pkg.activate("..")
	using HyperdimensionalComputing
	using Graphs
	using Statistics
	using CairoMakie
	using Random
	using Clustering
	using MolecularGraph
	using Random
end

# ╔═╡ 56b71817-8941-4ef3-a8be-ffa90612b3a1
testosterone = "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C"

# ╔═╡ 4d8f8c0c-282e-4b9e-8149-34b515468085
testosterone2 = "C[C@]12CC[C@H]3[C@H]([C@@H]1CCCCCCC(CCCC)C[C@@H]2O)CCC4=CC(=O)CC[C@]34C"

# ╔═╡ 432e5df6-1d4a-459f-ad2d-2d0104842571
m1 = smilestomol(testosterone);

# ╔═╡ a0c52178-5c30-429a-9cd0-85eb7e7c4713
m2 = smilestomol("FC(F)(F)c1cc(c(cc1)C(F)(F)F)NC(=O)[C@@H]3[C@]2(CC[C@H]4[C@H]([C@@H]2CC3)CC[C@H]5NC(=O)\\C=C/[C@]45C)C");

# ╔═╡ 24b2cb1c-c283-4132-9390-17ed7fea766b
m3 = smilestomol(testosterone2);

# ╔═╡ 3218500e-6a5f-437c-bfe0-150285890f54
m1, m2, m3

# ╔═╡ 3772dab2-2e40-40a6-972e-a4c1c54047fd
g1 = m1.graph

# ╔═╡ 4917c248-d581-4d4e-83bd-38d30899aa29
g2 = m2.graph

# ╔═╡ 55a55511-ad1c-4ff5-bcd1-1e79b2d231c5
g3 = m3.graph

# ╔═╡ 1bb3f71d-1655-4182-989e-dbbbeb70a48a
abstract type AbstractGraphEncoder end

# ╔═╡ 0f15f863-4baa-497e-80d5-3812205864c2
Base.@kwdef struct GraphHD{HV} <: AbstractGraphEncoder where {HV<:AbstractHV}
	centrality::Function = pagerank
end

# ╔═╡ 64942f56-1c7f-4d6a-a7ca-6d62be2ed9f7
struct GrapHD{HV} <: AbstractGraphEncoder where {HV<:AbstractHV} end

# ╔═╡ d2c32f3b-aa81-46f3-9541-1351254a28e6
HTML(drawsvg(m1, atomindex=true)), HTML(drawsvg(m2, atomindex=true)), HTML(drawsvg(m3, atomindex=true))

# ╔═╡ 74f0827a-a20b-45fb-8a30-b48d07dd1f51
HV = BinaryHV()

# ╔═╡ f5a49308-5a7b-4717-8406-645e08b509e7
function HyperdimensionalComputing.encode(e::GraphHD{HV}, g::SimpleGraph; D=1000) where {HV<:AbstractHV}
	C = e.centrality(g)
	lvl = LevelEncoder(HV, 0.0:0.001:1; D=D)
	H_n = encode.(Ref(lvl), C)
	H_g = multiset(map(e -> encode(e, H_n), edges(g)))
	return H_n, H_g
end

# ╔═╡ 0f913307-0610-44ad-b3bb-8b11ffdcae40
function HyperdimensionalComputing.shift(HV::AbstractHV; k::Int64=1)
	k == 0 ? HV : circshift(HV.v, k)
end

# ╔═╡ 64841c69-7afd-45f0-9e72-f0802c845dc7
function HyperdimensionalComputing.encode(
	e::Graphs.SimpleEdge,
	b::Vector{<:AbstractHV};
	directed::Bool = false
)
	# TODO: Add check for directed edges
	return b[e.src] * shift(b[e.dst], directed ? 1 : 0)
end

# ╔═╡ 1a301008-e06b-4b5b-bd8f-92ce0f761c3a
function HDGL(G, H_v, R)
	# Node neighborhoods
	H_N = Vector{eltype(H_v)}(undef, nv(G))
	for n in 1:nv(G) 
		r = [[] for _ in 1:R+1]
		map(neighborhood_dists(G, n, R)) do (n,d)
			push!(r[d+1], n)
		end
		H_N[n] = map(enumerate(r)) do (i,V)
			shift(normalize(multiset([H_v[j] for j in V])), i-1)
		end |> multibind
	end
	H_m = multiset(H_N)
	H_m, H_N
end

# ╔═╡ a3ac628a-671b-428e-af51-a1f844f3de33
function hdECFP(G, H_v, R)
	H = deepcopy(H_v)
	last = H
	if R > 0
		for r in 1:R
			_, H_r = HDGL(G, last, r)
			append!(H, H_r)
			last = H_r
		end
	end
	multiset(unique(H))
end

# ╔═╡ 4ef9e913-0910-4dca-b60f-dd210df30b79
function HyperdimensionalComputing.encode(
	e::GrapHD{HV},
	g::SimpleGraph,
	H_n::Vector{HV};
	D=1_000
) where {HV<:AbstractHV}
	H_g = multiset(map(e -> encode(e, H_n; D=D), edges(g)))
	return H_g
end

# ╔═╡ 255e1767-1c65-4f19-bd87-c778dee6e3a8
function HyperdimensionalComputing.encode(
	HV::Type{<:HyperdimensionalComputing.AbstractHV}, atom::MolecularGraph.SMILESAtom; D=10_000)
	invariants = [
		:symbol, 
		#:charge,
		#:multiplicity,
		#:isotope,
		:isaromatic,
		:stereo
	]
	props = [atom[p] for p in invariants]
	return multibind(encode.(HV, props; D=D))
end

# ╔═╡ d331b045-318b-418e-b81f-48bfbac4ba36
let
	HV = TernaryHV#TernaryHV{Int64}
	

	H_g1, H_N1 = HDGL(g1, encode.(HV, [m1.vprops[i] for i in 1:nv(g1)]), 2)
	H_g2, H_N2 = HDGL(g2, encode.(HV, [m2.vprops[i] for i in 1:nv(g2)]), 2)
	H_g3, H_N3 = HDGL(g3, encode.(HV, [m3.vprops[i] for i in 1:nv(g3)]), 2)
	
	similarity(H_g1, H_g2),	similarity(H_g1, H_g3),	similarity(H_g2, H_g3)
end

# ╔═╡ f334da59-f5ba-4189-8008-3bb04ba35f0a
let
	HV = TernaryHV

	H_g1 = hdECFP(g1, encode.(HV, [m1.vprops[i] for i in 1:nv(g1)]), 4)
	H_g2 = hdECFP(g2, encode.(HV, [m2.vprops[i] for i in 1:nv(g2)]), 4)
	H_g3 = hdECFP(g3, encode.(HV, [m3.vprops[i] for i in 1:nv(g3)]), 4)
	
	similarity(H_g1, H_g2),	similarity(H_g1, H_g3),	similarity(H_g2, H_g3)
end

# ╔═╡ 37e64a10-c590-447b-a98c-2ea48f6eaeef
let
	HV = TernaryHV{Int8}
	H_n1 = encode.(HV, [m1.vprops[i] for i in 1:nv(m1.graph)])
	H_n2 = encode.(HV, [m2.vprops[i] for i in 1:nv(m2.graph)])

	H_g1 = multiset(map(e -> encode(e, H_n1), edges(m1.graph))) |> normalize
	H_g2 = multiset(map(e -> encode(e, H_n2), edges(m2.graph))) |> normalize

	similarity(H_g1, H_g2)
end

# ╔═╡ 61fe6436-d416-4baf-855b-414e880ab3d7
let
	HV = TernaryHV{Int8}
	H_n1 = encode.(HV, [m1.vprops[i] for i in 1:nv(m1.graph)])
	H_n2 = encode.(HV, [m2.vprops[i] for i in 1:nv(m2.graph)])

	H_g1 = multiset(map(e -> encode(e, H_n1), edges(m1.graph)))
	H_g2 = multiset(map(e -> encode(e, H_n2), edges(m2.graph)))

	S = reduce(hcat, map(h -> similarity(h).(H_n2), H_n1))
	heatmap(S, colorrange=(-1,1), colormap=:redgreensplit)
end

# ╔═╡ 0fecadd1-476f-43d1-b8de-da2ab6df388f
mcis_result, mcis_status = connected_mcis(m1, m2)

# ╔═╡ f6854c6c-8012-4603-bd3b-6e655701ec93
let
	nodes1 = collect(keys(mcis_result))
	# induce MCS edges from the vertex mapping
	edges1 = induced_subgraph_edges(m1.graph, nodes1)
	svg1 = html_fixed_size(m1, 250, 250,
		atomhighlight=nodes1, bondhighlight=edges1, atomindex=true)

	nodes2 = collect(values(mcis_result))
	edges2 = induced_subgraph_edges(m2.graph, nodes2)
	svg2 = html_fixed_size(m2, 250, 250,
		atomhighlight=nodes2, bondhighlight=edges2, atomindex=true)

	[svg1, svg2]
end

# ╔═╡ 7df69ce9-ad0b-47e3-b542-3b874a6cf3bd
enc = GraphHD{TernaryHV}(centrality=betweenness_centrality)

# ╔═╡ 23dd16b2-395a-4765-b7fe-04ef2b8dc682
H_n, H_g = encode(enc, g1; D=100_000)

# ╔═╡ 7e31517f-5d03-46e3-9135-a6852de64e5b
H_n2, H_g2 = encode(enc, g2; D=100_000)

# ╔═╡ 229d29fe-953f-4560-b170-6d3ab13d2006
similarity(H_g, H_g2)

# ╔═╡ 749956d2-3b9f-4355-b39d-c698954875e8
let
	S = similarity(H_n)
	hc = hclust(1 .- S)
	heatmap(S[hc.order, hc.order], colorrange=(-1,1), colormap=:redgreensplit)
end

# ╔═╡ 6349e68e-5ce1-4632-b42e-abe7494b9ef9
let
	S = similarity(H_n2)
	hc = hclust(1 .- S,linkage=:average)
	heatmap(S[hc.order, hc.order], colorrange=(-1,1), colormap=:redgreensplit)
end

# ╔═╡ 90de689c-e4c8-419b-a1bb-1600f69e0c5a
let
	S = reduce(hcat, map(h -> similarity(h).(H_n2), H_n))
	heatmap(S, colorrange=(-1,1), colormap=:redgreensplit)
end

# ╔═╡ 113cf29d-e44e-419c-8ed4-fe81da2efc8a
similarity(H_n2)

# ╔═╡ 392c44d5-de05-4331-b95a-72a438378ff9
function double_edge_swap!(g::SimpleGraph{T}, nswap::Integer=1;
                           max_tries::Integer=100 * nswap,
                           rng::AbstractRNG=Random.default_rng()) where {T}
    nswap < 0 && throw(ArgumentError("nswap must be non-negative"))
    nswap > max_tries && throw(ArgumentError("max_tries must be at least nswap"))
    nswap == 0 && return 0
    m = ne(g)
    m < 2 && throw(ArgumentError("graph needs at least two edges to swap"))
 
    # Flat edge list we keep in sync with g so sampling is O(1).
    us = Vector{T}(undef, m)
    vs = Vector{T}(undef, m)
    for (k, e) in enumerate(edges(g))
        us[k] = src(e)
        vs[k] = dst(e)
    end
 
    swaps = 0
    tries = 0
    while swaps < nswap
        tries += 1
        tries > max_tries &&
            error("max_tries ($max_tries) exceeded after $swaps of $nswap swaps")
 
        i = rand(rng, 1:m)
        j = rand(rng, 1:m)
        i == j && continue
 
        u, v = us[i], vs[i]
        x, y = us[j], vs[j]
        rand(rng, Bool) && ((x, y) = (y, x))   # random orientation of 2nd edge
 
        # New edges: u–x and v–y.
        (u == x || v == y) && continue          # would create a self-loop
        (has_edge(g, u, x) || has_edge(g, v, y)) && continue  # multi-edge
 
        rem_edge!(g, u, v)
        rem_edge!(g, x, y)
        add_edge!(g, u, x)
        add_edge!(g, v, y)
 
        us[i], vs[i] = u, x
        us[j], vs[j] = v, y
        swaps += 1
    end
    return swaps
end

# ╔═╡ 81af4ba1-6e9a-4321-ae47-c0bc8d9333b0
let
	HV = TernaryHV
	H_g1, H_N1 = HDGL(g1, encode.(HV, [m1.vprops[i] for i in 1:nv(g1)]), 2)
	g1s = deepcopy(g1); double_edge_swap!(g1s, 10)
	H_g1s, H_N1s = HDGL(g1s, encode.(HV, [m1.vprops[i] for i in 1:nv(g1s)]), 2)

	similarity(H_g1, H_g1s)
end

# ╔═╡ Cell order:
# ╠═821010c2-b1be-11f1-b1ad-0541a2e225f7
# ╠═56b71817-8941-4ef3-a8be-ffa90612b3a1
# ╠═4d8f8c0c-282e-4b9e-8149-34b515468085
# ╠═432e5df6-1d4a-459f-ad2d-2d0104842571
# ╠═a0c52178-5c30-429a-9cd0-85eb7e7c4713
# ╠═24b2cb1c-c283-4132-9390-17ed7fea766b
# ╠═3218500e-6a5f-437c-bfe0-150285890f54
# ╠═3772dab2-2e40-40a6-972e-a4c1c54047fd
# ╠═4917c248-d581-4d4e-83bd-38d30899aa29
# ╠═55a55511-ad1c-4ff5-bcd1-1e79b2d231c5
# ╠═1bb3f71d-1655-4182-989e-dbbbeb70a48a
# ╠═0f15f863-4baa-497e-80d5-3812205864c2
# ╠═64942f56-1c7f-4d6a-a7ca-6d62be2ed9f7
# ╠═64841c69-7afd-45f0-9e72-f0802c845dc7
# ╠═1a301008-e06b-4b5b-bd8f-92ce0f761c3a
# ╠═d331b045-318b-418e-b81f-48bfbac4ba36
# ╠═a3ac628a-671b-428e-af51-a1f844f3de33
# ╠═f334da59-f5ba-4189-8008-3bb04ba35f0a
# ╠═81af4ba1-6e9a-4321-ae47-c0bc8d9333b0
# ╠═d2c32f3b-aa81-46f3-9541-1351254a28e6
# ╠═74f0827a-a20b-45fb-8a30-b48d07dd1f51
# ╠═f5a49308-5a7b-4717-8406-645e08b509e7
# ╠═0f913307-0610-44ad-b3bb-8b11ffdcae40
# ╠═4ef9e913-0910-4dca-b60f-dd210df30b79
# ╠═255e1767-1c65-4f19-bd87-c778dee6e3a8
# ╠═37e64a10-c590-447b-a98c-2ea48f6eaeef
# ╠═61fe6436-d416-4baf-855b-414e880ab3d7
# ╠═0fecadd1-476f-43d1-b8de-da2ab6df388f
# ╠═f6854c6c-8012-4603-bd3b-6e655701ec93
# ╠═7df69ce9-ad0b-47e3-b542-3b874a6cf3bd
# ╠═23dd16b2-395a-4765-b7fe-04ef2b8dc682
# ╠═7e31517f-5d03-46e3-9135-a6852de64e5b
# ╠═229d29fe-953f-4560-b170-6d3ab13d2006
# ╠═749956d2-3b9f-4355-b39d-c698954875e8
# ╠═6349e68e-5ce1-4632-b42e-abe7494b9ef9
# ╠═90de689c-e4c8-419b-a1bb-1600f69e0c5a
# ╠═113cf29d-e44e-419c-8ed4-fe81da2efc8a
# ╠═392c44d5-de05-4331-b95a-72a438378ff9
