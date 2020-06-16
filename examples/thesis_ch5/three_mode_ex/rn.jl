using DiffEqBiological
using LightGraphs, SimpleWeightedGraphs

function graphRandomWalk(dimX, dimY, rem_edges=[])
    grid = create_grid(dimX,dimY)
    L = LinearIndices(grid)
    graph = SimpleWeightedDiGraph(dimX*dimY)
    
    for g in grid
        for n in create_nextnodes(g)
            if (n in grid && !(Edge(L[g], L[n]) in collect(edges(graph))))
                add_edge!(graph, L[g], L[n])
            end
        end
    end

    for e in rem_edges
        rem_edge!(graph, L[e[1]], L[e[2]])
    end
    return graph
end

function create_nextnodes(g)
    up = g + CartesianIndex(0,1)
    down = g + CartesianIndex(0,-1)
    left = g + CartesianIndex(-1,0)
    right = g + CartesianIndex(1,0)

    return [up, down, left, right]
end

function create_grid(dimX, dimY)
    return CartesianIndices(collect(Iterators.product(1:dimX, 1:dimY)))
end

target1 = CartesianIndex(1,2)

rem1 = [[CartesianIndex(1,1), CartesianIndex(1,2)],
        [CartesianIndex(1,2), CartesianIndex(1,1)],
        ]

rem2 = [[CartesianIndex(2,1), CartesianIndex(1,1)],

        [CartesianIndex(2,2), CartesianIndex(2,1)],

        [CartesianIndex(1,2), CartesianIndex(2,2)],
        ]

dimX = 2 
dimY = 2 
m_grid = create_grid(dimX,dimY)
L = LinearIndices(m_grid)
graph1 = graphRandomWalk(dimX, dimY, rem1)
graph2 = graphRandomWalk(dimX, dimY, vcat(rem1, rem2))

rn_full = @empty_reaction_network
rn_mode_1 = @empty_reaction_network
rn_mode_2 = @empty_reaction_network

for node in m_grid
    addspecies!(rn_full, Symbol(:u, L[node]))
    addspecies!(rn_mode_1, Symbol(:u, L[node]))
    addspecies!(rn_mode_2, Symbol(:u, L[node]))
end
addspecies!(rn_full, Symbol(:z))
addparam!(rn_full, :rs)
addparam!(rn_mode_1, :rs)
addparam!(rn_mode_2, :rs)

addparam!(rn_full, :rm)
addparam!(rn_mode_1, :rm)
addparam!(rn_mode_2, :rm)

@reaction_func is_eq(z, val) = (z == val) ? 1.0 : 0.0

t1 = Symbol(:u, L[target1]) 

addreaction!(rn_full, :(is_eq(z, 0)*rs*$t1), :(0 --> z)) 
addreaction!(rn_full, :(is_eq(z, 1)*0.1*rs*$t1), :(0 --> z)) 

for edge in edges(graph1)
    if edge.weight != 0
        source = Symbol(:u, src(edge))
        dest = Symbol(:u, dst(edge))
        deg = sum(graph1.weights[:, src(edge)]) 
        weight = graph1.weights[dst(edge), src(edge)]
        addreaction!(rn_full, :(is_eq(z,0)*rm*$weight*(1/$deg)), :($source --> $dest))
        addreaction!(rn_mode_1, :(rm*$weight*(1/$deg)), :($source --> $dest))
        addreaction!(rn_full, :(is_eq(z,2)*rm*$weight*(1/$deg)), :($source --> $dest))
    end
end

for edge in edges(graph2)
    if edge.weight != 0
        source = Symbol(:u, src(edge))
        dest = Symbol(:u, dst(edge))
        deg = sum(graph2.weights[:, src(edge)]) 
        weight = graph2.weights[dst(edge), src(edge)]
        addreaction!(rn_full, :(is_eq(z,1)*rm*$weight*(1/$deg)), :($source --> $dest))
        addreaction!(rn_mode_2, :(rm*$weight*(1/$deg)), :($source --> $dest))
    end
end

addjumps!(rn_full)
addodes!(rn_mode_1)
addodes!(rn_mode_2)
addjumps!(rn_mode_1)
addjumps!(rn_mode_2)
