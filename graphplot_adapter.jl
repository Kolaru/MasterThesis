using Compose
using GraphPlot
import LightGraphs

# This adapter is broken for unkown reason
import GraphPlot: _nv, _ne, _vertices, _edges, _src_index, _dst_index,
                  _adjacency_matrix, _is_directed, _laplacian_matrix

_nv(g::Graph) = nv(g)
_ne(g::Graph) = ne(g)
_vertices(g::Graph) = vertices(g)
_edges(g::Graph) = edges(g)  #TODO
_src_index(e::Edge, g::Graph) = src(e)
_src_index(e::Edge) = src(e)
_dst_index(e::Edge, g::Graph) = dst(e)
_dst_index(e::Edge) = dst(e)
_adjacency_matrix(g::Graph) = adjacency_matrix(g)
_is_directed(g::Graph) = false
_laplacian_matrix(g::Graph) = laplacian_matrix(g)

function to_LG(g::Graph)
    G = LightGraphs.SimpleGraph(nv(g))
    for edge in edges(g)
        LightGraphs.add_edge!(G, LightGraphs.Edge(edge.src, edge.dst))
    end

    return G
end

function to_svg(g::Graph, filename)
    draw(SVG(filename, 16cm, 16cm), gplot(to_LG(g)))
end
