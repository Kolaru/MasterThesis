using Compose
using GraphPlot

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
