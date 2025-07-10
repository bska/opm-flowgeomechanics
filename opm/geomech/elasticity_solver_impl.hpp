//==============================================================================
//!
//! \file elasticity_upscale_impl.hpp
//!
//! \date Nov 9 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Elasticity upscale class - template implementations
//!
//==============================================================================
#ifndef OPM_ELASTICITY_SOLVER_IMPL_HPP
#define OPM_ELASTICITY_SOLVER_IMPL_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

namespace Opm::Elasticity
{

template <typename GridType>
std::vector<BoundaryGrid::Vertex>
ElasticitySolver<GridType>::extractFace(const Direction dir, const ctype coord)
{
    std::vector<BoundaryGrid::Vertex> result;
    const LeafVertexIterator itend = gv.leafGridView().template end<dim>();

    // make a mapper for codim dim entities in the leaf grid
    using LeafGridView = Dune::GridView<Dune::DefaultLeafGridViewTraits<GridType>>;
    Dune::MultipleCodimMultipleGeomTypeMapper<LeafGridView> mapper(gv.leafGridView(),
                                                                   Dune::mcmgVertexLayout());

    // iterate over vertices and find slaves
    LeafVertexIterator start = gv.leafGridView().template begin<dim>();
    for (LeafVertexIterator it = start; it != itend; ++it) {
        if (isOnPlane(dir, it->geometry().corner(0), coord)) {
            BoundaryGrid::Vertex v;
            v.i = mapper.index(*it);
            BoundaryGrid::extract(v.c, it->geometry().corner(0), log2(float(dir)));
            result.push_back(v);
        }
    }

    return result;
}

template <typename GridType>
BoundaryGrid
ElasticitySolver<GridType>::extractMasterFace(const Direction dir,
                                              const ctype coord,
                                              const SIDE side,
                                              const bool dc)
{
    static const int V1[3][4] = {{0, 2, 4, 6}, {0, 1, 4, 5}, {0, 1, 2, 3}};
    static const int V2[3][4] = {{1, 3, 5, 7}, {2, 3, 6, 7}, {4, 5, 6, 7}};
    const LeafIndexSet& set = gv.leafGridView().indexSet();

    int c = 0;
    int i = std::log2(float(dir));

    BoundaryGrid result;
    // we first group nodes into this map through the coordinate of lower left
    // vertex. we then split this up into pillars for easy processing later
    std::map<double, std::vector<BoundaryGrid::Quad>> nodeMap;
    for (LeafIterator cell = gv.leafGridView().template begin<0>();
         cell != gv.leafGridView().template end<0>();
         ++cell, ++c) {
        std::vector<BoundaryGrid::Vertex> verts;
        int idx = 0;

        if (side == LEFT) {
            idx = set.subIndex(*cell, V1[i][0], dim);
        } else if (side == RIGHT) {
            idx = set.subIndex(*cell, V2[i][0], dim);
        }

        Dune::FieldVector<double, 3> pos = gv.vertexPosition(idx);
        if (isOnPlane(dir, pos, coord)) {
            for (int j = 0; j < 4; ++j) {
                if (side == LEFT) {
                    idx = set.subIndex(*cell, V1[i][j], dim);
                }

                if (side == RIGHT) {
                    idx = set.subIndex(*cell, V2[i][j], dim);
                }

                pos = gv.vertexPosition(idx);
                if (!isOnPlane(dir, pos, coord)) {
                    continue;
                }

                BoundaryGrid::Vertex v;
                BoundaryGrid::extract(v, pos, i);

                v.i = idx;
                verts.push_back(v);
            }
        }

        if (verts.size() == 4) {
            BoundaryGrid::Quad q;
            q.v[0] = minXminY(verts);
            q.v[1] = maxXminY(verts);
            if (dc) {
                q.v[2] = minXmaxY(verts);
                q.v[3] = maxXmaxY(verts);
            } else {
                q.v[2] = maxXmaxY(verts);
                q.v[3] = minXmaxY(verts);
            }

            std::map<double, std::vector<BoundaryGrid::Quad>>::iterator it;
            for (it = nodeMap.begin(); it != nodeMap.end(); ++it) {
                if (fabs(it->first - q.v[0].c[0]) < 1.e-7) {
                    it->second.push_back(q);
                    break;
                }
            }

            if (it == nodeMap.end()) {
                nodeMap[q.v[0].c[0]].push_back(q);
            }

            result.add(q);
        }
    }

    int p = 0;
    std::map<double, std::vector<BoundaryGrid::Quad>>::const_iterator it;
    for (it = nodeMap.begin(); it != nodeMap.end(); ++it, ++p) {
        for (std::size_t ii = 0; ii < it->second.size(); ++ii) {
            result.addToColumn(p, it->second[ii]);
        }
    }

    return result;
}

template <typename GridType>
void
ElasticitySolver<GridType>::findBoundaries(double* const min, double* const max)
{
    max[0] = max[1] = max[2] = std::numeric_limits<double>::lowest();
    min[0] = min[1] = min[2] = std::numeric_limits<double>::max();

    const LeafVertexIterator itend = gv.leafGridView().template end<dim>();

    for (auto it = gv.leafGridView().template begin<dim>(); it != itend; ++it) {
        const auto& corner = it->geometry().corner(0);

        for (int i = 0; i < 3; ++i) {
            min[i] = std::min(min[i], corner[i]);
            max[i] = std::max(max[i], corner[i]);
        }
    }
}

template <typename GridType>
void
ElasticitySolver<GridType>::fixNodes(const std::vector<std::size_t>& fixed_nodes)
{
    using VertexLeafIterator = typename GridType::LeafGridView::template Codim<dim>::Iterator;
    const VertexLeafIterator itend = gv.leafGridView().template end<dim>();

    // make a mapper for codim 0 entities in the leaf grid
    using LeafGridView = typename GridType::LeafGridView;
    Dune::MultipleCodimMultipleGeomTypeMapper<LeafGridView> mapper(gv.leafGridView(),
                                                                   Dune::mcmgVertexLayout());

    NodeValue zerovec;
    zerovec = 0.0;

    // iterate over vertices
    for (VertexLeafIterator it = gv.leafGridView().template begin<dim>(); it != itend; ++it) {
        const int indexi = mapper.index(*it);
        const bool exist
            = std::find(fixed_nodes.begin(), fixed_nodes.end(), indexi) != fixed_nodes.end();
        if (exist) {
            A.updateFixedNode(indexi, std::make_pair(XYZ, zerovec));
        }
    }
}

template <typename GridType>
void
ElasticitySolver<GridType>::fixPoint(const Direction dir,
                                     const GlobalCoordinate coord,
                                     const NodeValue& value)
{
    using VertexLeafIterator = typename GridType::LeafGridView::template Codim<dim>::Iterator;
    const VertexLeafIterator itend = gv.leafGridView().template end<dim>();

    // make a mapper for codim 0 entities in the leaf grid
    using LeafGridView = Dune::GridView<Dune::DefaultLeafGridViewTraits<GridType>>;
    Dune::MultipleCodimMultipleGeomTypeMapper<LeafGridView> mapper(gv.leafGridView(),
                                                                   Dune::mcmgVertexLayout());

    // iterate over vertices
    for (VertexLeafIterator it = gv.leafGridView().template begin<dim>(); it != itend; ++it) {
        if (isOnPoint(it->geometry().corner(0), coord)) {
            const int indexi = mapper.index(*it);
            A.updateFixedNode(indexi, std::make_pair(dir, value));
        }
    }
}

template <typename GridType>
void
ElasticitySolver<GridType>::assemble(const Vector& pressure, const bool do_matrix, const bool do_vector)
{
    // int loadcase = -1;
    Dune::FieldVector<ctype, comp> eps0 = {1, 1, 1, 0, 0, 0};
    // eps0 = 0;
    Vector& b = A.getLoadVector();
    b = 0;
    A.getLoadVector() = 0;
    if (do_matrix) {
        A.getOperator() = 0;
    }

    for (int i = 0; i < 2; ++i) {
        if (color[1].size() && do_matrix) {
            std::cout << "\tprocessing " << (i == 0 ? "red " : "black ") << "elements" << std::endl;
        }

#pragma omp parallel for schedule(static)
        for (std::size_t j = 0; j < color[i].size(); ++j) {
            Dune::FieldMatrix<ctype, comp, comp> C;
            Dune::FieldMatrix<ctype, dim * bfunc, dim * bfunc> K;
            Dune::FieldMatrix<ctype, dim * bfunc, dim * bfunc>* KP = 0;
            Dune::FieldVector<ctype, dim * bfunc> ES;
            Dune::FieldVector<ctype, dim * bfunc>* EP = 0;

            if (do_matrix) {
                KP = &K;
            }

            if (do_vector) {
                EP = &ES;
            }

            for (std::size_t k = 0; k < color[i][j].size(); ++k) {
                LeafIterator it = gv.leafGridView().template begin<0>();
                for (int l = 0; l < color[i][j][k]; ++l) {
                    ++it;
                }

                const std::size_t cell_num = color[i][j][k];
                const bool valid = materials[color[i][j][k]]->getConstitutiveMatrix(C);

                if (!valid) {
                    OPM_THROW(std::runtime_error, "Not valid C matrix");
                }

                // determine geometry type of the current element and get the matching
                // reference element
                const Dune::GeometryType gt = it->type();

                Dune::FieldMatrix<ctype, dim * bfunc, dim * bfunc> Aq;
                K = 0;
                ES = 0;

                // get a quadrature rule of order two for the given geometry type
                const Dune::QuadratureRule<ctype, dim>& rule
                    = Dune::QuadratureRules<ctype, dim>::rule(gt, 2);
                for (const auto& r : rule) {
                    // compute the jacobian inverse transposed to transform the gradients
                    const auto jacInvTra = it->geometry().jacobianInverseTransposed(r.position());
                    const auto detJ = it->geometry().integrationElement(r.position());

                    if ((detJ <= 1.0e-5) && verbose) {
                        std::cout << "cell " << color[i][j][k] << " is (close to) degenerated, detJ "
                                  << detJ << std::endl;

                        double zdiff = 0.0;
                        for (int ii = 0; ii < 4; ++ii) {
                            zdiff = std::max(
                                zdiff, it->geometry().corner(ii + 4)[2] - it->geometry().corner(ii)[2]);
                        }

                        std::cout << " - Consider setting ctol larger than " << zdiff << std::endl;
                    }

                    Dune::FieldMatrix<ctype, comp, dim * bfunc> lB;
                    E.getBmatrix(lB, r.position(), jacInvTra);

                    if (do_matrix) {
                        E.getStiffnessMatrix(Aq, lB, C, detJ * r.weight());
                        K += Aq;
                    }

                    // load vector
                    if (EP) {
                        // body force
                        if (!body_force_.empty()) {
                            Dune::FieldVector<ctype, bfunc> Bvector;
                            Dune::FieldVector<ctype, dim * bfunc> lrhs;

                            // force piece wise constant over cell for now
                            const auto& force = body_force_[color[i][j][k]];
                            E.getBVector(Bvector, r.position());

                            for (std::size_t fi = 0; fi < bfunc; ++fi) {
                                for (std::size_t vd = 0; vd < dim; ++vd) {
                                    lrhs[dim * fi + vd] = Bvector[fi] * force[vd];
                                }
                            }

                            lrhs *= detJ * r.weight();
                            ES += lrhs;
                        }

                        // pressure force i.e. \div I*p
                        {
                            auto Ipressure = Dune::FMatrixHelp::mult(C, eps0);
                            Ipressure = eps0 * pressure[cell_num][0];

                            auto temp = Dune::FMatrixHelp::multTransposed(lB, Ipressure);
                            temp *= detJ * r.weight();

                            ES += temp;
                        }
                    }
                }

                A.addElement(KP, EP, it,
                             &b); // NULL is no static forse based on the itegration point??
            }
        }
    }
}

template <typename GridType>
void
ElasticitySolver<GridType>::solve()
{
    try {
        Dune::InverseOperatorResult r;
        Vector& rhs = A.getLoadVector();
        u.resize(rhs.size());
        u = 0;
        tsolver_->apply(u, rhs, r);
    } catch (Dune::ISTLError& e) {
        std::cerr << "exception thrown " << e << std::endl;
    }
}

} // namespace Opm::Elasticity

#endif
