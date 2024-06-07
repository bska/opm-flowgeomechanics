#include <cstddef>
#include <dune/common/exceptions.hh>
#include <dune/istl/solvercategory.hh>
#include <iostream>
#include <fstream>
#include <vector>
#include <assert.h>
#include <stdexcept>

#include <dune/istl/matrixmarket.hh>
#include <dune/istl/preconditioners.hh>
//#include <dune/istl/paamg/amg.hh>

@@ there must be a more correct way to ensure UMFpack is included here
#define HAVE_SUITESPARSE_UMFPACK 1
#include <dune/istl/umfpack.hh>

//#include "opm/geomech/Fracture.hpp"
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/common/mcmgmapper.hh> // mapper class
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/common/filledarray.hh> // needed for printSparseMatrix??
#include <dune/istl/io.hh> // needed for printSparseMatrix??


#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid.hh>

#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/setupPropertyTree.hpp>
#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>

#include "opm/geomech/DiscreteDisplacement.hpp"


namespace {
template<class M, class X, class Y> class ReducedMatrixAdapter; // forward declaration
  
using Grid = Dune::FoamGrid<2, 3>;
using Vector = Dune::BlockVector<double>;
using VectorHP = FieldVector<Vector, 2>;
using HTrans = std::tuple<size_t,size_t, double, double>;
  
using FullMatrix = Dune::DynamicMatrix<double>;
using SparseMatrix = Dune::BCRSMatrix<double>;
using EquationSystem = std::tuple<std::shared_ptr<FullMatrix>,    // aperture matrix
                                  std::shared_ptr<SparseMatrix>,  // pressure matrix
                                  std::vector<HTrans>,            // inter-cell transmissibility factors
                                  std::vector<double>>;           // leakoff factors

using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<Grid::LeafGridView>;
using RMAdapter = ReducedMatrixAdapter<SparseMatrix, Vector, Vector>;
// ----------------------------------------------------------------------------
std::unique_ptr<const Grid> readgrid(const char* const name)
// ----------------------------------------------------------------------------  
{
  std::unique_ptr<const Grid> result;
  try {
    result = Dune::GmshReader<Grid>::read(name); // unique_ptr
  } catch (const Dune::Exception& e) {
    std::cout << "Unable to read file: " << name << std::endl;
  }
  return result;
}

// ----------------------------------------------------------------------------
template<int N>
std::array<size_t, N> n_closest(const Grid& grid)
// ----------------------------------------------------------------------------
{
  using Elem = std::tuple<size_t, double>;
  std::vector<Elem> distances;
  const ElementMapper mapper(grid.leafGridView(), Dune::mcmgElementLayout());
  for (const auto& element : Dune::elements(grid.leafGridView())) {
    const int eIdx = mapper.index(element);
    const auto center = element.geometry().center();
    distances.push_back({eIdx, center.two_norm()});
  }

  std::sort(distances.begin(), distances.end(),
            [](Elem& a, Elem& b) {return std::get<1>(a) < std::get<1>(b);});
  
  std::array<size_t, N> result;
  for (int i = 0; i != N; ++i)
    result[i] = std::get<0>(distances[i]);
  
  return result;
}

// ----------------------------------------------------------------------------  
std::shared_ptr<FullMatrix> computeApertureMatrix(const Grid& G,
                                                  const double young,
                                                  const double poisson)
// ----------------------------------------------------------------------------  
{
  const int nc = G.leafGridView().size(0);
  std::shared_ptr<FullMatrix> result(new FullMatrix(nc, nc));
  ddm::assembleMatrix(*result, young, poisson, G);

  *result *= -1;
  return result;
}

// ----------------------------------------------------------------------------
  std::tuple<std::vector<HTrans>, std::vector<double>>
  computeHTrans(const Grid& grid, const double leakoff_fac)
// ----------------------------------------------------------------------------
{
  std::vector<HTrans> htransvec;
  std::vector<double> leakvec(grid.leafGridView().size(0), 0.0);

  const ElementMapper mapper(grid.leafGridView(), Dune::mcmgElementLayout());
  for (const auto& element : Dune::elements(grid.leafGridView())) {
    const int eIdx = mapper.index(element);
    const auto geom = element.geometry();
    for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
      if (is.boundary())
        continue;
      const int nIdx = mapper.index(is.outside());
      const auto igeom = is.geometry();
      if( eIdx < nIdx) {
        const auto iscenter = igeom.center();
        const auto ncenter = is.outside().geometry().center() - iscenter;
        const auto ecenter = geom.center() - iscenter;
        const double area = igeom.volume();
        const double h1 = ecenter.two_norm() / area;
        const double h2 = ncenter.two_norm() / area;

        htransvec.push_back({nIdx, eIdx, h1, h2});
      }
    }
    leakvec[eIdx] = geom.volume() * leakoff_fac;
  }
  return {htransvec, leakvec};
}

  
// ----------------------------------------------------------------------------
std::tuple<std::shared_ptr<SparseMatrix>, std::vector<HTrans>, std::vector<double>>
computePressureMatrix(const Grid& G, const double leakoff_fac)
// ----------------------------------------------------------------------------
{
  // initialize the sparsity pattern of a pressure matrix, with all entries set to zero.

  std::tuple<std::vector<HTrans>,
             std::vector<double>> hvec = computeHTrans(G, leakoff_fac);

  const size_t nc = G.leafGridView().size(0);
  std::shared_ptr<SparseMatrix> mat(new SparseMatrix());
  mat->setBuildMode(SparseMatrix::implicit);
  mat->setImplicitBuildModeParameters(3 * 6 - 3 - 2, 0.4);
  mat->setSize(nc, nc);

  for (const auto& elem : std::get<0>(hvec)) {
    const size_t ix(std::get<0>(elem));
    const size_t jx(std::get<1>(elem));
    mat->entry(ix, jx) = 0.0;
    mat->entry(jx, ix) = 0.0;
    mat->entry(jx, jx) = 0.0;
    mat->entry(ix, ix) = 0.0;
  }
  mat->compress();
  
  return std::tuple_cat(std::make_tuple(mat), hvec); 
}
  
  
// ----------------------------------------------------------------------------
EquationSystem computeEquationSystem(const Grid& G,
                                     const double young, const double poisson,
                                     const double leakoff_fac)
// ----------------------------------------------------------------------------
{
  return std::tuple_cat(std::make_tuple(computeApertureMatrix(G, young, poisson)),
                        computePressureMatrix(G, leakoff_fac));
}

// ----------------------------------------------------------------------------
template<class M, class X, class Y>
class ReducedMatrixAdapter : public Dune::LinearOperator<X, Y> 
// ----------------------------------------------------------------------------
{
public:
  ReducedMatrixAdapter(const M& mat, const std::vector<size_t> elim) :
    mat_(mat), elim_(sorted_(elim)), keep_(nonelim_(elim_, mat.M())),
    tmpX_(mat.M()), tmpY_(mat.M()) {}
  virtual ~ReducedMatrixAdapter() {}

  typedef M matrix_type;
  typedef X domain_type;
  typedef Y range_type;
  typedef typename X::field_type field_type;

  virtual void apply(const X& x, Y& y) const
  {
    expand(x, tmpX_); // expand x into tmpX_
    mat_.mv(tmpX_, tmpY_);
    contract(tmpY_, y); //contract tmpY_ into y
  }

  virtual void applyscaleadd(field_type alpha, const X& x, Y& y) const
  {
    expand(x, tmpX_); // expand x into tmpX_
    expand(y, tmpY_); // expand y into tmpY_
    mat_.usmv(alpha, tmpX_, tmpY_);
    contract(tmpY_, y); //contract tmpY_ into y    
  }
  
  virtual Dune::SolverCategory::Category category() const
  {
    return Dune::SolverCategory::sequential;
  }

  void expand(const X& source, X& target) const {
    target.resize(fullSize());
    target = 0;
    for (size_t i = 0; i != reducedSize(); ++i)
      target[keep_[i]] = source[i];
  }

  void contract(const Y& source, Y& target) const {
    target.resize(reducedSize());
    for (size_t i = 0; i != reducedSize(); ++i)
      target[i] = source[keep_[i]];
  }

  size_t fullSize() const {return mat_.M();}
  size_t reducedSize() const {return keep_.size();}

private:

  std::vector<size_t> sorted_(const std::vector<size_t>& vec)
  {
    std::vector<size_t> result(vec);
    std::sort(result.begin(), result.end());
    return result;
  }

  std::vector<size_t> nonelim_(const std::vector<size_t>& elim, size_t nc)
  {
    // return vector with all incides that are _not_ found in elim
    std::vector<size_t> all_ixs(nc);
    std::iota(all_ixs.begin(), all_ixs.end(), 0);
    std::vector<size_t> nonelim;
    std::set_difference(all_ixs.begin(), all_ixs.end(), elim.begin(), elim.end(),
                        std::back_inserter(nonelim));
    return nonelim;
  }
      
  const M& mat_;
  const std::vector<size_t> elim_;
  const std::vector<size_t> keep_;
  mutable X tmpX_;
  mutable Y tmpY_;
};

// ----------------------------------------------------------------------------  
template<typename T> inline T hmean(const T a, const T b) {return T(1) / ( T(1)/a + T(1)/b);};
// ----------------------------------------------------------------------------
  
// ----------------------------------------------------------------------------
void updateTrans(SparseMatrix& mat, const std::vector<HTrans>& htransvec,
                 const Vector& h, const std::vector<double>& leakvec)
// ----------------------------------------------------------------------------
{
  mat = 0; // reset matrix before starting to fill in values

  for (const auto& e : htransvec) {
    const size_t i = std::get<0>(e);
    const size_t j = std::get<1>(e);
    assert(i != j);
    const double t1 = std::get<2>(e);
    const double t2 = std::get<3>(e);
    const double h1 = h[i];
    const double h2 = h[j];

    const double trans1 = h1 * h1 * t1 / 12.0;
    const double trans2 = h2 * h2 * t2 / 12.0;
    
    const double T = hmean(trans1, trans2);

    // update diagonal
    mat[i][i] += T;
    mat[j][j] += T;

    // update off-diagonal terms
    mat[i][j] -= T;
    mat[j][i] -= T;
    
  }

  // add-in leakage term on diagonal
  for (size_t ix = 0; ix != mat.N(); ++ix)
    mat[ix][ix] += leakvec[ix];
  
}

// ----------------------------------------------------------------------------
// Note: 'elim' and 'keep' should both be sorted, complementary, and together
// contain all the indices from 0 to M.N() and/or M.M().
// The matrix 'M' will be reduced, and any columns extracted will be collected and
// retuerned. (Rows will just be removed).
SparseMatrix reduceMatrix(SparseMatrix& M, bool rows, bool cols,
                          const std::vector<size_t> elim,
                          const std::vector<size_t> keep)
// ----------------------------------------------------------------------------
{
  if (cols) assert(elim.size() + keep.size() == M.N());
  if (rows) assert(elim.size() + keep.size() == M.M());

  SparseMatrix Mreduced, cols;

  Mreduced.setBuildMode(SparseMatrix::implicit);
  Mreduced.setImplicitBuildModeParameters(3 * 6 - 3 - 2, 0.4);
  Mreduced.setSize(rows ? keep.size() : M.M(),
                   cols ? keep.size() : M.N());

  cols.setBuildMode(SparseMatrix::implicit);
  cols.setImplicitBuildModeParameters(3*6-3-2, 0.4);
  cols.setSize(rows ? elim.size() : M.M(),
               cols ? elim.size() : 0);

  vector<size_t> mapping(keep.size() + elim.size(), 0);
  for (size_t i = 0; i != keep.size(); ++i)  mapping[keep[i]] = i;
  for (size_t i = 0; i != elim.size(); ++i)  mapping[elim[i]] = -(i+1);
  
  for (auto rowIt = M.begin(); rowIt != M.end(); ++ rowIt) {
    const size_t i = rowIt->index();
    if (rows && mapping[i] < 0)
      continue; // this row should be eliminated, nothing further to do
    
    for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt) {
      const size_t j = colIt->index();
      const bool elim_col = cols && mapping[j] < 0;
      
      if (cols && mapping[j] < 0)
        // this column should be eliminated
        cols.entry(rows ? mapping[i] : i, -mapping[i] - 1) = M[i][j];
      else 
        // this column should be kept
        Mreduced.entry(i, mapping[j]) = M[i][j];
    }
  }
  
  Mreduced.compress();
  cols.compress();
  M.swap(Mreduced);
  return cols;
}

  
// ----------------------------------------------------------------------------
class FlowSystemMatrices
// ----------------------------------------------------------------------------
{
public:
  FlowSystemMatrices(const SparseMatrix& M,
                     const std::vector<HTrans>& htransvec,
                     const std::vector<double>& leakvec) :
    M_(M), C_(M), rhs_(), htransvec_(htransvec), eliminated_indices_(),
    kept_indices_(M.N()), remapping_(), fixed_bhp_(0), fixed_rate_(0), ratecells_()
  { std::iota(kept_indices_.begin(), kept_indices_.end(), 0);
    std::iota(remapping_.begin(), remapping_.end(), 0); }

  void setBHPCells(const std::vector<size_t>& bhp_cells, const double bhp)  {
    if (!eliminated_indices_.empty())
      throw std::runtime_error("System has already been reduced.  Cannot re-set BHP cells"); 

    eliminated_indices_ = sorted_(bhp_cells);
    kept_indices_ = noneliminated_(eliminated_indices_, M_.M());
    remapping_ = compute_remapping_(eliminated_indices_, kept_indices_, M_.M());

    reduceMatrix(M, true, true, eliminated_indices_, kept_indices_);
    reduceMatrix(C, true, false, eliminated_indices_, kept_indices_);

    rhs_.resize(M.N());
    fixed_bhp_ = bhp;
  }
  void setRateCells(const std::vector<size_t>& rate_cells, const double rate) {
    fixed_rate_ = rate;
    ratecells_ = rate_cells;
  }
  size_t pnum() const { return M_.N(); } // number of (remaining) pressure values
  size_t hnum() const { return C_.N(); } // number of aperture values
  void update(const VectorHP& hp) {
    update_M_and_rhs_(hp);
    update_C_(hp);
  } 
  const SparseMatrix& M() const {return M_;}
  const SparseMatrix& C() const {return C_;}
  Vector rhs() const { return rhs_;} // should only be called after update
  Vector reduceP(const Vector& p) {
    Vector result(pnum()); result = 0;
    for (size_t i = 0; i != pnum(); ++i)
      result[i] = p[kept_indices_[i]];
    return result;
  }
private:

  static std::vector<size_t> sorted_(const std::vector<size_t>& vec)  {
    std::vector<size_t> result(vec);
    std::sort(result.begin(), result.end());
    return result;
  }
  static std::vector<size_t> noneliminated_(const std::vector<size_t>& elim, size_t nc)  {
    // return vector with all incides that are _not_ found in elim
    std::vector<size_t> all_ixs(nc);
    std::iota(all_ixs.begin(), all_ixs.end(), 0);
    std::vector<size_t> nonelim;
    std::set_difference(all_ixs.begin(), all_ixs.end(), elim.begin(), elim.end(),
                        std::back_inserter(nonelim));
    return nonelim;
  }
  static std::vector<int> compute_remapping_(const std::vector<size_t>& eliminated,
                                                const std::vector<size_t>& kept,
                                                const size_t nc) {
    std::vector<int> result(nc);
    for (size_t i = 0; i != eliminated.size(); ++i) result[eliminated[i]] = -i;
    for (size_t i = 0; i != kept.size(); ++i) result[kept[i]] = i;
    return result;
  }
  void update_M_and_rhs_(const VectorHP& hp) {
    // update M_ and rhs_reduct_
    M_ = 0; rhs_ = 0;
    for (const auto& e : htransvec_) {
      const size_t i = std::get<0>(e);
      const size_t j = std::get<1>(e);
      assert(i != j);
      const double t1 = std::get<2>(e);
      const double t2 = std::get<3>(e);
      const double h1 = hp[0][i];
      const double h2 = hp[0][j];
      
      const double trans1 = h1 * h1 * t1 / 12.0;
      const double trans2 = h2 * h2 * t2 / 12.0;
      
      const double T = hmean(trans1, trans2);

      // remapped indices
      const int ir = remapping_[i];
      const int jr = remapping_[j];

      // update diagonal
      (ir > 0) && M_[ir][ir] += T; 
      (jr > 0) && M_[jr][jr] += T; 

      // update off-diagonal terms, and the reduction terms of rhs_reduct_
      if (ir > 0) if (jr > 0) M_[ir, jr] += T; else rhs_[ir] -= T * bhp_;
      if (jr > 0) if (ir > 0) M_[jr, ir] += T; else rhs_[jr] -= T * bhp_;
    }

    // // update the rest of the contribution of M on rhs
    // Vector p_reduced(M.N()); p_reduced = 0;
    // for (size_t i = 0; i != p_reduced.size(); ++i)
    //   p_reduced[i] = hp[1][kept_indices_[i]];
    
    // M_.mmv(p_reduced, rhs_); // rhs_ = rhs_ - M * p_reduced
    // // NB: note that update_C_ must also be called before rhs_ is complete
  }
  void update_C_(const VectorHP& hp) {
    C_ = 0;
    // intermediary matrices for computing the partial derivatives
    SparseMatrix Q(C_), R(C_), dQ(C_), dR(C_);
    Q = 0; R = 0; dQ = 0; dR = 0;

    for (const auto& e : htransvec_) {
      const size_t i = std::get<0>(e);
      const size_t j = std::get<1>(e); assert(i != j);

      const int ir = remapping_[i];
      if (ir < 0) continue; // this equation was eliminated
      
      const double t1 = std::get<2>(e);
      const double t2 = std::get<3>(e);
      const double h1 = hp[0][i];
      const double h2 = hp[0][j];
      const double p1 = hp[1][i];
      const double p2 = hp[1][j];

      const double q   = h1 * h1 * h2 * h2 * t1 * t2;
      const double d1q = 2 * h1 * h2 * h2 * t1 * t2;
      const double d2q = h1 * h1 * 2 * h2 * t1 * t2;

      const double r   = 12 * (h1 * h1 * t1 + h2 * h2 * t2);
      const double d1r = 24 * h1 * t1;
      const double d2r = 24 * h2 * t2;
      
      C_[i][j] = (d1q * r - q * d1r) / (r * r) * (p1 - p2);
      C_[j][i] = (d2q * r - q * d2r) / (r * r) * (p2 - p1);
      
    }

    // // add contribution to rhs
    // C_.mmv(hp[0], rhs_); // rhs_ = rhs_ - C_ * hp[0]
  }

  SparseMatrix M_;
  SparseMatrix C_;
  Vector rhs_;
  const std::vector<HTrans>& htransvec_;
  const std::vector<double>& leakvec_;
  std::vector<size_t> eliminated_indices_;
  std::vector<size_t> kept_indices_;
  std::vector<int> remapping_ ;
  double fixed_bhp_;
  double fixed_rate_;
  std::vector<size_t> ratecells_;
};

// ----------------------------------------------------------------------------
SparseMatrix makeIdentity(size_t num)
// ----------------------------------------------------------------------------
{
  //build a sparse matrix to represent the identity matrix of a given size
  SparseMatrix M;
  M.setBuildMode(SparseMatrix::implicit);
  M.setImplicitBuildModeParameters(1, 0.1);
  M.setSize(num, num);
  for (size_t i = 0; i != num; ++i)
    M.entry(i,i) = 1;
  M.compress();
  return M;
}
  
// ----------------------------------------------------------------------------
VectorHP& computeIncrementAndResidual(const FullMatrix& A,
                                      const FlowSystemMatrices& SFM, // will be updated
                                      const VectorHP& hp,
                                      VectorHP& dhp, // increment (to be computed)
                                      VectorHP& residual) // residual (to be computed)
// ----------------------------------------------------------------------------
{
  using SystemMatrix =
    Dune::MultiTypeBlockMatrix<Dune::MultiTypeBlockVector<FullMatrix,   SparseMatrix>,
                               Dune::MultiTypeBlockVector<SparseMatrix, SparseMatrix>>;

  SFM.update(hp);

  SystemMatrix S { { A      , makeIdentity(SFM.pnum()) },
                   { SFM.C(), SFM.M()                  } };

  // we use the 'residual' vector also to represent the right-hand-side
  residual[0] = 0;
  residual[1] = SFM.rhs();

  Vector hp_reduced = hp;
  hp_reduced[1] = SFM.reduceP(hp[1]);

  S.mmv(hp_reduced, residual); // residual = residual - S * hp_reduced

  // solve system equations
  Dune::UMFPack psolver(S);
  Dune::InverseOperatorResult r;  
  psolver.apply(dhp, residual, r);

  // compute residual
  VectorHP tmp(hp);
  tmp += dhp;
  S.mmv(tmp, residual);
  
  return residual;
}

// ----------------------------------------------------------------------------
int solveCoupledFull(Vector& p, Vector&h, const EquationSystem& eqsys,
                       const std::vector<size_t> bhpcells, const double bhp,
                       const std::vector<size_t> ratecells, const double rate, 
                       const double convergence_tol=1e-4, const int max_nonlin_iter=40)
// ----------------------------------------------------------------------------
{
  // We consider the nonlinear system:
  // | A          -I |  | h |   | 0 |   (mechanics equation)
  // |               |  |   | = |   |
  // | C(h, p)  M(h) |  | p |   | q |   ( flow equation)

  const FullMatrix&          A         = *std::get<0>(eqsys);
  const SparseMatrix&        M         = *std::get<1>(eqsys);
  const std::vector<HTrans>& htransvec =  std::get<2>(eqsys);
  const std::vector<double>& leakvec   =  std::get<3>(eqsys);

  VectorHP hp(h, p); // solution vector, grouping aperture 'h' and pressure 'p'

  FlowSystemMatrices FSM(hp, M, htransvec, leakvec);
  FSM.setBHPCells(bhpcells, bhp);
  FSM.setRateCells(ratecells, rate);

  auto converged = [&](const VectorHP& res) { return res.infinity_norm() < convergence_tol; };
  
  VectorHP residual, dhp;
  int iter = 0;

  // nonlinear solve loop
  while (!converged( computeIncrementAndResidual(A, FSM, hp, dhp, residual) ) &&
         ++iter < max_nonlin_iter) 
    hp += dhp;

  if (iter == max_nonlin_iter) {
      std::cout << "System did not converge in max allowed number of iterations." << std::endl;
      return -1;
  }

  // system converged, unpacking results  
  h = hp[0];
  p = hp[1];

  return 0;
};    
  

// ----------------------------------------------------------------------------
void solveCoupledSplit(Vector& p, Vector& h, const EquationSystem& eqsys,
                       const std::vector<size_t> bhpcells, const double bhp,
                       const std::vector<size_t> ratecells, const double rate, 
                       const double convergence_tol=1e-4, const int max_nonlin_iter=400)
// ----------------------------------------------------------------------------
{
  // NB: there should be no overlap between the cells in 'bhpcells' and 'ratecells'
  // reduce pressure system

  //Dune::MatrixAdapter<FullMatrix, Vector, Vector> MA_h(*std::get<0>(eqsys));
  const FullMatrix&          M_h       = *std::get<0>(eqsys);
  SparseMatrix&              M_p       = *std::get<1>(eqsys);
  const std::vector<HTrans>& htransvec =  std::get<2>(eqsys);
  const std::vector<double>& leakvec   =  std::get<3>(eqsys);

  const RMAdapter MA_p(M_p, bhpcells); // allows reducing the system without creating new matrices

  // define convergence criterion
  auto converged = [&](const Vector& v1, const Vector& v2) {
    auto diff = v1; diff -= v2; // diff = v1 - v2
    const double delta = diff.infinity_norm();
    const double denom = std::max(v1.infinity_norm(), v2.infinity_norm());
    std::cout << "Delta: " << delta << " Denom: " << denom;
    std::cout << " Delta/denom: " << delta/denom << std::endl;
    return delta/denom < convergence_tol;
  };

  // helper function to set vector values
  auto setvals = [](Vector& v, const std::vector<size_t>& ixs, const double val) {
    for (auto i : ixs) v[i] = val;
  };
    
  // initialize pressure
  p = 0;
  if (bhpcells.size() > 0)
    setvals(p, bhpcells, bhp);
  else
    p = 1e6; // we need something to get started without h being 0.
      
  // solve for aperture
  M_h.solve(h, p); // solve for aperture (h) given pressure (p)

  Vector htmp(h.N()); htmp = 0;
  Vector rhs_full_rate(p.N()), rhs_full(p.N());
  Vector rhs;
  Dune::InverseOperatorResult iores;  
  rhs_full_rate = 0; setvals(rhs_full_rate, ratecells, rate);

  int i;
  for (i = 0; i != max_nonlin_iter; ++i) {
  
    // update pressure matrix entries
    updateTrans(M_p, htransvec, h, leakvec);

    // determine right-hand side modifications from eliminated degrees of freedom
    
    M_p.mv(p, rhs_full); // only imposed values of p should be nonzero here.
    rhs_full *= -1;
    rhs_full += rhs_full_rate;
    
    MA_p.contract(rhs_full, rhs); // remove entries corresponding to eliminated equations
  
    // solve for pressure
    Dune::Richardson<Vector, Vector> precond(1); // "no" preconditioner 
    //Dune::SeqJac<SparseMatrix, Vector, Vector> precond(M_p, 3, 0.3);
    // using Smoother = Dune::SeqSSOR<SparseMatrix, Vector, Vector>;
    // Dune::Amg::SmootherTraits<Smoother>::Arguments smootherArgs;
    // smootherArgs.iterations = 3;
    // smootherArgs.relaxationFactor = 1;

    // auto criterion =
    //   Dune::Amg::CoarsenCriterion<Dune::Amg::UnSymmetricCriterion<SparseMatrix,
    //                                                               Dune::Amg::FirstDiagonal>>(15, 50);
    // criterion.setDefaultValuesIsotropic(2);
    // criterion.setAlpha(.67);
    // criterion.setBeta(1.0e-4);
    // criterion.setGamma(1);
    // criterion.setDebugLevel(2);
    // using AMG = Dune::Amg::AMG<RMAdapter, Vector, Smoother>;
    // AMG precond(MA_p, criterion, smootherArgs);

    auto psolver = Dune::CGSolver<Vector>(MA_p,
                                          precond, 
                                          1e-7, // desired residual reduction factor
                                          100, // max number of iterations
                                          1); // verbose
                                          
    Vector p_reduced(MA_p.reducedSize());
    psolver.apply(p_reduced, rhs, iores);
    MA_p.expand(p_reduced, p);
    setvals(p, bhpcells, bhp);

  // solve for aperture again
    M_h.solve(htmp, p); // solve for aperture (h) given pressure (p)

    if (converged(h, htmp))
      break;

    //h = htmp;
    h *= 0.0;
    htmp *= 1.0;
    h += htmp;
    
    
  }
    
  if (i == max_nonlin_iter)
    std::cout << "Warning, did not converge in max number of nonlinear iterations." << std::endl;
  else
    std::cout << "Converged in: " << i << " iterations." << std::endl;
}
  
}; // end anonymous namespace

// ============================================================================
int main(int varnum, char** vararg)
// ============================================================================
{
  // check argument list
  if (varnum < 2) {
    std::cout << "Grid file name must be given as argument." << std::endl;
    return -1;
  }
  
  // read grid
  const auto grid = readgrid(vararg[1]);
  if (!grid)
    return -1;

  // identify well cells
  const auto wellcells = n_closest<2>(*grid);

  std::cout << "Identified wellcells: ";
  std::copy(wellcells.begin(), wellcells.end(), std::ostream_iterator<size_t>(std::cout, " "));
  std::cout << std::endl;
  
  // compute equation system
  const double young = 1e9; // Young's modulus
  const double poisson = 0.25; // Poisson's ratio
  const double leakoff_fac = 1e-6; //1e-13; // a bit heuristic; conceptually rock perm divided by distance  
  const auto eqsys = computeEquationSystem(*grid, young, poisson, leakoff_fac);
  
  // solve fixed pressure system
  const size_t nc = grid->leafGridView().size(0);
  Vector pressure(nc), aperture(nc);
  pressure = 0;
  const double bhp = 1e7; // in Pascal
  const double rate = -0.1;

  // solve fixed pressure system  
  // solveCoupledSplit(pressure, aperture, eqsys,
  //              std::vector<size_t>(wellcells.begin(), wellcells.end()), bhp,
  //              std::vector<size_t>(), rate);

  // solve fixed rate system
  solveCoupledSplit(pressure, aperture, eqsys,
               std::vector<size_t>(), bhp,
               std::vector<size_t>(wellcells.begin(), wellcells.end()), rate);


  // write output
  Dune::VTKWriter<Grid::LeafGridView> vtkwriter(grid->leafGridView(), Dune::VTK::nonconforming);
  vtkwriter.addCellData(aperture, "aperture");
  vtkwriter.addCellData(pressure, "pressure");
  vtkwriter.write("output");
  
}; // end main
