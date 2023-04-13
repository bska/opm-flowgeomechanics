#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif

#include <dune/grid/yaspgrid.hh>

int main(int argc, char** argv){
    constexpr int dim = 3;
    using Grid = Dune::YaspGrid<dim>;
    Dune::FieldVector<double,dim> L = {4.0, 4,0, 4.0};
    std::array<int,dim> s = {1,1,1};
    Grid grid( L , s);
    //auto gv = grid.leafGridView();
    //auto set = gv.indexSet();
    //asmhandler_impl.hpp
    //elasticity_upscale_impl.hpp
    for(const auto& cell:elements(gv)){
        //auto gm = cell->geometry();
        //int index2 = set.subIndex(*cell,j,dim);
        for (const auto& is: Dune::intersections(gv,cell)){
            //std::cout << is << std::endl;
        }
        //for (const auto& ed: Dune::vertices(gv,cell)){
        //}
        
    }
}
