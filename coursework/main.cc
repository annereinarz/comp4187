#include <vector>
#include <iostream>
#include <cmath>

class Grid {
    public:
        Grid(int levels_) : levels(levels_) {
            n.resize(levels);
        }

        inline int multiindex(int i, int j, int level) const {
            return i*n[level] + j;
        }

        inline void setn(int n_, int level) {
            n[level] = n_;
        }

        inline int getn(int level) {
            return n[level];
        }

        inline int getlevels() {
            return levels;
        }

    private:
        int levels;
        std::vector<int> n;
};

// Gauss - Seidel relaxation
void gaussseidel(Grid& grid, int level, std::vector<double>& sol, std::vector<double>& rhs) {
    // Iterate over all nodes that are not boundary nodes
    for(int i = 1; i < grid.getn(level)-1; i++) {
        for(int j = 1; j < grid.getn(level)-1; j++) {
            double invh2 = pow(grid.getn(level)-1,2);
            double denom = 1.0/(4.0 * invh2);
            sol[grid.multiindex(i, j, level)] = (rhs[grid.multiindex(i, j, level)] +
                    invh2 * (sol[grid.multiindex(i+1, j, level)] 
                           + sol[grid.multiindex(i-1, j, level)]
                           + sol[grid.multiindex(i, j+1, level)]
                           + sol[grid.multiindex(i, j-1, level)] )) * denom;
        }
    }
}

void restrict(Grid& grid, int level, std::vector<double>& fine, std::vector<double>& coarse) {
    // loop over coarse grid points
    int n = grid.getn(level);
    for(int i = 1; i < n - 1; i++) {
        int fi = 2 * i;
        for( int j =1; j < n -1; j++) {
            int fj = 2 * j;
            coarse[grid.multiindex(i, j, level)] =  0.25  *
                ( fine[grid.multiindex(fi,fj,level)] + fine[grid.multiindex(fi -1 , fj, level)] 
                  + fine[grid.multiindex(fi , fj -1, level)] + fine[grid.multiindex(fi -1 , fj -1,level )] );
        }
    }
}

double residual(Grid& grid, int level, std::vector<double>& sol, std::vector<double>& rhs) {
    double res = 0.0;
    double invh2 = pow(grid.getn(level)-1,2);
    auto apply_matrix = [=](int i, int j) {
        return invh2 * ( sol[grid.multiindex(i+1,j, level)]
                + sol[grid.multiindex(i,j+1, level)]
                + sol[grid.multiindex(i,j-1, level)]
                + sol[grid.multiindex(i-1,j, level)]
                - 4.0 * sol[grid.multiindex(i, j, level)]);
    };
    double rf;
    for(int i = 1; i < grid.getn(level)-1; i++) {
        for(int j = 1; j < grid.getn(level)-1; j++) {
            rf = rhs[grid.multiindex(i,j,level)] + apply_matrix(i, j);
            res += rf * rf;
        }
    }
    return sqrt(res)/invh2;
}

void restrict_residual(Grid& grid, int level, std::vector<std::vector<double>>& sol, std::vector<std::vector<double>>& rhs) {
    double invh2 = pow(grid.getn(level)-1,2);
    auto get_residual = [=](int i, int j) {
        return rhs[level][grid.multiindex(i,j,level)] + invh2 * ( sol[level][grid.multiindex(i+1,j, level)]
                + sol[level][grid.multiindex(i,j+1, level)]
                + sol[level][grid.multiindex(i,j-1, level)]
                + sol[level][grid.multiindex(i-1,j, level)]
                - 4.0 * sol[level][grid.multiindex(i, j, level)]);
    };

    // loop over coarse grid points
    for(int i=1; i < grid.getn(level+1)-1; i++) {
        int fi = 2 * i;
        for(int j =1; j < grid.getn(level+1)-1; j++) {
            int fj = 2 * j;
            rhs[level+1][grid.multiindex(i , j, level+1)] = 0.25  *
                (get_residual(fi , fj) + get_residual(fi -1 , fj) +
                  get_residual(fi , fj -1) + get_residual(fi -1 , fj -1));
        }
    }
}

void interpolate(Grid& grid, int level, std::vector<double>& fine, std::vector<double>& coarse) {
    double v ;
    // loop over coarse grid points
    for(int i=1; i < grid.getn(level)-1; i++) {
        int fi=2*i;
        for(int j =1; j < grid.getn(level)-1; j++) {
            int fj = 2*j;
            v = coarse[grid.multiindex(i, j, level)];
            fine[grid.multiindex( fi , fj , level-1)] += v ;
            fine[grid.multiindex( fi -1 , fj, level-1)] += v ;
            fine[grid.multiindex( fi , fj -1, level-1)] += v ;
            fine[grid.multiindex( fi -1 , fj -1, level-1)] += v ;
        }
    }
}

void VCycle (Grid& grid, int level, std::vector<std::vector<double>>& sol,
        std::vector<std::vector<double>>& rhs, int npre, int npost, int ncoarse) {
    // on coarsest level solve using smoother
    if(level == grid.getlevels()-1) {
        for(int i = 0; i < ncoarse; i++) {
            gaussseidel(grid, level, sol[level], rhs[level]);
        }
    }
    // or recursively perform V-cycle
    else {
        // pre-smoothing steps
        for(int i = 0; i < npre; i++) {
            gaussseidel(grid, level, sol[level], rhs[level]);
        }
        // compute and restrict the residual
        restrict_residual(grid, level, sol, rhs);
        // initialize the coarse solution to zero and recursively call V-cycle
        std::fill((sol[level+1]).begin(), (sol[level+1]).end(), 0);
        VCycle(grid, level+1, sol, rhs, npre , npost , ncoarse);
        // interpolate error and correct fine solution
        interpolate(grid, level+1, sol[level], sol[level+1]);
        // post-smoothing steps
        for(int i = 0; i < npost; i++) {
            gaussseidel(grid, level, sol[level], rhs[level]);
        }
    }
}

int main (int argc , char ** argv) {
    // Read input
    if(argc != 7) {
        std::cout << "Required arguments: " << std::endl
            << "n: 2^n + 1 is number of points in x and y direction," <<std::endl
            << "levels: number of multigrid levels to user" << std::endl
            << "iters: number of multigrid cycles to perform " << std::endl 
            << "npre: number of pre-smoother steps" << std::endl
            << "npost: number of post-smoother steps" << std::endl
            << "ncoarse: number of smoother steps to solve coarse system " << std::endl ;
        exit(0);
    }
    int n = pow(2,atoi(argv[1])), maxlevel = atoi(argv[2]), iters = atoi(argv[3]);
    int npre = atoi(argv[4]), npost = atoi(argv[5]), ncoarse = atoi(argv[6]);

    // Initialise grid
    Grid grid(maxlevel);

    // Allocate memory for solution and right hand side vectors on each level
    std::vector<std::vector<double>> sol(maxlevel);
    std::vector<std::vector<double>> rhs(maxlevel);
    for (int i = 0; i < maxlevel; i++) {
        sol[i].resize(n*n); rhs[i].resize(n*n);
        grid.setn(n,i);
        std::cout << "Size " << grid.getn(i) << std::endl;
        n = n/2;
    }

    // Initialise right hand side on finest level and restrict to coarser levels
    auto f = [](double x, double y) {
        return 2*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
    };

    double h = 1.0/(grid.getn(0)-1);
    for(int i = 1; i < grid.getn(0)-1; i++) {
        for(int j = 1; j < grid.getn(0)-1; j++) {
            rhs[0][grid.multiindex(i,j,0)] = f(i*h,j*h);
        }
    }
    for(int l = 1; l < maxlevel; l++)
        restrict(grid, l, rhs[l-1] , rhs[l]);

    // Call V-Cycle
    for(int i = 0; i < iters; i++) {
        VCycle(grid,0,sol,rhs, npre, npost, ncoarse);
        std::cout << "Residual after " << i << " iterations of VCycle is: " 
            << residual(grid,0,sol[0],rhs[0]) << std::endl;
    }
    return 0;
}

