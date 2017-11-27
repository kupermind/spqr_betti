unset LD_LIBRARY_PATH
module load gcc/4.9.0
g++ spqr_betti.cpp -I/home/ovcharen/spack_install/viz/linux-rhel6-x86_64/gcc-4.9.0/suite-sparse-4.5.6-d7jhkrppfturpqitwamax2bx5yjtcbcl/include -L/home/ovcharen/spack_install/viz/linux-rhel6-x86_64/gcc-4.9.0/suite-sparse-4.5.6-d7jhkrppfturpqitwamax2bx5yjtcbcl/lib -lspqr -lcholmod -fopenmp -o spqr_betti
