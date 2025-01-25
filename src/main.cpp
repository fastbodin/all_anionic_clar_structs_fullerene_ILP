#include "include.h"

int main(int argc, char *argv[]) {
  // want to solve for the p-anionic Clar number
  int p = atoi(argv[1]);
  // isomer
  Fullerene F;
  // out files
  string out_file_names[NFILE] = {"output/pp_anionic_clar_num",
                                  "output/pp_r_pent", "output/pp_r_hex",
                                  "output/pp_match_e", "output/pp_graph_id",
                                  "output/pp_full"};
  ofstream out_files_ptr[NFILE];
  open_out_file(p, out_file_names, out_files_ptr);

  // define gurobi solve environments
  GRBEnv grb_env = GRBEnv(true);
  grb_env.set(GRB_IntParam_OutputFlag, DEBUG_GUROBI);
  grb_env.start();

  F.id = 1;
  // while there are isomers to read in
  while (read_fullerene(F, p)) {
    // construct planar dual graph
    construct_planar_dual(F, p);
    // determine all p-anionic Clar structures
    p_anionic_clar_lp(F, p, grb_env, out_files_ptr);
#if DEBUG
    cout << "Graph number " << graph_num << endl;
    print_primal(F.n, F.primal);
    print_dual(F.dual_n, F.dual);
#endif
    F.id++;
  }
  close_files(out_files_ptr);
}
