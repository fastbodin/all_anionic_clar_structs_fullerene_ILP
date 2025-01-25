#include "include.h"


int check_if_sol_valid(const Fullerene(&F), const int p,
                       const vector<GRBVar> fvars, const vector<GRBVar> evars) {
  int num_res_faces = 0, res_pents = 0;
  string msg;
  // for each vertex in the graph
  for (int i = 0; i < F.n; i++) {
    // they should be covered by the p-anionic Clar structure exactly once
    int covered = 0;
    for (int j = 0; j < 3; j++) {
      // note the tolerance given to the variable assignment, i.e. > 0.99
      // covered by matching edge
      if (evars[F.primal[i].edges[j]].get(GRB_DoubleAttr_X) > 0.99)
        covered++;
      // covered by resonant face
      if (fvars[F.primal[i].faces[j]].get(GRB_DoubleAttr_X) > 0.99)
        covered++;
    }
    if (covered != 1) {
      msg = "\nVertex " + to_string(i) + " is covered " + to_string(covered) +
            " times.";
      throw_error(F.n, p, F.id, msg);
    }
  }
  // for each face in graph
  for (int i = 0; i < F.dual_n; i++) {
    // if face is resonant
    if (fvars[i].get(GRB_DoubleAttr_X) > 0.99) {
      num_res_faces++;
      if (F.dual[i].size == 5)
        res_pents += 1;
    }
  }
  if (res_pents != p) {
    msg = "\nIncorrect # of res. pents: " + to_string(res_pents);
    throw_error(F.n, p, F.id, msg);
  }
  return num_res_faces;
}

int assess_rec_solve(const Fullerene(&F), const int p, GRBModel(&model),
                     vector<GRBVar>(&fvars), vector<GRBVar>(&evars),
                     ofstream out_files_ptr[NFILE], int opt_val) {
  const int optimstatus = model.get(GRB_IntAttr_Status);
  // if optimum is attained
  if (optimstatus == GRB_OPTIMAL) {
    // check solution and grab number of resonant faces
    const int num_res_faces = check_if_sol_valid(F, p, fvars, evars);
    // only save the solution if it has correct number of resonant faces
    if (num_res_faces == opt_val) {
      save_sol(F, p, num_res_faces, fvars, evars, out_files_ptr);
    }
#if DEBUG_CLAR
    print_sol(F, num_res_faces, fvars, evars);
#endif
    return num_res_faces;
    // if there is no solution
  } else if (optimstatus == GRB_INFEASIBLE) {
    return 0;
    // error handling
  } else {
    const string msg = "\nStatus of solve: " + to_string(optimstatus) +
                       "\nCheck Gurobi Optimization Status Codes";
    throw_error(F.n, p, F.id, msg);
    return -1;
  }
}

int assess_solve(const Fullerene(&F), const int p, GRBModel(&model),
                 vector<GRBVar>(&fvars), vector<GRBVar>(&evars),
                 ofstream out_files_ptr[NFILE]) {
  int optimstatus = model.get(GRB_IntAttr_Status);
  // if optimum is attained
  if (optimstatus == GRB_OPTIMAL) {
    // check solution and grab number of resonant faces
    int num_res_faces = check_if_sol_valid(F, p, fvars, evars);
    save_sol(F, p, num_res_faces, fvars, evars, out_files_ptr);
#if DEBUG_CLAR
    print_sol(F, num_res_faces, fvars, evars);
#endif
    return num_res_faces;
    // if there is no solution
  } else if (optimstatus == GRB_INFEASIBLE) {
    // there are 0 resonant faces since no valid solution
    save_sol(F, p, 0, fvars, evars, out_files_ptr);
#if DEBUG_CLAR
      print_sol(F, 0, fvars, evars);
#endif
    return 0;
  } else {
    const string msg = "\nStatus of solve: " + to_string(optimstatus) +
                       "\nCheck Gurobi Optimization Status Codes";
    throw_error(F.n, p, F.id, msg);
    return -1;
  }
}

void exlude_previous_sol(const Fullerene(&F), const int p, GRBModel(&model),
                         vector<GRBVar>(&fvars), vector<GRBVar>(&evars),
                         vector<int> res_faces, vector<int> match_edges) {
  // record solution
  int vec_index = 0;
  for (int f = 0; f < F.dual_n; f++) {
    // if face is resonat
    if (fvars[f].get(GRB_DoubleAttr_X) > 0.99) {
      res_faces[vec_index++] = f;
    }
  }
  vec_index = 0;
  for (int e = 0; e < F.num_edges; e++) {
    // if edge is matching edge
    if (evars[e].get(GRB_DoubleAttr_X) > 0.99) {
      match_edges[vec_index++] = e;
    }
  }
  // reset the model to add a new constraint
  model.reset();
  // constraint to exclude old solutions
  GRBLinExpr cons = 0;
  for (int f = 0; f < res_faces.size(); f++) {
    cons += fvars[res_faces[f]];
  }
  for (int e = 0; e < match_edges.size(); e++) {
    cons += evars[match_edges[e]];
  }
  model.addConstr(cons <= res_faces.size() + match_edges.size() - 1);
}

void add_Clar_num_cons(const Fullerene(&F), const int opt_val, GRBModel(&model),
                       vector<GRBVar>(&fvars)) {
  // add known p-anionic Clar number as constraint
  GRBLinExpr cons = 0;
  for (int f = 0; f < F.dual_n; f++) {
    cons += fvars[f];
  }
  // add constraint to model
  model.addConstr(cons == opt_val);
}

int solve_all_structs(const Fullerene(&F), const int p, GRBModel(&model),
                      vector<GRBVar>(&fvars), vector<GRBVar>(&evars),
                      ofstream out_files_ptr[NFILE], const int opt_val) {

  const int num_match_e = (F.n-5*p-6*(opt_val-p))/2;
  // vectors to hold previous solutions
  vector<int> res_f(opt_val);
  vector<int> match_e(num_match_e);

  // first remove the previous solution
  exlude_previous_sol(F, p, model, fvars, evars, res_f, match_e);
  // use the fact that we know the p-anionic Clar number as a constraint
  add_Clar_num_cons(F, opt_val, model, fvars);
  // begin the solve for all p-anionic Clar structures
  int cur_val = opt_val;
  // while correct number of resonant faces
  while (cur_val == opt_val) {
    // solve for next solution
    model.optimize();
    // check validity of new solution and whether it has the right
    // number of resonant faces
    cur_val = assess_rec_solve(F, p, model, fvars, evars, out_files_ptr, 
                               opt_val);
    // if the solution achieves the correct number of resonat faces,
    // remove it an prepare to solve the next
    if (cur_val == opt_val) {
      exlude_previous_sol(F, p, model, fvars, evars, res_f, match_e);
    }
  }
  return opt_val;
}

void add_cons(const Fullerene(&F), const int p, GRBModel(&model),
              vector<GRBVar>(&fvars), vector<GRBVar>(&evars)) {
  // each vertex is either in a resonant face or is the endpoint of
  // a matching edge
  for (int i = 0; i < F.n; i++) {
    GRBLinExpr cons1 = 0;
    // each vertex lies on exactly three faces and is the endpoint of
    // exactly three edges
    for (int j = 0; j < 3; j++) {
#if DEBUG_CLAR
      cout << i << " is endpoint of edge " << F.primal[i].edges[j] << endl;
      cout << i << " lies on face " << F.primal[i].faces[j] << endl;
#endif
      cons1 += evars[F.primal[i].edges[j]] + fvars[F.primal[i].faces[j]];
    }
    // add the constraint to the model
    model.addConstr(cons1 == 1);
  }
  // need p resonant pentagons
  GRBLinExpr cons2 = 0;
  for (int f = 0; f < F.dual_n; f++) {
    if (F.dual[f].size == 5)
      cons2 += fvars[f];
  }
  // add constraint to model
  model.addConstr(cons2 == p);
}

void add_vars(const Fullerene(&F), const int p, GRBModel(&model),
              vector<GRBVar>(&fvars), vector<GRBVar>(&evars)) {
  // make face variables
  for (int f = 0; f < F.dual_n; f++) {
    // lower bound, upper bound, objective coeff, type
    fvars[f] = model.addVar(0.0, 1.0, 1.0, GRB_BINARY);
  }
  // make edge variables
  for (int i = 0; i < F.num_edges; i++) {
    // lower bound, upper bound, objective coeff, type
    evars[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
  }
}

int p_anionic_clar_lp(const Fullerene(&F), const int p, GRBEnv grb_env,
                      ofstream out_files_ptr[NFILE]) {
#if DEBUG_CLAR
  cout << "n = " << F.n << ", p = " << p << ", graph num = " << F.id << endl;
  cout << "Solving LP" << endl;
#endif
  try {
    // create an empty model
    GRBModel model = GRBModel(grb_env);
    // The objective is to maximize number of resonant faces
    model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
    // get the solver to focus on finding integer solutions
    //model.set(GRB_IntParam_IntegralityFocus, 1);

    // fvars[f] = 1 if the face f is resonant and 0 otherwise
    vector<GRBVar> fvars(F.dual_n);
    // evars[i] = 1 if the edge i is matching edge and 0 otherwise
    vector<GRBVar> evars(F.num_edges);

    // add variables to model
    // note that the obj. coefficients are set during this step 
    add_vars(F, p, model, fvars, evars);
    // add constraints to model
    add_cons(F, p, model, fvars, evars);
    // run model
    model.optimize();
    //asses solve
    int opt_val = assess_solve(F, p, model, fvars, evars, out_files_ptr);
    // if there exists no p-anionic Clar structure, stop here
    if (opt_val < 1) return 0;
    // otherwise, find all other p-anionic Clar structures
    return solve_all_structs(F, p, model, fvars, evars, out_files_ptr, opt_val);
  // error handling
  } catch (GRBException e) {
    const string msg = "\nCode: " + to_string(e.getErrorCode()) +
                       "\nMessage: " + e.getMessage();
    throw_error(F.n, p, F.id, msg);
  } catch (runtime_error e) {
    throw runtime_error(e);
  } catch (...) {
    throw_error(F.n, p, F.id, "\nUnknown error during optimization");
  }
  return -1;
}
