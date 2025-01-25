#include "include.h"

// For error reporting
void throw_error(const int n, const int p, const int graph_id,
                 string error_message) {
  throw runtime_error("\nError: n = " + to_string(n) + ", p = " + to_string(p) +
                      ", graph num = " + to_string(graph_id) + error_message);
}

// read in fullerene graph and populate some default values into the
// data structure
bool read_fullerene(Fullerene(&F), const int p) {
  int degree, n;
  // read in number of vertices of isomer
  if (!(cin >> n))
    return false;
  // check that number of vertices is valid
  // n should be an even number at least 20 and not equal to 22
  if (n < 20 || n == 22 || n % 2 != 0) {
    throw_error(F.n, p, F.id, "\nInvalid fullerene size " + to_string(n));
  }
  // resize the fullerene
  F.Resize(n);
  // for each vertex in the graph
  for (int i = 0; i < n; i++) {
    // read in degree of vertex i
    if (!(cin >> degree)) {
      throw_error(F.n, p, F.id,
                  "\nFailed reading vertex " + to_string(i) + "'s degree");
    }
    // check vertex degree
    if (degree != 3) {
      const string msg = "\nVertex " + to_string(i) +
                         " has invalid degree: " + to_string(degree);
      throw_error(F.n, p, F.id, msg);
    }
    // for each neighbour of i
    for (int j = 0; j < 3; j++) {
      // update adjacency list of vertex i
      if (!(cin >> F.primal[i].adj_v[j])) {
        const string msg = "\nFailed reading neighbour " + to_string(j) +
                           " of vertex " + to_string(i);
        throw_error(F.n, p, F.id, msg);
      }
      // each vertex lies on 3 faces, at this time, we do not know what their
      // ids are, we will therefore assign then as -1 to represent 'unassigned'
      F.primal[i].faces[j] = -1;
    }
  }
  return true;
}

void print_primal(const int n, const vector<vertex> primal) {
  cout << "Primal graph" << endl;
  cout << "Vert:  Neighbours         Faces           Edges" << endl;
  for (int i = 0; i < n; i++) {
    cout << setw(4) << i << ": ";
    for (int j = 0; j < 3; j++) {
      cout << setw(4) << primal[i].adj_v[j];
      if (j == 2)
        cout << " || ";
    }
    for (int j = 0; j < 3; j++) {
      cout << setw(4) << primal[i].faces[j];
      if (j == 2)
        cout << " || ";
    }
    for (int j = 0; j < 3; j++) {
      cout << setw(4) << primal[i].edges[j];
    }
    cout << endl;
  }
}

void print_dual(const int dual_n, const vector<face> dual) {
  cout << "Dual graph" << endl << "Face:       Neighbours" << endl;
  for (int i = 0; i < dual_n; i++) {
    cout << setw(4) << i << ": ";
    for (int j = 0; j < dual[i].size; j++) {
      cout << setw(4) << dual[i].adj_f[j];
    }
    cout << endl;
  }
}

void print_sol(const Fullerene(&F), const int num_res_faces,
               const vector<GRBVar> fvars, const vector<GRBVar> evars) {
  cout << endl
       << "n = " << F.n << ", graph num = " << F.id << endl
       << "Solution summary:" << endl;
  // if there is no solution
  if (num_res_faces == 0) {
    cout << "No solution" << endl;
    return;
  }
  cout << "Num resonant faces: " << num_res_faces << endl;
  cout << "Resonant face ids: ";
  for (int i = 0; i < F.dual_n; i++) {
    if (fvars[i].get(GRB_DoubleAttr_X) < 0.99)
      continue;
    cout << i << " ";
  }
  cout << endl << "Matching edges: ";
  for (int i = 0; i < F.num_edges; i++) {
    if (evars[i].get(GRB_DoubleAttr_X) < 0.99)
      continue;
    cout << "(" << F.edges[i].vertices[0] << "," << F.edges[i].vertices[1]
         << ") ";
  }
  cout << endl << endl;
}

void save_sol(const Fullerene(&F), const int p, const int num_res_faces,
              const vector<GRBVar> fvars, const vector<GRBVar> evars,
              ofstream out_files_ptr[NFILE]) {
  // record graph number
  out_files_ptr[4] << F.id << endl;

  out_files_ptr[5] << F.n;
  for (int i = 0; i < F.dual_n; i++) {
    if (F.dual[i].size == 5) {
      out_files_ptr[5] << " " << i;
    }
  }
  out_files_ptr[5] << endl;

  // if there is no solution
  if (num_res_faces == 0) {
    out_files_ptr[0] << 0 << endl;
    out_files_ptr[1] << 0 << endl;
    out_files_ptr[2] << 0 << endl;
    out_files_ptr[3] << 0 << endl;
    return;
  }
  // number of resonant faces
  out_files_ptr[0] << num_res_faces << endl;
  // print resonant faces
  out_files_ptr[1] << p;
  out_files_ptr[2] << num_res_faces - p;
  for (int i = 0; i < F.dual_n; i++) {
    if (fvars[i].get(GRB_DoubleAttr_X) < 0.99)
      continue;
    if (F.dual[i].size == 5) {
      out_files_ptr[1] << " " << i;
    } else {
      out_files_ptr[2] << " " << i;
    }
  }
  out_files_ptr[1] << endl;
  out_files_ptr[2] << endl;
  // matching edges
  out_files_ptr[3] << F.n - 6 * (num_res_faces - p) - 5 * p;
  for (int i = 0; i < F.num_edges; i++) {
    if (evars[i].get(GRB_DoubleAttr_X) < 0.99)
      continue;
    out_files_ptr[3] << " " << F.edges[i].vertices[0] << " "
                     << F.edges[i].vertices[1];
  }
  out_files_ptr[3] << endl;
}

void get_out_name(const int p, string &fname) {
  // determine first digit of p
  int p1 = p / 10;
  // determine last digit of p
  int p2 = p % 10;
  // update the file name characters
  fname[7] = (char)p1 + '0';
  fname[8] = (char)p2 + '0';
}

void open_out_file(const int p, string (&out_file_names)[NFILE],
                   ofstream out_files_ptr[NFILE]) {
  for (int i = 0; i < NFILE; i++) {
    get_out_name(p, out_file_names[i]);
    out_files_ptr[i].open(out_file_names[i], ios::app);
    if (!out_files_ptr[i].is_open()) {
      throw runtime_error("\nError: Could not open file " + out_file_names[i]);
    }
  }
}

void close_files(ofstream out_files_ptr[NFILE]) {
  for (int i = 0; i < NFILE; i++) {
    out_files_ptr[i].close();
  }
}
