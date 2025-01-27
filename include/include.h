#include "gurobi_c++.h"
#include <fstream>
#include <iomanip>

using namespace std;

// For debugging purposes
#define DEBUG 0
#define DEBUG_DUAL 0
#define DEBUG_CLAR 0
#define DEBUG_GUROBI 0

//------ DO NOT CHANGE BELOW ------
//-------- HERE BE DRAGONS --------

// number of out files
constexpr int NFILE = 4;

// information on each vertex
struct vertex {
  // vertices adjacent to it
  int adj_v[3];
  // faces it lies on
  int faces[3];
  // edges it is an endpoint of
  int edges[3];
};

// information on each face
struct face {
  // length of face (will be 5 or 6)
  int size;
  // faces it shares an edge with
  int adj_f[6];
  // vertices on face, note they will be recorded in counter clockwise order
  int vertices[6];
};

// information on each edge
struct edge {
  // vertices that form its endpoints
  // vertices[0] < vertices[1]
  int vertices[2];
};

// information on each fullerene isomer
class Fullerene {
public:
  Fullerene() : n(0), dual_n(0), num_edges(0), id(0) {}
  // resize fullerene
  void Resize(int num_vertices) {
    n = num_vertices;
    // Fullerenes are 3-regular graphs, therefore a fullerene on n vertices has
    // 3n/2 edges
    num_edges = (3 * num_vertices / 2);
    // By Euler's formula, the number of faces in a fullerene is 3n/2 - n + 2
    dual_n = (num_vertices / 2 + 2);
    // resize vectors appropriately
    primal.resize(num_vertices);
    dual.resize(dual_n);
    edges.resize(num_edges);
  }
  // Attributes
  int n, dual_n, num_edges,
      id;                // size of primal, size of dual, # of edges, and id
  vector<vertex> primal; // planar graph information
  vector<face> dual;     // planar dual graph information
  vector<edge> edges;    // edge information
};

// From read_and_print.cpp
void throw_error(const int n, const int p, const int graph_id,
                 string error_message);
bool read_fullerene(Fullerene(&F), const int p);
void print_primal(const int n, const vector<vertex> primal);
void print_dual(const int dual_n, const vector<face> dual);
void print_sol(const Fullerene(&F), const int num_res_faces,
               const vector<GRBVar> fvars, const vector<GRBVar> evars);
void save_sol(const Fullerene(&F), const int p, const int num_res_faces,
              const vector<GRBVar> fvars, const vector<GRBVar> evars,
              ofstream out_files_ptr[NFILE]);
void get_out_name(const int p, string &fname);
void open_out_file(const int p, string (&out_file_names)[NFILE],
                   ofstream out_files_ptr[NFILE]);
void close_files(ofstream out_files_ptr[NFILE]);

// From dual.cpp
int find_position(const int v, const int u, int v_adj[3]);
int counter_clockwise_walk(const int face_id, int u, int v, const int n,
                           vector<vertex>(&primal), face(&cur_face));
void construct_planar_dual(Fullerene(&F), const int p);

// From lp.cpp
int check_if_sol_valid(const Fullerene(&F), const int p,
                       const vector<GRBVar> fvars, const vector<GRBVar> evars);
int p_anionic_clar_lp(const Fullerene(&F), const int p, GRBEnv grb_env,
                      ofstream out_files_ptr[NFILE]);
