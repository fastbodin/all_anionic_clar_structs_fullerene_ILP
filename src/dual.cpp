#include "include.h"
#include <stdexcept>

int find_position(const int v, const int u, int v_adj[3]) {
#if DEBUG_DUAL
  cout << "Looking for neighbour " << u << " of vertex " << v << endl;
#endif
  for (int i = 0; i < 3; i++) {
    if (v_adj[i] == u)
      return i;
  }
  throw runtime_error("\nCould not find position of vertex " + to_string(u) +
                      " in neighbourhood of vertex " + to_string(v));
}

int counter_clockwise_walk(const int face_id, int u, int v, const int n,
                           vector<vertex>(&primal), face(&cur_face)) {
  // find the position of v in u's adj list, remember this list is in clockwise
  // order in a planar embedding
  int pos = find_position(u, v, primal[u].adj_v);
  // initialize the face size
  int face_size = 0;
  // w will be used to walk around the face
  int w;

  // while there exists a vertex on the face that has not recorded the face id
  while (primal[u].faces[pos] == -1) {
    // record u as a vertex of face we are walking the boundary of counter
    // clockwise
    cur_face.vertices[face_size++] = u;
    // record the face id for vertex u, note the position is unique identified
    // by v's position in u's adj list
    primal[u].faces[pos] = face_id;
    // we want to grab the next vertex on the face during our counter clockwise
    // walk, first we grab the position of u in v's adj list.
    try {
      pos = find_position(v, u, primal[v].adj_v);
    } catch (runtime_error e) {
      throw runtime_error(e);
    }
    // move to the next position (modulo 3) in v's adj list to get the next
    // vertex
    pos = (pos + 1) % 3;
    // record w as this next vertex
    w = primal[v].adj_v[pos];
    // shift along vertices and repeat loop
    u = v;
    v = w;
    // this will terminate once v becomes the original u handed to the function
  }
  // lets check that face is pentagon or hexagon
  if (face_size != 5 && face_size != 6) {
    throw runtime_error("\nError: face " + to_string(face_id) + " has size " +
                        to_string(face_size));
  }
#if DEBUG_DUAL
  cout << face_id << " has size " << face_size << endl;
#endif
  return face_size;
}

void construct_planar_dual(Fullerene(&F), const int p) {
  // When we read in the fullerene, we had not yet assigned face ids. Therefore,
  // we set the face ids at each vertex as -1 (to represent unassigned). We
  // will now assign them by constructing each planar face
#if DEBUG_DUAL
  cout << "n = " << F.n << ", p = " << p << ", graph num = " << F.id << endl;
  cout << "Constructing planar dual" << endl;
#endif

  int face_id = 0, edge_id = 0, u;
  // for each vertex v
  for (int v = 0; v < F.n; v++) {
    // and each face it lies on
    for (int j = 0; j < 3; j++) {
      // u is the jth neighbour of vertex v in primal
      u = F.primal[v].adj_v[j];
      // record the edge exactly once
      if (v < u) {
        // u is the jth neighbour of v
        F.primal[v].edges[j] = edge_id;
        // find position of v in u's neighbourhood
        try {
          F.primal[u].edges[find_position(u, v, F.primal[u].adj_v)] = edge_id;
        } catch (runtime_error e) {
          throw_error(F.n, p, F.id, e.what());
        }
        // record vertices of given edge
        F.edges[edge_id].vertices[0] = v;
        F.edges[edge_id].vertices[1] = u;
        edge_id++;
      }
      // if face is unassigned
      if (F.primal[v].faces[j] == -1) {
        // lets walk the face containing u and v
        try {
          F.dual[face_id].size = counter_clockwise_walk(
              face_id, v, u, F.n, F.primal, F.dual[face_id]);
        } catch (runtime_error e) {
          throw_error(F.n, p, F.id, e.what());
        }
        face_id++;
      }
    }
  }
  // record the number of edges in the graph
  if (edge_id != 3 * F.n / 2) {
    throw_error(F.n, p, F.id, "\nIncorrect # of edges: " + to_string(edge_id));
  }
  F.num_edges = edge_id;
  // record the number of faces in the planar dual
  F.dual_n = face_id;

  // by this point, we have determined which vertices and in which face by face
  // id we now want to determine the adjacency between faces in the planar dual
  int v, face_g, face_f_size, pos;
  // for each face f
  for (int f = 0; f < face_id; f++) {
    face_f_size = F.dual[f].size;
    // for each vertex on face i
    for (int j = 0; j < face_f_size; j++) {
      // uv is an edge on face i
      u = F.dual[f].vertices[j];
      v = F.dual[f].vertices[(j + 1) % face_f_size];
      // get the position of u in v's adj list
      try {
        pos = find_position(v, u, F.primal[v].adj_v);
      } catch (runtime_error e) {
        throw_error(F.n, p, F.id, e.what());
      }
      // face_g is the face that v and u lie on that is not equal to face f
      face_g = F.primal[v].faces[pos];
      // update adj list of face f
      F.dual[f].adj_f[j] = face_g;
#if DEBUG_DUAL
      cout << "Face " << setw(4) << f << " is adj to face " << setw(4) << face_g
           << endl;
#endif
    }
  }
#if DEBUG_DUAL
  cout << endl;
#endif
}
