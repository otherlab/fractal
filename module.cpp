//#####################################################################
// Module fractal
//#####################################################################
#include <other/core/python/module.h>
#include <other/core/mesh/TriangleMesh.h>
#include <other/core/openmesh/TriMesh.h>
#include <other/core/python/stl.h>
#include <other/core/structure/Hashtable.h>
#include <other/core/structure/Tuple.h>
#include <other/core/structure/UnionFind.h>
#include <other/core/utility/curry.h>
#include <other/core/utility/Log.h>
#include <other/core/vector/Frame.h>
#include <other/core/vector/Matrix4x4.h>
#include <other/core/vector/normalize.h>
#include <other/core/solver/powell.h>
#include <tr1/unordered_map>
#include <vector>
using namespace other;

typedef real T;
typedef Vector<T,2> TV2;
typedef Vector<T,3> TV3;
using std::tr1::unordered_map;
using std::vector;
using Log::cout;
using std::endl;

static Array<int> boundary_edges_to_faces(const TriangleMesh& mesh, RawArray<const Vector<int,2>> edges) {
  Array<int> face(edges.size(),false);
  auto incident = mesh.incident_elements();
  for (int s=0;s<edges.size();s++) {
    OTHER_ASSERT(incident.valid(edges[s][0]));
    for (int t : incident[edges[s][0]])
      if (mesh.elements[t].contains(edges[s][1])) {
        face[s] = t;
        goto found;
      }
    OTHER_ASSERT(false);
    found:;
  }
  return face;
}

static vector<Array<Vector<T,2>>> iterate_L_system(real start_angle_per_level, real shrink_factor, const string& axiom, const unordered_map<char,string>& rules, const unordered_map<char,double>& turns, const int levels) {
  OTHER_ASSERT(levels>=0);
  vector<Array<TV2>> curves;
  string pattern = axiom;
  for (int level : range(levels+1)) {
    // Trace curve
    TV2 X;
    const double step = pow(shrink_factor,level);
    double angle = level*start_angle_per_level;
    Array<TV2> curve;
    curve.append(X);
    for (char c : pattern) {
      if (c=='f') {
        X += step*polar(angle);
        curve.append(X);
      } else {
        auto it = turns.find(c);
        if (it != turns.end())
          angle += it->second;
      }
    }
    if (pattern[pattern.size()-1]=='c') {
      OTHER_ASSERT(magnitude(curve[0]-curve.back())<.01*step);
      curve.pop();
    }
    curves.push_back(curve);
    if (level == levels)
      break;

    // Refine
    string fine;
    for (char c : pattern) {
      auto it = rules.find(c);
      if (it != rules.end())
        fine += it->second;
      else
        fine += c;
    }
    swap(pattern,fine);
  }
  return curves;
}

static Ref<TriangleMesh> branching_mesh(const int branching, const int levels, const int base, bool closed) {
  OTHER_ASSERT(branching>=2);
  OTHER_ASSERT(levels>=1);
  const int64_t count = (base-!closed)*(1+branching)*(1-(int64_t)pow((double)branching,levels))/(1-branching);
  OTHER_ASSERT(0<count && count<(1<<30));
  Array<Vector<int,3>> tris;
  tris.preallocate(count);
  int n = base-!closed;
  int lo = 0;
  for (int level=0;level<levels;level++) {
    const int hi = lo+n+1-closed,
              lom = n+1-closed,
              him = branching*n+1-closed;
    for (int i : range(n)) {
      const int a = lo+i,
                ap = lo+(i+1)%lom;
      const int b = hi+branching*i;
      const int mid = branching/2;
      tris.append_assuming_enough_space(vec(ap,a,b+mid));
      for (int j : range(mid))
        tris.append_assuming_enough_space(vec(a,b+j,hi+(branching*i+j+1)%him));
      for (int j : range(mid,branching))
        tris.append_assuming_enough_space(vec(ap,b+j,hi+(branching*i+j+1)%him));
    }
    lo = hi;
    n *= branching;
  }
  OTHER_ASSERT(tris.size()==count);
  return new_<TriangleMesh>(tris);
}

namespace {
struct PatchInfo {
  bool boundary;
  Array<Vector<int,3>> tris;
  Array<TV3> X;
  Array<T> thick;
};
}

typedef vector<Tuple<Ref<TriangleMesh>,Array<TV3>,Array<T>,Array<Matrix<T,4>>>> Instances;

static void add_rotated_neighbors(RawArray<const int> neighbors, Hashtable<int,int>& vert_map, Array<int>& verts) {
  int start = -1;
  int score = numeric_limits<int>::max();
  for (int a=0;a<neighbors.size();a++) {
    int* s = vert_map.get_pointer(neighbors[a]);
    if (s && score>*s) {
      start = a;
      score = *s;
    }
  }
  OTHER_ASSERT(start>=0);
  for (int j : range(neighbors.size())) {
    const int a = neighbors[(start+j)%neighbors.size()];
    if (!vert_map.contains(a)) {
      vert_map.set(a,vert_map.size());
      verts.append(a);
    }
  }
}

// Generate one vertex per (triangle,vertex) pair, merging vertices according to edge connectivity
static Tuple<Ref<TriangleMesh>,Array<TV3>,Array<T>> make_manifold(const TriangleMesh& mesh, RawArray<const TV3> X, RawArray<const T> thick) {
  const auto adjacent_elements = mesh.adjacent_elements();
  UnionFind union_find(3*mesh.elements.size());
  for (int t0 : range(mesh.elements.size())) {
    const auto nodes0 = mesh.elements[t0];
    for (int i=0;i<3;i++) {
      const int t1 = adjacent_elements[t0][i];
      if (t1>=0) {
        const auto nodes1 = mesh.elements[t1];
        const int j = nodes1.find(nodes0[(i+1)%3]);
        OTHER_ASSERT(j>=0);
        union_find.merge(3*t0+i,3*t1+(j+1)%3);
        union_find.merge(3*t0+(i+1)%3,3*t1+j);
      }
    }
  }
  Array<int> map(union_find.size(),false);
  Array<TV3> X2;
  Array<T> thick2;
  for (int t : range(mesh.elements.size()))
    for (int i : range(3))
      if (union_find.is_root(3*t+i)) {
        map[3*t+i] = X2.append(X[mesh.elements[t][i]]);
        thick2.append(thick[mesh.elements[t][i]]);
      }
  Array<Vector<int,3>> tris2(mesh.elements.size(),false);
  for (int t : range(mesh.elements.size()))
    for (int i : range(3))
      tris2[t][i] = map[union_find.find(3*t+i)];
  return tuple(new_<TriangleMesh>(tris2),X2,thick2);
}

static Tuple<Array<int>,Instances,Instances> classify_loop_patches(const TriangleMesh& mesh, RawArray<const TV3> X, RawArray<const T> thickness, const int count) {
  const T tolerance = 1e-4;
  OTHER_ASSERT(count>=3);
  OTHER_ASSERT(mesh.nodes()==X.size());
  OTHER_ASSERT(mesh.nodes()==thickness.size());
  const int patches = mesh.elements.size()/count;
  OTHER_ASSERT(mesh.elements.size()==count*patches);
  vector<Tuple<PatchInfo,Array<Matrix<T,4>>>> reps;
  Array<int> names;
  const auto sorted_neighbors = mesh.sorted_neighbors();
  const auto adjacent_elements = mesh.adjacent_elements();
  const auto incident_elements = mesh.incident_elements();
  int boundary_count = 0;
  for (int p : range(patches)) {
    PatchInfo info;
    // Determine transform from the first triangle
    info.boundary = false;
    for (int i=0;i<count;i++)
      if (adjacent_elements[count*p+i].min()<0) {
        info.boundary = true;
        break;
      }
    const auto tri = mesh.elements[count*p];
    const TV3 x0 = X[tri.y];
    TV3 dx = X[tri.x]-x0;
    const T scale = normalize(dx),
            inv_scale = 1/scale;
    const TV3 dy = normalized((X[tri.z]-x0).projected_orthogonal_to_unit_direction(dx)),
              dz = cross(dx,dy);
    const auto transform = Matrix<T,4>::translation_matrix(x0)*Matrix<T,4>::from_linear(Matrix<T,3>(dx,dy,dz))*Matrix<T,4>::scale_matrix(scale);
    const auto inv_transform = transform.inverse();
    // Collect our vertices
    Hashtable<int> ours;
    Hashtable<int,int> vert_map;
    Array<int> verts;
    for (int t : count*p+range(count))
      for (int i : mesh.elements[t])
        if (ours.set(i))
          vert_map.set(i,verts.append(i));
    // Collect one ring vertices
    for (int i : range(verts.size()))
      add_rotated_neighbors(sorted_neighbors[verts[i]],vert_map,verts);
    // Collect two ring vertices adjacent to extraordinary vertices to account for the larger stencil of modified Loop subdivision
    for (int i : range(ours.size(),verts.size()))
      if (sorted_neighbors[verts[i]].size()!=6)
        add_rotated_neighbors(sorted_neighbors[verts[i]],vert_map,verts);
    // Compute signature
    for (int i : verts) {
      info.X.append(inv_transform.homogeneous_times(X[i]));
      info.thick.append(inv_scale*thickness[i]);
    }
    // Check for signature matches
    for (int r : range(reps.size())) {
      auto& rep = reps[r];
      if (rep.x.X.size()!=info.X.size())
        goto next;
      for (int i : range(info.X.size()))
        if (sqr_magnitude(rep.x.X[i]-info.X[i])>sqr(tolerance))
          goto next;
      rep.y.append(transform);
      names.append(r);
      goto found;
      next:;
    }
    {
      OTHER_ASSERT(reps.size()<10000);
      boundary_count += info.boundary;
      Array<Matrix<T,4>> transforms(1,false);
      transforms[0] = transform;
      names.append(reps.size());
      // Fill in triangles
      Array<int> tris;
      Hashtable<int> tri_set;
      for (int i=0;i<count;i++) {
        const int t = count*p+i;
        tris.append(t);
        tri_set.set(t);
      }
      for (int v : verts)
        for (int t : incident_elements[v])
          if (tri_set.set(t))
            tris.append(t);
      for (int t : tris) {
        Vector<int,3> tri;
        for (int i=0;i<3;i++) {
          const int* p = vert_map.get_pointer(mesh.elements[t][i]);
          if (!p)
            goto skip;
          tri[i] = vert_map.get(mesh.elements[t][i]);
        }
        info.tris.append(tri);
        skip:;
      }
      reps.push_back(tuple(info,transforms));
    }
    found:;
  }
  cout << "patches = "<<patches<<endl;
  cout << format("representatives = %d interior + %d boundary = %d total",reps.size()-boundary_count,boundary_count,reps.size())<<endl;

  // Extract one ring meshes for each representative
  Instances interior, boundary;
  for (int r : range(reps.size())) {
    const PatchInfo& info = reps[r].x;
    const auto fixed = make_manifold(new_<TriangleMesh>(info.tris),info.X,info.thick);
    OTHER_ASSERT(fixed.x->nodes()==fixed.y.size());
    OTHER_ASSERT(!fixed.x->nonmanifold_nodes(true).size());
    auto inst = tuple(fixed.x,fixed.y,fixed.z,reps[r].y);
    (info.boundary?boundary:interior).push_back(inst);
  }

  // All done
  return tuple(names,interior,boundary);
}

static T min_instance_dot(RawArray<const TV3> X, RawArray<const Matrix<T,4>> transforms, const TV3 up) {
  OTHER_ASSERT(X.size());
  T min_dot = inf;
  for (auto A : transforms) {
    const T upt = dot(up,A.translation());
    const TV3 Bup = A.linear().transpose_times(up);
    for (auto x : X)
      min_dot = min(min_dot,upt+dot(Bup,x));
  }
  return min_dot;
}

static TV3 shift_up(TV3 up, RawArray<const T> dup) {
  OTHER_ASSERT(dup.size()==2);
  auto r = Rotation<TV3>::from_rotated_vector(TV3(0,0,1),up);
  return (r*Rotation<TV3>::from_rotation_vector(TV3(dup[0],dup[1],0)))*TV3(0,0,1);
}

static T settling_energy(const vector<Tuple<Ref<TriangleMesh>,Array<const TV3>,Array<const Matrix<T,4>>>>* instances, TV3 up, RawArray<const T> dup) {
  up = shift_up(up,dup);
  // Slice as far down as possible
  T ground = inf;
  for (auto& inst : *instances)
    ground = min(ground,min_instance_dot(inst.y,inst.z,up));
  // Which way would we fall about the center?
  T energy = 0;
  for (auto& inst : *instances) {
    auto areas = inst.x->vertex_areas(inst.y);
    for (auto t : inst.z)
      for (int p=0;p<areas.size();p++)
        energy += areas[p]*(dot(up,t.homogeneous_times(inst.y[p]))-ground);
  }
  return energy;
}

static TV3 settle_instances(const vector<Tuple<Ref<TriangleMesh>,Array<const TV3>,Array<const Matrix<T,4>>>>& instances, const TV3 up, const T step) {
  Array<T> dup(2);
  powell(curry(settling_energy,&instances,up),dup,step,1e-5,0,20);
  return shift_up(up,dup);
}

static Tuple<Ref<TriangleMesh>,Array<TV3>> torus_mesh(const T R, const T r, const int N, const int n) {
  Array<TV3> X(N*n,false);
  Array<Vector<int,3>> tris(2*N*n,false);
  const T dA = 2*pi/N,
          da = 2*pi/n;
  for (int i : range(N))
    for (int j : range(n)) {
      const T u = dA*i, v = da*j;
      const T s = R+r*cos(v);
      const int I = i*n+j;
      X[I] = vec(s*cos(u),s*sin(u),r*sin(v));
      const int ii = (i+1)%n, jj = (j+1)%n;
      tris[2*I+0] = vec(i*n+j,ii*n+j,ii*n+jj);
      tris[2*I+1] = vec(i*n+j,ii*n+jj,i*n+jj);
    }
  return tuple(new_<TriangleMesh>(tris),X);
}

OTHER_PYTHON_MODULE(other_fractal) {
  OTHER_FUNCTION(boundary_edges_to_faces)
  OTHER_FUNCTION(iterate_L_system)
  OTHER_FUNCTION(branching_mesh)
  OTHER_FUNCTION(classify_loop_patches)
  OTHER_FUNCTION(min_instance_dot)
  OTHER_FUNCTION(settle_instances)
  OTHER_FUNCTION(torus_mesh)
  OTHER_FUNCTION(make_manifold)
}
