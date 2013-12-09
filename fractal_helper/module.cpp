//#####################################################################
// Module fractal
//#####################################################################
#include <geode/array/NdArray.h>
#include <geode/array/IndirectArray.h>
#include <geode/geometry/ParticleTree.h>
#include <geode/geometry/Ray.h>
#include <geode/geometry/SimplexTree.h>
#include <geode/geometry/traverse.h>
#include <geode/geometry/Triangle3d.h>
#include <geode/mesh/SegmentSoup.h>
#include <geode/mesh/TriangleSoup.h>
#include <geode/openmesh/TriMesh.h>
#include <geode/python/Class.h>
#include <geode/python/module.h>
#include <geode/python/stl.h>
#include <geode/python/wrap.h>
#include <geode/solver/powell.h>
#include <geode/structure/Hashtable.h>
#include <geode/structure/Tuple.h>
#include <geode/structure/UnionFind.h>
#include <geode/utility/curry.h>
#include <geode/utility/Log.h>
#include <geode/utility/tr1.h>
#include <geode/vector/Frame.h>
#include <geode/vector/Matrix4x4.h>
#include <geode/vector/normalize.h>
#include <vector>
using namespace geode;

typedef real T;
typedef Vector<T,2> TV2;
typedef Vector<T,3> TV3;
using std::vector;
using Log::cout;
using std::endl;

static Array<int> boundary_edges_to_faces(const TriangleSoup& mesh, RawArray<const Vector<int,2>> edges) {
  Array<int> face(edges.size(),false);
  auto incident = mesh.incident_elements();
  for (int s=0;s<edges.size();s++) {
    GEODE_ASSERT(incident.valid(edges[s][0]));
    for (int t : incident[edges[s][0]])
      if (mesh.elements[t].contains(edges[s][1])) {
        face[s] = t;
        goto found;
      }
    GEODE_ASSERT(false);
    found:;
  }
  return face;
}

static vector<Array<TV2>> iterate_L_system(NdArray<const T> start_angle_per_level, const T shrink_factor, const string& axiom, const unordered_map<char,string>& rules, const unordered_map<char,double>& turns, const int levels) {
  GEODE_ASSERT(levels>=0);
  GEODE_ASSERT(start_angle_per_level.rank()<=1 && start_angle_per_level.flat.size());
  vector<Array<TV2>> curves;
  string pattern = axiom;
  double start_angle = 0;
  for (int level : range(levels+1)) {
    // Trace curve
    TV2 X;
    const double step = pow(shrink_factor,level);
    double angle = start_angle;
    start_angle += start_angle_per_level.flat[level%start_angle_per_level.flat.size()];
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
      GEODE_ASSERT(magnitude(curve[0]-curve.back())<.01*step);
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

static Ref<TriangleSoup> branching_mesh(const int branching, const int levels, const int base, bool closed) {
  GEODE_ASSERT(branching>=2);
  GEODE_ASSERT(levels>=1);
  const int64_t count = (base-!closed)*(1+branching)*(1-(int64_t)pow((double)branching,levels))/(1-branching);
  GEODE_ASSERT(0<count && count<(1<<30));
  Array<Vector<int,3>> tris;
  tris.preallocate(int(count));
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
  GEODE_ASSERT(tris.size()==count);
  return new_<TriangleSoup>(tris);
}

namespace {
struct PatchInfo {
  bool boundary;
  Array<Vector<int,3>> tris;
  Array<TV3> X;
  Array<T> thick;
};
}

typedef vector<Tuple<Ref<TriangleSoup>,Array<TV3>,Array<T>,Array<Matrix<T,4>>>> Instances;

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
  GEODE_ASSERT(start>=0);
  for (int j : range(neighbors.size())) {
    const int a = neighbors[(start+j)%neighbors.size()];
    if (!vert_map.contains(a)) {
      vert_map.set(a,vert_map.size());
      verts.append(a);
    }
  }
}

// Generate one vertex per (triangle,vertex) pair, merging vertices according to edge connectivity
static Tuple<Ref<TriangleSoup>,Array<TV3>,Array<T>> make_manifold(const TriangleSoup& mesh, RawArray<const TV3> X, RawArray<const T> thick) {
  const auto adjacent_elements = mesh.adjacent_elements();
  UnionFind union_find(3*mesh.elements.size());
  for (int t0 : range(mesh.elements.size())) {
    const auto nodes0 = mesh.elements[t0];
    for (int i=0;i<3;i++) {
      const int t1 = adjacent_elements[t0][i];
      if (t1>=0) {
        const auto nodes1 = mesh.elements[t1];
        const int j = nodes1.find(nodes0[(i+1)%3]);
        GEODE_ASSERT(j>=0);
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
  return tuple(new_<TriangleSoup>(tris2),X2,thick2);
}

static Tuple<Array<int>,Instances,Instances> classify_loop_patches(const TriangleSoup& mesh, RawArray<const TV3> X, RawArray<const T> thickness, const int count, const bool two_ring) {
  const T tolerance = 1e-4;
  GEODE_ASSERT(count>=3);
  GEODE_ASSERT(mesh.nodes()==X.size());
  GEODE_ASSERT(mesh.nodes()==thickness.size());
  const int patches = mesh.elements.size()/count;
  GEODE_ASSERT(mesh.elements.size()==count*patches);
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
      if (two_ring || sorted_neighbors[verts[i]].size()!=6)
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
      GEODE_ASSERT(reps.size()<10000);
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
    const auto fixed = make_manifold(new_<TriangleSoup>(info.tris),info.X,info.thick);
    GEODE_ASSERT(fixed.x->nodes()==fixed.y.size());
    GEODE_ASSERT(!fixed.x->nonmanifold_nodes(true).size());
    auto inst = tuple(fixed.x,fixed.y,fixed.z,reps[r].y);
    (info.boundary?boundary:interior).push_back(inst);
  }

  // All done
  return tuple(names,interior,boundary);
}

static T min_instance_dot(RawArray<const TV3> X, RawArray<const Matrix<T,4>> transforms, const TV3 up) {
  GEODE_ASSERT(X.size());
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
  GEODE_ASSERT(dup.size()==2);
  auto r = Rotation<TV3>::from_rotated_vector(TV3(0,0,1),up);
  return (r*Rotation<TV3>::from_rotation_vector(TV3(dup[0],dup[1],0)))*TV3(0,0,1);
}

static T settling_energy(const vector<Tuple<Ref<TriangleSoup>,Array<const TV3>,Array<const Matrix<T,4>>>>* instances, TV3 up, RawArray<const T> dup) {
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

static TV3 settle_instances(const vector<Tuple<Ref<TriangleSoup>,Array<const TV3>,Array<const Matrix<T,4>>>>& instances, const TV3 up, const T step) {
  Array<T> dup(2);
  powell(curry(settling_energy,&instances,up),dup,step,1e-5,0,20);
  return shift_up(up,dup);
}

static Tuple<Ref<TriangleSoup>,Array<TV3>> torus_mesh(const T R, const T r, const int N, const int n) {
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
  return tuple(new_<TriangleSoup>(tris),X);
}

static Array<const int> boundary_curve_at_height(const TriangleSoup& mesh, RawArray<const TV3> X, const T z) {
  const T tolerance = 1e-5;
  Array<Vector<int,2>> curve;
  for (const auto s : mesh.boundary_mesh()->elements)
    if (abs(X[s.x].z-z)<tolerance && abs(X[s.y].z-z)<tolerance)
      curve.append(s);
  const auto curves = new_<SegmentSoup>(curve)->polygons();
  GEODE_ASSERT(curves.x.size()==0 && curves.y.size()==1);
  return curves.y.flat;
}

static Array<TV3> unit_spring_energy_gradient(RawArray<const Vector<int,2>> edges, RawArray<const T> rest,
                                              RawArray<const TV3> X) {
  GEODE_ASSERT(edges.size()==rest.size());
  Array<TV3> gradient(X.size());
  for (const int s : range(edges.size())) {
    const auto edge = edges[s];
    TV3 dX = X[edge.y]-X[edge.x];
    const T len = normalize(dX);
    const TV3 dE = (len-rest[s])*dX;
    gradient[edge.x] -= dE;
    gradient[edge.y] += dE;
  }
  return gradient;
}

namespace {
struct SimpleCollisions : public Object {
  GEODE_DECLARE_TYPE(GEODE_NO_EXPORT)

  const Array<TV3> X;
  const Ref<ParticleTree<TV3>> points;
  const Ref<SimplexTree<TV3,1>> edges;
  const Ref<SimplexTree<TV3,2>> faces;
  const T close;
  const bool include_faces;

private:
  SimpleCollisions(const TriangleSoup& mesh, RawArray<const TV3> X0, const T close, const bool include_faces)
    : X(X0.copy())
    , points(new_<ParticleTree<TV3>>(X,1))
    , edges(new_<SimplexTree<TV3,1>>(mesh.segment_soup(),X,1))
    , faces(new_<SimplexTree<TV3,2>>(mesh,X,1))
    , close(close)
    , include_faces(include_faces) {}
public:

  T distance_energy(const T d) const {
    return d<close ? .5*sqr(d-close) : 0;
  }

  T distance_gradient(const T d) const {
    return d<close ? d-close : 0;
  }

  typedef Vector<int,2> IV2;
  typedef Vector<int,3> IV3;
  typedef Segment<TV3> Seg;
  typedef Triangle<TV3> Tri;

  template<class Visit> struct EEVisitor {
    const SimplexTree<TV3,1>& edges;
    RawArray<const Vector<int,2>> elements;
    RawArray<const TV3> X;
    const Visit& visit;

    EEVisitor(const SimplexTree<TV3,1>& edges, RawArray<const TV3> X, const Visit& visit)
      : edges(edges), elements(edges.mesh->elements), X(X), visit(visit) {}

    template<class... A> bool cull(A...) const {
      return false;
    }

    void leaf(const int n0) const {} // Edges do not intersect themselves

    void leaf(const int n0, const int n1) const {
      const auto ei = elements[edges.prims(n0)[0]],
                 ej = elements[edges.prims(n1)[0]];
      if (!(ei.contains(ej.x) || ei.contains(ej.y)))
        visit(ei,ej,simplex(X[ei.x],X[ei.y]),simplex(X[ej.x],X[ej.y]));
    }
  };

  template<class Visit> struct PFVisitor {
    const ParticleTree<TV3>& points;
    const SimplexTree<TV3,2>& faces;
    RawArray<const Vector<int,3>> elements;
    RawArray<const TV3> X;
    const Visit& visit;

    PFVisitor(const ParticleTree<TV3>& points, const SimplexTree<TV3,2>& faces, RawArray<const TV3> X, const Visit& visit)
      : points(points), faces(faces), elements(faces.mesh->elements), X(X), visit(visit) {}

    bool cull(const int n0, const int n1) const {
      return false;
    }

    void leaf(const int n0, const int n1) const {
      const int p = points.prims(n0)[0]; 
      const auto tri = elements[faces.prims(n1)[0]];
      if (!tri.contains(p)) {
        const auto T = Tri(X[tri.x],X[tri.y],X[tri.z]);
        visit(p,tri,X[p],T);
      }
    }
  };

  template<class EE,class PF> void traverse(RawArray<const TV3> X, const EE& edge_edge, const PF& point_face) const {
    GEODE_ASSERT(X.size()==faces->mesh->nodes());
    this->X.copy(X);
    edges->update();
    double_traverse(*edges,EEVisitor<EE>(edges,X,edge_edge),close);
    if (include_faces) {
      points->update();
      faces->update();
      double_traverse(*points,*faces,PFVisitor<PF>(points,faces,X,point_face),close);
    }
  }

  T closest(RawArray<const TV3> X) const {
    T closest = inf;
    traverse(X,
      [&](const IV2 ei, const IV2 ej, const Seg si, const Seg sj) {
        closest = min(closest,segment_segment_distance(si,sj));
      },
      [&](const int p, const IV3 tri, const TV3 Xp, const Tri T) {
        closest = min(closest,magnitude(Xp-T.closest_point(Xp).x));
      });
    return closest;
  }

  T energy(RawArray<const TV3> X) const {
    T sum = 0;
    traverse(X,
      [&](const IV2 ei, const IV2 ej, const Seg si, const Seg sj) {
        sum += distance_energy(segment_segment_distance(si,sj));
      },
      [&](const int p, const IV3 tri, const TV3 Xp, const Tri T) {
        sum += distance_energy(magnitude(Xp-T.closest_point(Xp).x));
      });
    return sum;
  }

  Array<TV3> gradient(RawArray<const TV3> X) const {
    Array<TV3> grad(X.size());
    traverse(X,
      [&](const IV2 ei, const IV2 ej, const Seg si, const Seg sj) {
        const auto I = segment_segment_distance_and_normal(simplex(X[ei.x],X[ei.y]),simplex(X[ej.x],X[ej.y]));
        const T d = distance_gradient(I.x);
        grad[ei.x] -= d*(1-I.z.x)*I.y;
        grad[ei.y] -= d*   I.z.x *I.y;
        grad[ej.x] += d*(1-I.z.y)*I.y;
        grad[ej.y] += d*   I.z.y *I.y;
      },
      [&](const int p, const IV3 tri, const TV3 Xp, const Tri Xtri) {
        const auto I = Xtri.closest_point(Xp);
        auto N = Xp-I.x;
        const T d = distance_gradient(normalize(N));
        grad[p] += d*N;
        grad[tri.x] -= d*I.y.x*N;
        grad[tri.y] -= d*I.y.y*N;
        grad[tri.z] -= d*I.y.z*N;
      });
    return grad;
  }

  struct CollisionVisitor : public boost::noncopyable {
    const SimplexTree<TV3,1>& edges;
    const SimplexTree<TV3,2>& faces;
    RawArray<const Vector<int,2>> segs;
    RawArray<const Vector<int,3>> tris;
    RawArray<const TV3> X;
    int count;

    CollisionVisitor(const SimplexTree<TV3,1>& edges, const SimplexTree<TV3,2>& faces, RawArray<const TV3> X)
      : edges(edges), faces(faces), segs(edges.mesh->elements), tris(faces.mesh->elements), X(X), count(0) {}

    bool cull(const int n0, const int n1) const {
      return false;
    }

    void leaf(const int n0, const int n1) {
      const auto seg = segs[edges.prims(n0)[0]];
      const auto tri = tris[faces.prims(n1)[0]];
      if (!tri.contains(seg.x) && !tri.contains(seg.y)) {
        const Seg S(X[seg.x],X[seg.y]);
        const Tri T(X[tri.x],X[tri.y],X[tri.z]);
        Ray<TV3> ray(S);
        count += T.lazy_intersection(ray);
      }
    }
  };

  int collisions(RawArray<const TV3> X) const {
    GEODE_ASSERT(X.size()==faces->mesh->nodes());
    this->X.copy(X);
    edges->update();
    faces->update();
    CollisionVisitor visit(edges,faces,X);
    double_traverse(*edges,*faces,visit,close);
    return visit.count;
  }

  Array<TV3> strain_limit(RawArray<const T> restlengths, RawArray<const TV3> X, const T alpha) const {
    GEODE_ASSERT(restlengths.size()==edges->mesh->elements.size());
    const auto XL = X.copy();
    const Array<bool> frozen(X.size());
    for (const int s : range(restlengths.size())) {
      const auto e = edges->mesh->elements[s];
      TV3 v = X[e.y]-X[e.x];
      const T L = normalize(v);
      const TV3 dx = .5*alpha*(restlengths[s]-L)*v;
      XL[e.x] -= dx;
      XL[e.y] += dx;
    }
    traverse(X,
      [&](const IV2 ei, const IV2 ej, const Seg si, const Seg sj) {
        const auto I = segment_segment_distance_and_normal(simplex(X[ei.x],X[ei.y]),simplex(X[ej.x],X[ej.y]));
        if (I.x < close) {
          const TV3 dx = .5*alpha*(close-I.x)*I.y;
          XL[ei.x] -= (1-I.z.x)*dx;
          XL[ei.y] -=    I.z.x *dx;
          XL[ej.x] += (1-I.z.y)*dx;
          XL[ej.y] +=    I.z.y *dx;
        }
      },
      [&](const int p, const IV3 tri, const TV3 Xp, const Tri Xtri) {
        const auto I = Xtri.closest_point(Xp);
        auto N = Xp-I.x;
        const T d = normalize(N);
        if (d < close) {
          const TV3 dx = .5*alpha*(close-d)*N;
          XL[p] += dx;
          XL[tri.x] -= I.y.x*dx;
          XL[tri.y] -= I.y.y*dx;
          XL[tri.z] -= I.y.z*dx;
        }
      });
    return XL;
  }
};

GEODE_DEFINE_TYPE(SimpleCollisions)
}

GEODE_PYTHON_MODULE(fractal_helper) {
  GEODE_FUNCTION(boundary_edges_to_faces)
  GEODE_FUNCTION(iterate_L_system)
  GEODE_FUNCTION(branching_mesh)
  GEODE_FUNCTION(classify_loop_patches)
  GEODE_FUNCTION(min_instance_dot)
  GEODE_FUNCTION(settle_instances)
  GEODE_FUNCTION(torus_mesh)
  GEODE_FUNCTION(make_manifold)
  GEODE_FUNCTION(boundary_curve_at_height)
  GEODE_FUNCTION(unit_spring_energy_gradient)

  typedef SimpleCollisions Self;
  Class<Self>("SimpleCollisions")
    .GEODE_INIT(const TriangleSoup&,RawArray<const TV3>,const T,bool)
    .GEODE_METHOD(closest)
    .GEODE_METHOD(energy)
    .GEODE_METHOD(gradient)
    .GEODE_METHOD(collisions)
    .GEODE_METHOD(strain_limit)
    ;
}
