from __future__ import absolute_import

from .fractal_helper import *
from geode import *

def circle(r,n):
  "Delaunay triangulation of a circle with n sides"
  X = r*polar(2*pi/n*arange(n))
  return delaunay_points(X),X

def annulus(r0,r1,n0,n1):
  "Delaunay triangulation of an annulus with inner radius r0, outer radius r1, and counts n0, n1"
  X = concatenate([r0*polar(2*pi/n0*arange(n0)),r1*polar(2*pi/n1*arange(n1))])
  mesh = delaunay_points(X).mutate()
  for f in mesh.faces():
    if mesh.face_vertices(f).max()<n0:
      mesh.erase_face(f,0)
  return mesh,X

def variable_sor(radius,height,counts,periodic=False):
  "A closed surface of revolution where we allow a variable number of points at each height"
  n = len(radius)
  radius,height,counts = map(asarray,(radius,height,counts))
  assert radius.shape==height.shape==counts.shape==(n,)
  offsets = concatenate([[0],cumsum(counts)])
  def with_z(xy,z):
    X = empty((len(xy),3))
    X[:,:2] = xy
    X[:,2] = z
    return X
  X = [with_z(r*polar(2*pi/k*arange(k)),z) for r,z,k in zip(radius,height,counts)]
  tris = []
  if not periodic:
    tris.extend((delaunay_points(X[0][:,:2].copy()).elements()[:,::-1],
                 delaunay_points(X[-1][:,:2].copy()).elements()+offsets[-2]))
  for i in xrange(n-1+periodic):
    tris.append((annulus(1,2,counts[i],counts[(i+1)%n])[0].elements()[:,::-1]+offsets[i])%offsets[-1])
  return TriangleSoup(concatenate(tris)),concatenate(X)

def extrude((mesh,xy),z):
  "Extrude a flat mesh"
  xy = asarray(xy)
  z = asarray(z)
  assert 0<=z.ndim<=1
  if not z.ndim:
    z = array([-z,z])/2
  else:
    assert len(z)==2
  assert xy.ndim==2 and xy.shape[-1]==2
  n = len(xy)
  X = empty((2,n,3))
  X[:,:,:2] = xy
  X[:,:,2] = z[:,None]
  tris = []
  for e in mesh.boundary_edges():
    i = mesh.src(e)
    j = mesh.dst(e)
    tris.append((j,i,j+n))
    tris.append((i,i+n,j+n))
  tris = concatenate([tris,mesh.elements()[:,::-1],n+mesh.elements()])
  return TriangleSoup(tris.astype(int32)),X.reshape(-1,3)
