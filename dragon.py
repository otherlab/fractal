#!/usr/bin/env python

from __future__ import division
import os
import sys
import traceback
from other.core import *
from other.gui import *
from other.gui.show_tree import *
from other.core.value import parser
from other.core.openmesh import *
from other.fractal import *
from other.fractal import mitsuba

# Examples:
#
# 1. The version we printed: ./dragon.py --level 11 --smooth 3 --size 150 --thickness .7 --alpha .8 --z-scale .25
# 2. Bowl with missing base: ./dragon.py --level 10 --smooth 4 --size 150 --thickness .9 --closed 1 --closed-base 0 --alpha .8 --z-scale .25
# 3. Version for Eugene: ./dragon.py --level 14 --smooth 2 --size 150 --thickness .3 --z-scale .5

props = PropManager()
ftype = props.add('type','dragon').set_allowed('dragon terdragon koch gosper sierpinski'.split()).set_category('fractal')
levels = props.add('level',4).set_category('fractal')
scale_level = props.add('scale_level',-1).set_category('fractal')
smooth = props.add('smooth',0).set_category('fractal')
corner_shift = props.add('corner_shift',1/20).set_category('fractal')
size = props.add('size',150.).set_category('fractal')
thickness = props.add('thickness',0.).set_category('fractal')
thickness_alpha = props.add('thickness_alpha',-1.).set_category('fractal')
output = props.add('output','').set_category('fractal').set_abbrev('o')
closed = props.add('closed',False).set_category('fractal')
closed_base = props.add('closed_base',True).set_category('fractal')
z_scale = props.add('z_scale',.5).set_category('fractal')
sharp_corners = props.add('sharp_corners',False).set_category('fractal')
colorize = props.add('colorize',False).set_category('fractal')
color_seed = props.add('color_seed',184811).set_category('fractal')
instance = props.add('instance',True).set_category('fractal')
border_crease = props.add('border_crease',True).set_category('fractal')
border_layers = props.add('border_layers',1).set_category('fractal')
flip = props.add('flip',False).set_category('fractal')
curve_debug = props.add('curve_debug',-1).set_category('fractal')
two_ring = props.add('two_ring',False).set_category('fractal')
rearrange = props.add('rearrange',zeros(2)).set_category('fractal').set_hidden(1)

ground = props.add('ground',False).set_category('render')
min_dot_override = props.add('min_dot_override',inf).set_category('render')
settle_step = props.add('settle_step',.01).set_category('render')
mitsuba_dir = props.add('mitsuba_dir','').set_category('render')
origin = props.add('origin',(0,0,0)).set_category('render')
target = props.add('target',(0,0,0)).set_category('render')
rotation = props.add('rotation',Rotation.identity(3)).set_category('render')
console = props.add('console',False).set_category('render').set_help('skip gui')

extra_mesh_name = props.add('extra_mesh','').set_category('extra').set_help('draw an extra mesh for comparison')

@cache
def system():
  # start angle/level, shrink factor, axiom, rules, turns
  if ftype()=='dragon':
    return pi/4,sqrt(1/2),'fx',{'x':'x-yf','y':'fx+y'},{'+':pi/2,'-':-pi/2}
  elif ftype()=='terdragon':
    return pi/6,sqrt(1/3),'f',{'f':'f-f+f'},{'+':2*pi/3,'-':-2*pi/3}
  elif ftype()=='koch':
    return 0,1/3,'f+f+fc',{'f':'f-f+f-f'},{'+':2*pi/3,'-':-pi/3}
  elif ftype()=='gosper':
    a = angle(complex(sqrt(25/28),sqrt(3/28)))
    return -a,sqrt(1/7),'fx',{'x':'x+yf++yf-fx--fxfx-yf+','y':'-fx+yfyf++yf+fx--fx-y'},{'+':pi/3,'-':-pi/3}
  elif ftype()=='sierpinski':
    return [pi/3,-pi/3],1/2,'xf',{'x':'yf-xf-y','y':'xf+yf+x'},{'+':pi/3,'-':-pi/3}
  assert 0

def heights_helper(levels):
  heights = []
  z = 0
  alpha = system()[1]
  for level in xrange(levels+1):
    heights.append(z)
    z += z_scale()*alpha**level
  return asarray(heights)
heights = cache(lambda:heights_helper(levels()))

def attach_z(xy,z):
  X = empty((len(xy),3))
  X[:,:2] = xy
  X[:,2] = z
  return X

def curves_helper(levels):
  closed = is_closed()
  debug = curve_debug()
  with Log.scope('curves'):
    # Generate fractal
    start,shrink,axiom,rules,turns = system()
    curves = iterate_L_system(start,shrink,axiom,rules,turns,levels)
    if debug>=0:
      for level,curve in enumerate(curves):
        dx = curve[-1]-curve[0]
        Log.write('level %d, dist %g, dx %s, angle %s'%(level,magnitude(dx),repr(dx),repr(vector.angle(dx))))
    # Flip if desired
    if flip():
      for curve in curves:
        curve[:,1] *= -1
    # Shift corners inwards 
    for level,curve in enumerate(curves):
      Log.write('level %d, vertices %d'%(level,len(curve)))
      full = concatenate([[curve[-1]],curve,[curve[0]]]) if closed else curve
      partial = curve if closed else curve[1:-1]
      partial[:] += corner_shift()*(full[2:]+full[:-2]-2*full[1:-1])
    if debug>=0:
      import pylab
      curve = curves[debug]
      pylab.plot(curve[:,0],curve[:,1])
      pylab.show()
    return curves
curves = cache(lambda:curves_helper(levels()))

@cache
def is_closed():
  return system()[2][-1]=='c'

@cache
def branching():
  assert levels()>0
  open = not is_closed()
  a,b = (len(curves()[1])-open),(len(curves()[0])-open)
  branching = a//b
  assert a==branching*b
  Log.write('branching = %d'%branching)
  return branching

@cache
def mesh():
  curves_ = curves()
  shrink = system()[1]
  closed = is_closed()
  base = len(curves_[0])
  with Log.scope('mesh'):
    # Assemble mesh
    mesh = branching_mesh(branching(),levels(),base,closed)
    X = concatenate([attach_z(curve,height) for curve,height in zip(curves_,heights())])
    assert mesh.nodes()==len(X)
    assert not len(mesh.nonmanifold_nodes(True))
    # Rescale
    if scale_level()>=0:
      sX = concatenate([attach_z(curve,height) for curve,height in zip(curves_helper(scale_level()),heights_helper(scale_level()))])
    else:
      sX = X
    Xmin = sX.min(axis=0)
    Xmax = sX.max(axis=0)
    sizes = Xmax-Xmin
    scale = size()/sizes.max()
    Log.write('scale = %g'%scale)
    Log.write('sizes = %g %g %g'%(tuple(scale*sizes)))
    center = .5*(Xmin+Xmax)
    X = scale*(X-center)
    # Compute thickness field
    alpha = thickness_alpha()
    if alpha<0:
      alpha = shrink
    thick = hstack([repeat(alpha**level,(base+closed-1)*branching()**level+1-closed) for level in xrange(levels()+1)])
    thick *= thickness()/thick[-1]
    if scale_level()>=0:
      thick *= alpha**(levels()-scale_level())
    # Label patches
    patch = arange(len(mesh.elements)//(1+branching())).repeat(1+branching())
    return mesh,X,thick,patch,(scale,center)

@cache
def instances():
  m,X,thick,_,_ = mesh()
  with Log.scope('classify'):
    _,interior,boundary = classify_loop_patches(m,X,thick,1+branching(),two_ring())
    if any(rearrange()):
      i,b = map(int,rearrange())
      interior = interior[:i]+boundary[:b]+interior[i:]+boundary[b:]
      boundary = []
    return interior,boundary

def closed_mesh():
  assert not instance()
  m,X,_,patch,_ = mesh()
  # Four rotated and translated copies of the dragon curve form a closed curve
  translations = (0,0),(1,0),(1,1),(0,1)
  scale = magnitude(X[0]-X[1])
  cX = vstack([Frames(scale*hstack([t,0]),Rotation.from_angle_axis(pi/2*i,(0,0,1)))*X for i,t in enumerate(translations)])
  n = len(X)
  tris = vstack([m.elements+n*i for i in xrange(4)]+[asarray([(1,0,2*n),(0,2*n+1,2*n)],dtype=int32)]*closed_base())
  patch = concatenate([patch]*4+[[patch[-1]+1]]*2)
  # Weld meshes together
  p = ParticleTree(cX,10).remove_duplicates(1e-8)
  ip = empty(p.max()+1,int32)
  ip[p] = arange(len(p),dtype=int32)
  cX = cX[ip]
  tris = p[tris]
  return TriangleMesh(tris),cX,patch

def smoothed_mesh(mesh,X,thick,sharp):
  for _ in xrange(smooth()):
    sub = TriangleSubdivision(mesh)
    if sharp:
      sub.corners = array([0,1],dtype=int32)
    mesh = sub.fine_mesh
    X = sub.loop_subdivide(X)
    if thick is not None:
      thick = sub.loop_subdivide(thick)
  return mesh,X,thick

@cache
def smoothed():
  m,X,thick,patch = closed_mesh() if closed() else mesh()[:4] 
  m,X,thick = smoothed_mesh(m,X,thick,sharp=sharp_corners())
  patch = patch.repeat(4**smooth(),axis=0)
  Log.write('smoothed: triangles = %d, vertices = %d'%(len(m.elements),len(X)))
  return m,X,thick,patch

@cache
def smoothed_instances():
  def smooth_instance((mesh,X,thick,frames),sharp):
    sm,sX,st = smoothed_mesh(mesh,X,thick,sharp)
    norm = sm.vertex_normals(sX)
    cut = TriangleMesh(sm.elements[:(1+branching())*4**smooth()])
    # Rearrange so that all interior vertices come first
    interior = unique(cut.elements)
    is_interior = repeat(False,len(sX))
    is_interior[interior] = True
    vmap = empty(len(sX),dtype=int32)
    vmap[interior] = arange(len(interior),dtype=int32)
    vmap[logical_not(is_interior)] = arange(len(interior),len(sX),dtype=int32)
    inv = vmap.copy()
    inv[vmap] = arange(len(vmap),dtype=int32)
    # Apply mapping
    cut = TriangleMesh(vmap[cut.elements])
    assert cut.nodes()==len(interior)
    sm = TriangleMesh(vmap[sm.elements])
    return cut,sm,sX[inv],st[inv],norm[inv],frames
  interior,boundary = instances()
  with Log.scope('smoothing'):
    interior = [smooth_instance(inst,False) for inst in interior]
    boundary = [smooth_instance(inst,sharp_corners() and i==0) for i,inst in enumerate(boundary)]
  return interior,boundary

def thicken_mesh(mesh,X,thick,normals,border=None):
  layers = border_layers()
  n = len(X)
  m = n
  offset = .5*thick[...,None]*normals
  X = vstack([X+offset,X-offset])
  if border is None:
    border = mesh.boundary_mesh().elements
  border_nodes = unique(border)
  if border_crease():
    m = len(border_nodes)
    inv_border_nodes = empty(len(X),dtype=int32)
    inv_border_nodes[border_nodes] = arange(len(border_nodes),dtype=int32)
    X = vstack([X]+[(1-k/layers)*X[border_nodes]+k/layers*X[border_nodes+n] for k in xrange(layers+1)])
    border = inv_border_nodes[border]+2*n
    border = [border+k*m for k in xrange(layers+1)]
    border_nodes = arange(2*n,len(X),dtype=int32)
  else:
    assert layers==1
    border = [border,m+border]
    border_nodes = hstack([border_nodes,border_nodes+n])
  tris = [mesh.elements,mesh.elements[:,::-1]+n]
  for k in xrange(layers):
    tris.append(vstack([border[k].T[::-1],border[k+1][:,1]]).T)
    tris.append(vstack([border[k+1].T,border[k][:,0]]).T)
  mesh = TriangleMesh(ascontiguousarray(vstack(tris)))
  if 0:
    assert not len(mesh.nonmanifold_nodes(border_crease()))
  return mesh,X,border_nodes

@cache
def thicken():
  flat,X,thick,patch = smoothed()
  if not thickness():
    assert not len(flat.nonmanifold_nodes(True))
    return flat,X,patch
  mesh,X,_ = thicken_mesh(flat,X,thick,flat.vertex_normals(X))
  if patch is not None:
    border_face = boundary_edges_to_faces(flat,flat.boundary_mesh().elements)
    border_patch = patch[border_face]
    patch = concatenate([patch]*2+[border_patch]*(border_layers()+1))
  return mesh,X,patch

def write_mesh(filename,mesh,X,normals=None):
  trimesh = TriMesh()
  trimesh.add_vertices(X)
  trimesh.add_faces(mesh.elements)
  if normals is None:
    trimesh.write(filename)
  else:
    trimesh.set_vertex_normals(normals)
    trimesh.write_with_normals(filename)

@cache
def thicken_instances():
  layers = border_layers()
  def thicken_instance((mesh,full_mesh,X,thick,normals,frames)):
    if thickness():
      # Thicken representative with neighbors
      full_nodes = full_mesh.nodes()
      full_mesh,full_X,full_border = thicken_mesh(full_mesh,X,thick,normals)
      full_normals = full_mesh.vertex_normals(full_X)
      # Thicken representative without neighbors
      nodes = mesh.nodes()
      border = full_mesh.boundary_mesh().elements
      border = border[all(border<nodes,axis=-1)]
      mesh,X,border = thicken_mesh(mesh,X[:nodes],thick[:nodes],normals[:nodes],border=border)
      # Combine full normals with trimmed mesh
      normals = zeros((len(X),3))
      normals[:nodes] = full_normals[:nodes]
      normals[nodes:][:nodes] = full_normals[full_nodes:][:nodes]
      b = len(border)//(layers+1)
      fb = len(full_border)//(layers+1)
      for k in xrange(layers+1):
        normals[border[b*k:b*k+b]] = full_normals[full_border[fb*k:fb*k+b]]
      assert allclose(sqr_magnitudes(normals),1)
      # All done
      return mesh,X,normals,frames
    else:
      nodes = mesh.nodes()
      return mesh,X[:nodes],normals[:nodes],frames
  smoothed = smoothed_instances()
  with Log.scope('thicken'):
    return [map(thicken_instance,instances) for instances in smoothed]

@cache
def instance_lists():
  lists = []
  interior,boundary = thicken_instances()
  for mesh,X,normals,frames in interior+boundary:
    lists.append(cache_render(lambda:gl_triangles(X[mesh.elements]))())
  return lists

@cache_render
def render():
  try:
    if not instance():
      mesh,X,patch = thicken()
      if colorize():
        m = patch.max()+1
        random.seed(color_seed())
        color = wheel_color(random.uniform(0,1,size=m))[patch]
        gl_colored_triangles(color,X[mesh.elements])
      else:
        GL.glColor(1,0,0)
        gl_triangles(X[mesh.elements])
    else:
      interior,boundary = thicken_instances()
      lists = instance_lists()
      random.seed(color_seed())
      for (mesh,X,normals,frames),list in zip(interior+boundary,instance_lists()):
        GL.glColor(*(wheel_color(random.uniform(0,1)) if colorize() else (1,0,0)))
        gl_instances(frames,list)
  except:
    traceback.print_exc()

@cache
def ground_mesh():
  X = 10*size()*asarray([[-1,-1,0],[1,-1,0],[1,1,0],[-1,1,0]])
  mesh = TriangleMesh([[0,1,2],[0,2,3]])
  return mesh,X

@cache
def ground_frame():
  if not ground():
    return Frame.identity(3)
  up = up_frame().r.inverse()*(0,0,1)
  min_dot = min_dot_override()
  if not isfinite(min_dot):
    if instance():
      min_dot = inf
      for instances in thicken_instances():
        for _,X,_,frames in instances:
          min_dot = min(min_dot,min_instance_dot(X,frames,up))
    else:
      _,X,_,_,_ = mesh()
      min_dot = dots(X,up).min()
  Log.write('min dot = %g'%min_dot)
  return Frames(up*min_dot,up_frame().r.inverse())

def settle():
  interior,boundary = thicken_instances()
  instances = [(mesh,X,frames) for mesh,X,normals,frames in interior+boundary]
  with Log.scope('settle'):
    up = up_frame().r.inverse()*(0,0,1)
    new_up = settle_instances(instances,up,settle_step())
    up_frame.set(Frames(zeros(3),Rotation.from_rotated_vector(up,new_up)*up_frame().r.inverse()).inverse())

class DragonScene(Scene):
  def render(self,*args):
    if instance():
      instance_lists()
    with gl_scope():
      if ground():
        mesh,X = ground_mesh()
        GL.glColor(array([1,.775,.5431]))
        gl_triangles(X[mesh.elements])
        GL.glMultMatrixf(ground_frame().inverse().matrix().T)
      render().call()

  @staticmethod
  @cache
  def bounding_box():
    _,X,_,_,_ = mesh()
    return Box(X.min(axis=0),X.max(axis=0))

def raw_extra_mesh():
  name = extra_mesh_name()
  assert name
  tm = TriMesh()
  tm.read(name) 
  return TriangleMesh(tm.elements()),tm.X()

def extra_mesh():
  # If true, use known values for the transforms
  hack = True

  # Look up extra mesh information
  emesh,eX = raw_extra_mesh()
  if not hack:
    ebase = boundary_curve_at_height(emesh,eX,eX[:,2].min())
  eX = eX.copy()
  eX[:,2] *= -1 # Reflect
  if not hack:
    # Trim off extra bits on either end
    n = len(ebase)
    assert popcount(n)==3
    trim = min_bit(n-min_bit(n))//2
    ebase = ebase[trim:-trim]
    Log.write('ebase = %d'%len(ebase.shape))
    ebase_scale = magnitude(eX[ebase[-1]]-eX[ebase[0]])
    Log.write('ebase scale = %r'%ebase_scale)
  else:
    ebase_scale = 72.002998352050781

  # Look up information about Loop version
  if hack:
    frame = Frame.from_reals((-32.473187150735164,-13.043779937542258,79.059441919669837,0.99999997105133731,0,0,-0.00024061862871912305))
    base_scale = 77.689804371717031
    base_height = 150
  else:
    base = curves()[-1]
    Log.write('base = %s, height = %g'%(base.shape,heights()[-1]))
    _,_,_,_,(scale,center) = mesh()
    base_height = scale*abs(heights()[0]-heights()[-1])
    base = scale*(attach_z(base,heights()[-1])-center)
    base_scale = magnitude(base[-1]-base[0])
    Log.write('extra base scale = %r'%base_scale)
    Log.write('extra base height = %r'%base_height)

  # Fix length
  eX *= base_scale/ebase_scale

  # Fix height
  ratio = base_height/(eX[:,2].max()-eX[:,2].min())
  Log.write('extra height ratio = %r'%ratio)
  eX[:,2] *= ratio

  # Compute frame
  if not hack:
    frame = rigid_register(eX[ebase],base[::-1].copy())
    Log.write('extra frame = %r'%(tuple(frame.reals()),))
  if ground():
    frame = ground_frame().inverse()*frame
  return emesh,frame*eX

def save_mesh():
  try:
    mesh,X,_ = thicken()
    tm = TriMesh()
    tm.add_vertices(X)
    tm.add_faces(mesh.elements)
    with Log.scope('mesh write'):
      Log.write('filename = %s'%output())
      tm.write(output())
  except:
    traceback.print_exc()

def save_mitsuba():
  try:
    if instance():
      interior,boundary = thicken_instances()
    else:
      mesh,X,patch = thicken()
    with Log.scope('mitsuba write'): 
      dir = mitsuba_dir()
      Log.write('dir = %s'%dir)
      assert dir
      if not os.path.exists(dir):
        os.makedirs(dir)
      if instance():
        open(os.path.join(dir,'count'),'w').write('%d\n'%len(interior+boundary))
        for i,(mesh,X,normals,frames) in enumerate(interior+boundary):
          filename = 'rep-%d.obj'%i
          write_mesh(os.path.join(dir,filename),mesh,X,normals)
          open(os.path.join(dir,'instances-%d.xml'%i),'w').write(
            mitsuba.scene_version('0.4.1',
              mitsuba.instances('group-%d'%i,ground_frame().inverse().matrix()*frames)))
      else:
        open(os.path.join(dir,'count'),'w').write('%d\n'%0)
        filename = 'single.obj'
        write_mesh(os.path.join(dir,filename),mesh,ground_frame().inverse()*X)
      if extra_mesh_name():
        write_mesh(os.path.join(dir,'extra.obj'),*extra_mesh())
  except:
    traceback.print_exc()

def save_misc():
  try:
    with Log.scope('misc'):
      cam = main.view.cam
      f = cam.frame
      origin.set(f*(0,0,cam.distance))
      Log.write('old up = %s'%(cam.frame.r*(0,1,0)))
      target.set(f.t)
      rotation.set(up_frame().r)
      Log.write('scene information:')
      Log.write('  command = %s'%parser.command(props))
      Log.write('  command all = %s'%parser.command(props,0))
  except:
    traceback.print_exc()

def load_misc(props_set):
  up_frame.set(Frames(zeros(3),rotation()))
  if not console():
    if any([s in props_set for s in 'origin target rotation'.split()]):
      cam = main.view.cam
      d = target()-origin()
      cam.distance = magnitude(d)
      d = normalized(d)
      up = normalized(projected_orthogonal_to_unit_direction((0,0,1),d))
      cam.frame = Frames(target(),Rotation.from_matrix(vstack([cross(d,up),up,-d]).T.copy()))
    else:
      main.view.show_all(True)

def dump():
  thicken_instances.dump(0)

if __name__=='__main__':
  Log.configure('fractal',0,0,100)
  props_set = parser.parse(props,'Dragon curve visualizer')
  if console():
    up_frame = Prop('up_frame',Frame.identity(3))
    load_misc(props_set)
    if output():
      save_mesh()
    if mitsuba_dir():
      save_mitsuba()
  else:
    from OpenGL import GL
    app = QEApp(sys.argv,True)
    main = MainWindow(props)
    main.view.add_scene('dragon',DragonScene())
    extra = MeshScene(props,cache(lambda:extra_mesh()[0]),cache(lambda:extra_mesh()[1]),(.1,.3,.9),(.1,.9,.2))
    extra.active = cache(lambda:bool(extra_mesh_name()))
    main.view.add_scene('extra',extra)
    up_mode = main.view.add_interaction_mode(PlaneInteractionMode('up'),False,'Ctrl+u')
    up_frame = up_mode.frame
    main.add_menu_item('File','Save mesh',save_mesh,'Ctrl+s')
    main.add_menu_item('File','Save camera and command',save_misc,'Ctrl+c')
    main.add_menu_item('File','Save mitsuba data',save_mitsuba,'Ctrl+m')
    main.add_menu_item('Edit','Settle (fall down)',settle,'Ctrl+f')
    main.add_menu_item('Edit','Dump dependencies',dump,'Ctrl+d')
    main.add_menu_item('Edit','Show dependencies',curry(show_tree,main,True),'Ctrl+t')
    load_misc(props_set)
    main.init()
    app.run()
