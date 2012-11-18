#!/usr/bin/env python

from __future__ import division
from other.core import *
from numpy import *

def indent(lines):
  return ['  '+line for line in lines]

def header():
  return ['<?xml version="1.0" encoding="utf-8"?>','<!-- Automatically generated: do not edit! -->']

def concat(lists):
  s = []
  for x in lists:
    s += x
  return s

def tag(tagname,children=[],**fields):
  s = '<%s '%tagname + ' '.join('%s="%s"'%(k,v) for k,v in fields.iteritems())
  if not children:
    return [s+'/>']
  else:
    return [s+'>']+indent(concat(children))+['</%s>'%tagname]

def scene_version(version,*children):
  return '\n'.join(header()+tag('scene',version=version,children=children))

def scene(*children):
  return scene_version('0.3.0',*children)

field_types = {int:'integer',float:'float',str:'string',bool:'boolean'}
def fields(**fields):
  def convert(v):
    if type(v)==bool:
      return str(v).lower()
    else:
      return v
  return concat(tag(field_types[type(v)],name=n,value=convert(v)) for n,v in fields.iteritems())

def typed(name):
  def f(type,*children):
    return tag(name,type=type,children=children)
  return f
integrator = typed('integrator')
sampler = typed('sampler')
film = typed('film')
rfilter = typed('rfilter')
shape = typed('shape')
luminaire = typed('luminaire')
emitter = typed('emitter')
camera = typed('camera')
sensor = typed('sensor')
phase = typed('phase')
rfilter = typed('rfilter')

def transform(name,*children):
  return tag('transform',name='toWorld',children=children)

def scale(**fields):
  return tag('scale',**fields)

def translate(v):
  x,y,z = v
  return tag('translate',x=x,y=y,z=z)

def number(name,v):
  return tag('float',name=name,value=v)

def vector(name,v):
  x,y,z = v
  return tag('vector',name=name,x=x,y=y,z=z)

def point(name,v):
  x,y,z = v
  return tag('point',name=name,x=x,y=y,z=z)

def rotate(r):
  angle,axis = r.angle_axis()
  return tag('rotate',x=axis[0],y=axis[1],z=axis[2],angle=180/pi*angle)

def matrix(m):
  assert m.shape==(4,4)
  return tag('matrix',value=' '.join(map(str,m.ravel())))

def comma_sep(v):
  return ','.join(map(str,v))

def translate(vec):
  return tag('translate', x = vec[0], y = vec[1], z = vec[2])

def scale(vec):
  return tag('scale', x = vec[0], y = vec[1], z = vec[2])

def lookAt(origin,target,up):
  return tag('lookAt',origin=comma_sep(origin),target=comma_sep(target),up=comma_sep(up))

def bsdf(id,type,*children):
  return tag('bsdf',id=id,type=type,children=children)

def medium(id,type,*children):
  return tag('medium',id=id,name=id,type=type,children=children)

def spectrum(name,value):
  return tag('spectrum',name=name,value=value)

def blackbody(name, temp):
  return tag('blackbody', name=name, temperature=temp)

def ref(id, name = None):
  if name is None:
    return tag('ref', id=id)
  else:
    return tag('ref', name=name, id=id)

def shapegroup(id,*children):
  return tag('shape',type='shapegroup',id=id,children=children)

def instance(id,frame):
  if isinstance(frame,Frames):
    trans = transform('toWorld',rotate(frame.r),translate(frame.t))
  else:
    trans = transform('toWorld',matrix(frame))
  return shape('instance',
    ref(id),
    trans)

def instances(id,frames):
  return concat(instance(id,f) for f in frames)

def rgb(name,color):
  return tag('rgb',name=name,value=comma_sep(color))

def skip(*args,**kwargs):
  return []

def skipif(criterion, *children):
  if criterion:
    return []
  else:
    return concat(children)

def select(criterion, list1, list2):
  if criterion:
    return list1
  else:
    return list2

def include(filename):
  return tag('include',filename=filename)
