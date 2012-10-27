Various fractal scripts
=======================

`other/fractal` contains scripts for a few different types of fractals.  The code started as a tutorial
for the Python programming language inside Otherlab.  In addition to very simple examples (L-system curves, 
noise curves, and a poker-based fractal), it contains code to generated so-called "developing fractal curves"
illustrating the refinement of L-system curves as a smooth surface:

  Geoffrey Irving and Henry Segerman, "Developing fractal curves" (in review).

The license is standard three-clause BSD (see the included `LICENSE` file or
[LICENSE](https://github.com/otherlab/fractal/blob/master/LICENSE)).

### Dependencies

The simple fractals (everything except for developing fractal curves), depend on

* [python >= 2.6](http://python.org): A scripting language
* [numpy >= 1.5](http://numpy.scipy.org): Efficient multidimensional arrays for Python
* [scipy](http://www.scipy.org): Scientific computation for Python
* [matplotlib](http://matplotlib.sourceforge.net): Python plotting

Developing fractal curves (`dragon.py` and `render-dragon`) additionally depend directly on

* [other/core](https://github.com/otherlab/core): Otherlab core library
* other/gui: Otherlab gui library
* [mitsuba >= 0.4.1](http://www.mitsuba-renderer.org): Physically based rendering

`other/core` has a few more indirect dependencies (boost, scons).

Unfortunately, other/gui is not yet open source, so the interactive features of `dragon.py`
will not work outside of Otherlab.  It can still be run in console mode, however, and then
visualized and rendered via Mitsuba.

### Setup

The simple scripts can be run immediately if their dependencies are available.  `dragon.py`
depends on a few C++ routines, which are built as follows:

1. Install `other/core` inside a root directory `other`.
2. Place `fractal` inside `other`.
3. Build

        cd other/fractal
        scons -u

### Usage

The simple scripts are run as follows:

    cd other/fractal
    ./basics
    ./poker
    ./laplace
    ./noise
    ./l-system # List available kinds
    ./l-system koch # Visualize a Koch snowflake

If `other/gui` is available, `dragon.py` can be run without arguments and all parameters
adjusted interactively.  Otherwise, it can run in console mode with `--console 1` to dump
out data for visualization in other programs.  There are two modes: instanced (the default)
and single mesh, controlled with `--instanced <flag>`.  Single mesh mode generates a single
smooth surface (suitable for 3D printing with an appropriate `--thickness` value).  Instanced
mode splits the surface into self-similar patches for efficient rendering in Mitsuba, with
the benefit that each patch isomophism class can be given a different pretty color.  See
`render-dragon` for the commands used in the "Developing fractal curves" paper.  Here is one
example:

    # Generate a level 13 dragon curve, writing mitsuba data to gen-dragon, and view it with Mitsuba
    ./dragon.py --type dragon --level 14 --smooth 4 --border-layers 4 --thickness 0.2 --ground 1 --rotation=0,-1,0,0 --mitsuba-dir gen-dragon --console 1
    ./render-dragon --gui 1 --view front --data gen-dragon

    # Generate a level 11 dragon curve, writing a single mesh to dragon.stl
    ./dragon.py --type dragon --level 12 --smooth 3 --border-crease 0 --thickness 0.8 --thickness-alpha .8 --instance 0 -o dragon.stl --console 1
