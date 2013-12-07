Various fractal scripts
=======================

`fractal` contains scripts for a few different types of fractals, plus a bit of non-fractal
mathematical exploration.  The repo started as a tutorial for the Python programming language
inside Otherlab.  In addition to very simple examples (L-system curves, noise curves, and a
poker-based fractal), it contains code to generated so-called "developing fractal curves"
illustrating the refinement of L-system curves as a smooth surface:

  Geoffrey Irving and Henry Segerman, "Developing fractal curves" (in review).

The license is standard three-clause BSD (see the included `LICENSE` file or
[LICENSE](https://github.com/otherlab/fractal/blob/master/LICENSE)).

### Dependencies

The simple pure Python fractals depend on

* [python >= 2.6](http://python.org): A scripting language
* [numpy >= 1.5](http://numpy.scipy.org): Efficient multidimensional arrays for Python
* [scipy](http://www.scipy.org): Scientific computation for Python
* [matplotlib](http://matplotlib.sourceforge.net): Python plotting

Mixed Python/C++ code such as developing fractal curves (`dragon.py` and `render-dragon`)
additionally depend directly on

* [geode](https://github.com/otherlab/geode): Otherlab computational geometry library
* other/gui: Otherlab gui library
* [mitsuba >= 0.4.1](http://www.mitsuba-renderer.org): Physically based rendering

`geode` has a few more indirect dependencies (boost, scons).

Unfortunately, `other/gui` is not yet open source, so the interactive features of `dragon.py`
and other 3D interactive scripts will not work outside of Otherlab.  We will be releasing
and open source version soon, at which point `fractal` will be updated accordingly.  Note that
`dragon.py` can be run in console mode without `other/gui` and then visualized and rendered
via Mitsuba.

### Setup

The simple scripts can be run immediately if their dependencies are available.  For scripts
which use C++ (anything that fails on `import geode` or `import fractal_helper`), first
install `geode` via the instructions at https://github.com/otherlab/geode.  Then build the
C++ components of `fractal` via

    git clone https://github.com/otherlab/fractal.git
    cd fractal
    ln -s <path-to-geode>/config.py # If geode required configuration
    scons -j 5

If `geode` was built in place and not installed, add the following lines to `config.py` to
tell `fractal` where to find it:

    # In config.py
    geode_dir = '<path-to-geode>'
    geode_include = [geode_dir]
    geode_libpath = [geode_dir+'/build/$arch/$type/lib']

### Simple scripts

The simple scripts are run as follows:

    cd other/fractal
    ./basics
    ./poker
    ./laplace
    ./noise
    ./l-system # List available kinds
    ./l-system koch # Visualize a Koch snowflake

### Developing fractal curves

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
    ./dragon.py --type dragon --level 13 --smooth 4 --border-layers 4 --thickness 0.2 --ground 1 --rotation=0,-1,0,0 --mitsuba-dir gen-dragon --console 1
    ./render-dragon --gui 1 --view front --data gen-dragon

    # Generate a level 11 dragon curve, writing a single mesh to dragon.stl
    ./dragon.py --type dragon --level 11 --smooth 3 --border-crease 0 --thickness 0.8 --thickness-alpha .8 --instance 0 -o dragon.stl --console 1

### Hyperbolic tesselations

To see a triangular tiling of part of hyperbolic space represented in the Poincare disk model
(https://en.wikipedia.org/wiki/Poincare_disk), run

    ./poincare --help # For available options
    ./poincare
    ./poincare --degree 8 --level 4

To (attempt to) discretely isometrically embed the tiling in 3D Euclidean space, run

    ./poincare --mode flop <options>

This will use `other/gui` if available and fall back to matplotlib otherwise (in which case it will be
impossible to change options after startup).  To save the mesh to a file on startup, use

    ./poincare --mode flop --autosave poincare.obj
