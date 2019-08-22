Generic Makefile for C, C++ and Fortran
=======================================

This is a generic makefile for compiling programs and libraries from Fortran,
C and C++ source files with minimal configuration.

Features
--------

* easily switch between compilers/optimisation levels/configurations using
  a diferent `ARCH` ``make.inc`` file (see below).
* compilation rules allow for both single- and mixed-language compilation.
* dependencies are automatically generated using the C or C++ compiler for
  C and C++ source code and using sfmakedepend (written by Kate Hedstrom) for
  Fortran code.
* produce executables and/or libraries.
* objects are placed in a configuration-specific directory and executables and
  libraries are given configuration-specific names to facilitate fast
  compilation using different settings (e.g. temporarily switching to a debug
  configuration and then back to an optimised configuration without having to
  recompile all source files each time).

This allows for rapid and easy comparison between different compile-time
settings.

File structure
--------------

``Makefile`` will create the directories as needed, with the following default
names (easily changed--see comments in ``Makefile``) relative to the working
directory:

bin/
    Contains compiled executables.
lib/
    Contains compiled libraries.
dest/depend/
    Contains files which provide the dependency rules governing the compilation
    order and which files need to be recompiled if a given file changes.
    Automatically generated and updated.
dest/CONFIG/OPT/
    Contains compiled object files.  CONFIG and OPT are variables set in
    a ``make.inc`` file.

.. warning::

    The dest directory is removed by the ``cleanall`` target and so should not
    contain **any** files which are not produced by a target in the makefile..
    Similarly **all** files in the dest/CONFIG/OPT/ directory are deleted by
    the ``clean`` target.

``Makefile`` can be used to create:

PROG_NAME.x
    A symbolic link to the most recently compiled PROG_NAME.CONFIG.OPT.x
    binary.  PROG_NAME is a variable set in ``Makefile``.
PROG_NAME.CONFIG.OPT.x
    Project executable (configuration-specific name).
libPROG_NAME.a
    A symbolic link to the most recently compiled PROG_NAME.CONFIG.OPT.a
    library.
libPROG_NAME.CONFIG.OPT.a
    Project library (configuration-specific name).

Usage
-----

#. Set the PROG_NAME and VPATH variables and (if desired or relevant) the MAIN
   and FORCE_REBUILD_FILES variables at the top of ``Makefile``.
#. Set the relevant variables in ``make.inc`` or ``make.inc.XXX`` (where XXX is
   an arbitrary variable).  Variables need not be set if they are not
   used---e.g.  the Fortran and C++ flags can be left unset for a pure
   C project.  An example ``make.inc`` file for the GCC compiler suite is
   included.
#. Run::

       $ make help

   to see available targets.  Running::

       $ make

   will compile the program.

   Using::

       $ make ARCH=XXX

   will use the settings in ``make.inc.XXX`` rather than ``make.inc``.

If Fortran is used, then sfmakedepend must be placed in a subdirectory entitled
``tools`` or the rule for generating $(F_DEPEND) must be modified accordingly.

Compilation
-----------

All files ending in ``.F``, ``.f``, ``.F90``, ``.f90``, ``.c`` and ``.cpp`` in
the directories in VPATH are automatically detected and compiled:

``.F``, ``.F90``
    passed through a pre-processor and then compiled using the Fortran compiler.
``.f``, ``.f90``
    compiled using the Fortran compiler.
``.c``
    passed through a pre-processor and then compiled using the C compiler.
``.cpp``
    passed through a pre-processor and then compiled using the C++ compiler.

Pre-processing is actually performed by the relevant compiler; this is standard
functionality for C and C++ compilers and very common for Fortran compilers.

Executables are produced by linking together all compiled object files using
the specified linker (i.e. ``ld`` or one of the compilers).

Libraries are produced by creating an archive (usually using ``ar``) from all
source files with the exception of (if relevant) the source file containing the
entry point to the program (i.e. the file containing the ``main`` procedure or
equivalent.

Compatibility
-------------

Only tested with GNU Make 3.81 and 3.82.  Unlikely to work with other versions
of Make and very unlikely to (fully) work with earlier versions of GNU Make.

mkconfig
--------

``tools/mkconfig`` can generate ``make.inc`` files from simple ini-style
configuration files.  An example configuration file for the GCC compiler suite,
``config/gnu``, is included.  See::

    tools/mkconfig --help

and::

    tools/mkconfig --help-long

for more details.


License
-------

MIT.  See comments in ``Makefile`` for more details. 

Acknowledgements
----------------

sfmakedepend (included for convenience in the tools subdirectory) was written
by Kate Hedstrom and is available under a MIT-based license.  See
http://www.myroms.org and tools/License_ROMS.txt for more details.
