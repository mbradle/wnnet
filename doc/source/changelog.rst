Changelog
=========

All notable changes to this project will be documented in this file.  This
project adheres to `Semantic Versioning <http://semver.org/spec/v2.0.0.html>`_.

Version 2.1.2
-------------

New:

  * Project now uses pyproject configuration file.

Version 2.1.1
-------------

Fix:

  * Some documentation formatting has been fixed.

Version 2.1.0
-------------

New:

  * Importation of the gslconsts math constants has been added.

Fix:

  * Documentation on the importation of the gslconsts namespaces has been
    clarified.

Version 2.0.3
-------------

Fix:

  * The selection of nuclides and reactions has been changed to be more
    consistent with expected behavior.

Version 2.0.2
-------------

Fix:

  * The Zenodo badge has been fixed.

Version 2.0.1
-------------

Fix:

  * An edge key has been fixed.

Internal:

  * The Zenodo doi has been updated to resolve to most recent version.

Version 2.0.0
-------------

New:

  * A number of method and function prototypes have been altered.  These are
    backwards-incompatible changes.
  * Linting and testing have been added.
  * It is now possible to add special arcs to graphs.

Internal:

  * A number of internal changes have been made to improve code operation and
    readability.

Version 1.2.6
-------------

Fix:

  * A deprecated math function has been updated.
  * Initialization of an array for graph anchors has been fixed.

Version 1.2.5
-------------

Fix:

  * The graph node positioning has been changed for better rendering of graphs.

Version 1.2.4
-------------

Fix:

  * The Q value for beta+ decay has been fixed.

Version 1.2.3
-------------

Fix:

  * Missing sphinx themes have been added for proper documentation building.

Version 1.2.2
-------------

New:

  * A configuration file has been added for proper documentation building
    on readthedocs.

Version 1.2.1
-------------

Fix:

  * An error in selecting both reaction directions on flows has been fixed.

Version 1.2.0
-------------

New:

  * The reaction validity check now ensures baryon number, charge,
    and lepton number conservation.

Version 1.1.0
-------------

New:

  * The ability to add user-rate functions has been added to the API.

Version 1.0.3
-------------

Fix:

  * An error in negative values for current graphs has been fixed.
  * Some tutorial text has been corrected.

Internal:

  * An execution of black has been added to the build script.

Version 1.0.2
-------------

Fix:

  * The tutorial notebook installation of graphviz and libgraphviz-dev has been     fixed.
  * Some tutorial text has been updated and corrected.

Version 1.0.1
-------------

Fix:

  * The tutorial notebook installation of pygraphviz has been fixed.

Version 1.0.0
-------------

New:

  * Initial release.

