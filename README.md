JCF
===

__JCF__ is a companion library for __Jcup__, former Jcup/c.white.

JCF is a free software, you can use it under the terms of The MIT
License (MIT). See `LICENSE` file and User's guide for more details.
                                                              

JCF and Jcup
============

JCF is intended to be used to handle two or more ''grid'' or ''mesh''
to create mappingtable between them.  And more, for example, you can
use this to visualize the relation between grid points or control
volumes, etc.

When you create your own coupler with Jcup, you have to prepare the
''mappingtable'', which describe reation and weigts between grids of
coupled components by yourself.

JCF may help you to create this mappingtable.


JCF ''core'' and ''app''
========================

The `core` directory contains JCF core modules. The `app` directory
contains ''application'' using JCF core modules.

''Application'' modules should be implemented for each ''grid''
Jcup-used coupler handles. For exmaple, there is a NICAM-COCO coupled
model powered by Jcup, there is an application module for NICAM grid
and the one for COCO grid.

These modules are to be used not only mappingtable creation,
but also tools to visualize relation of two grids, etc.

Currently, application modules for NICAM grid, NICAM-IO grid are
implemented.  The `example` directory contains mappingtable-creation
program for these two grids, as a sample.

