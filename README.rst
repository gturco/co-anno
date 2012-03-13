CNS Pipeline
============
:Author: Gina Turco (`gturco <https://github.com/gturco>`_), Brent Pedersen (`brentp <http://github.com/brentp>`_)
:Email: gturco88@gmail.com
:License: MIT
.. contents ::


Description
===========
Coding sequences of one species are blast to the noncoding sequences of the other. Blastn is ran at a word size of 20 and E-Value < 0.001. Blast 
hits that hit the same coding region are summed by length. Groups with a sum greater then 100 are recorded as a missed exon strand.

.. image:: http://upload.wikimedia.org/wikipedia/commons/b/b4/Coanno.png

Installation
============
 + Python version >= 2.7

 + `blast <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/>`_
   (download latest and run)

 + `lastz <http://www.bx.psu.edu/~rsharris/lastz/newer/>`_
   (download latest .tar.gz; configure; make; make install) and adjust path in quota.sh)

 + `bpbio <http://pypi.python.org/pypi/processing#downloads>`_ (``svn checkout http://bpbio.googlecode.com/svn/trunk/ bpbio-read-only``)
   (run biostuff,coanno and bblast ``sudo python setup.py install``)


Run
===
Inputs
-------

 + Fasta File (it is recommended to run `50x mask repeat <http://code.google.com/p/bpbio/source/browse/trunk/scripts/mask_genome/mask_genome.py>`_)
 + Bed File (supports `UCSC bed format <http://genome.ucsc.edu/FAQ/FAQformat#format1>`_)
 + Converting GFF to Bed
  ``BCBio`` module required::
      
      python scripts/gff_to_bed.py rice_v6.gff >rice_v6.bed


 + If you have access to Coge the fasta and bed file for each organism can be obtained using export_to_bed.pl e.g.::

    perl scripts/export_to_bed.pl \
                          -fasta_name rice_v6.fasta \
                          -dsg 8163 \
                          -name_re "^Os\d\dg\d{5}$" > rice_v6.bed

   where ``dsg`` is from CoGe OrganismView and the prefix for the .bed and
   .fasta file **must be the same** (in this case ``rice_v6``).
   You likely need to run this on new synteny and then copy the .bed and
   .fasta files to the ``data/`` directory.
   The -name_re regular expression is not required, but in this case, it will
   prefer the readable Os01g101010 names over the names like m103430.


Editing Run File
::::::::::::::::

 + *Once only*: edit quota.sh to correct path for ``quota-alignment``
 + edit quota.sh to the correct `ORGA`, `ORGB`, `QUOTA`
 + run cmd: ``sh run.sh`` #that will call quota.sh (this will take a long time as it's doing a full blast (lastz) and then all of quota align, then cns pipeline).
 + this will create png's for the dotplots. check those to make sure the quota-blocks look correct.

Output files
::::::::::::

 + Query and subject CNS position
 + Missing Exons from ORGA ORGB blast
 + CNS blast to  RNA file
 + CNS blast to proteins file
 + CNS assigned to nearest Ortholog
