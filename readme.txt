
        -------------------
        --- TURBOTOPICS ---
        -------------------

---------------------------------------------------------------------------

(C) Copyright 2009, David M. Blei (blei@cs.princeton.edu)

This file is part of TURBOTOPICS.

TURBOTOPICS is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

TURBOTOPICS is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 59
Temple Place, Suite 330, Boston, MA 02111-1307 USA

---------------------------------------------------------------------------


This code contains python scripts for running the TURBOTOPICS method
on either a corpus of documents or a corpus of documents combined with
the output of LDA-C.  for both scripts, the corpus is a file of the
original text of documents, one per line.

Note that neither script requires specifying how large N should be in
the N-grams.  For more information about the method, see the paper at

    http://arxiv.org/abs/0907.1013

The two scripts are

  compute_ngrams.py:

  Compute recursive multi-word expressions from a corpus.  This will
  write out a file of vocabulary (including multi-word expressions)
  and their counts.

  lda_topics.py:

  Compute multi-word expressions per-topic from a corpus and LDA-C
  fit.  (Note: the argument --ntopics is the same as K in LDA-C.)
  This will write out a file for each topic with the expressions and
  counts.  Again, see the paper for details.

Any questions/comments about this code should be posted to the topic
models mailing list.  Subscribe at

  https://lists.cs.princeton.edu/mailman/listinfo/topic-models
