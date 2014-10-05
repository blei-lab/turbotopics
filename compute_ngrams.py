#! /usr/bin/python

# (C) Copyright 2009, David M. Blei (blei@cs.princeton.edu)

# This file is part of TURBOTOPICS.

# TURBOTOPICS is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your
# option) any later version.

# TURBOTOPICS is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA

import turbotopics as tt
import sys
from pprint import *

def main(corpus_file, pvalue, use_perm, out_filename, min_count=5,
         min_bigram_count=5, min_char_count=3):

    """recursively find collocations for a given corpus.  writes out
    the marginal counts to a specified file"""

    sys.stdout.write("computing n-grams from %s\n" % corpus_file)

    ### read corpus
    corpus = file(corpus_file).readlines()

    ### set up recursive hypothesis tests
    lr = tt.LikelihoodRatio(pvalue=pvalue, use_perm=use_perm)
    def iter_gen(): return(corpus)
    # note: some hidden defaults here, e.g., no numbers
    char_filter = tt.make_char_filter(min_char_count)
    def my_filter(w):
        char_filter(w) and tt.stop_filter(w) and tt.digit_filter(w)
    def update_fun(count, doc):
        count.update_counts(doc, root_filter=tt.stop_filter)

    ### compute significant n-grams
    cnts = tt.nested_sig_bigrams(iter_gen, update_fun, lr, min_count)

    ### write n-grams to file
    sys.stdout.write("writing to %s\n" % out_filename)
    f = file(out_filename, 'w')
    # this can be adjusted to write out any information you need
    [f.write('%s|%g\n' % (term, count)) for (term, count)
      in sorted(cnts.marg.items(), key=lambda x:-x[1])]
    f.close()

    return(cnts)


if __name__ == "__main__":

    from optparse import *

    parser = OptionParser()
    parser.add_option("--corpus", type="string", dest="corpus")
    parser.add_option("--perm",action="store_true", dest="use_perm")
    parser.add_option("--pval", type="float", dest="p_value")
    parser.add_option("--out", type="string", dest="out_filename")
    parser.add_option("--min-count", type="float", dest="min_count")
    parser.set_defaults(min_count=5)

    (opt, args) = parser.parse_args()

    main(corpus=opt.corpus, p_value=opt.p_value,
         use_perm_test=opt.use_perm, out_filename=opt.out_filename,
         min_count = opt.min_count)
