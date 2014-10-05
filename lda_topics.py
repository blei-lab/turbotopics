#! /usr/bin/python

# -------------------------------------------------------------------------
# This code is for computing significant n-grams within topics that
# are obtained with LDA-C.  It can be used as a module or as a script.
# -------------------------------------------------------------------------

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

# create turbo topics from the output of LDA C

import turbotopics as tt, itertools, sys, re

def read_vocab(vocab_fname):

    "reads a vocabulary and returns a map of words to indices"

    sys.stdout.write('reading vocabulary from %s\n' % vocab_fname)
    terms = file(vocab_fname).read().split('\n')
    return(terms)


def parse_word_assignments(assigns_fname, vocab):
    """
    given a word assignments file and a list of words,
    returns a list of dictionaries mapping words to topics
    """
    results = []
    for assign in file(assigns_fname):
        wordmap = {}
        for (term, topic) in [x.split(':') for x in assign.split(' ')[1:]]:
            wordmap[vocab[int(term)]] = int(topic)
        results.append(wordmap)
    return(results)


def update_counts_from_topic(doc, topicmap, topic, counts_obj):
    """
    updates the counts of a counts object from a
    - doc : line of text
    - topicmap : mapping of words to topics
    - topic : integer of the topic to focus on
    - counts_obj : counts object to update
    """
    # the word filter looks if the first word is assigned to the topic
    if (topic not in topicmap.values()): return
    root_filter = lambda w: topicmap.get(w.split()[0], -1)==topic
    counts_obj.update_counts(doc, root_filter=root_filter)


def turbo_topic(corpus, assigns, topic, use_perm=False, pvalue=0.1, min=25):

    def iter_gen(): return(itertools.izip(corpus, assigns))
    def update_fun(counts, doc):
        update_counts_from_topic(doc[0], doc[1], topic, counts)

    test = tt.LikelihoodRatio(pvalue, use_perm=use_perm)
    cnts = tt.nested_sig_bigrams(iter_gen, update_fun, test, min)

    return(cnts)


if (__name__ == "__main__"):

    from optparse import *

    parser = OptionParser()
    parser.add_option("--corpus", type="string", dest="corpus")
    parser.add_option("--assign", type="string", dest="assignments")
    parser.add_option("--vocab", type="string", dest="vocab")
    parser.add_option("--perm",action="store_true", dest="use_perm")
    parser.add_option("--pval", type="float", dest="pvalue")
    parser.add_option("--out", type="string", dest="out")
    parser.add_option("--min-count", type="float", dest="min_count")
    parser.add_option("--ntopics", type="int", dest="ntopics")
    parser.set_defaults(min_count=25, use_perm=False, pval=0.001)

    (opt, args) = parser.parse_args()

    vocab = read_vocab(opt.vocab)
    assigns = parse_word_assignments(opt.assignments, vocab)
    corpus = file(opt.corpus).readlines()

    for topic in range(opt.ntopics):
        sys.stdout.write('writing topic %d\n' % topic)
        sig_bigrams = turbo_topic(corpus, assigns, topic,
                                  use_perm=opt.use_perm,
                                  min=opt.min_count,
                                  pvalue=opt.pvalue)
        tt.write_vocab(sig_bigrams.marg,
                       '%stopic%03d.txt' % (opt.out, topic))


# python lda_topics.py --assign=word-assignments.dat --corpus=corpus.txt --vocab=vocab.dat --out=tt --pval=0.001 --min-count=25 --perm
