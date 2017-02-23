#!/usr/bin/env python

# Taken from Chris Gulvik's repository at 
#  https://github.com/chrisgulvik/genomics_scripts/blob/f6cb0d056f252714135ddb710903bc83bf7fbb8d/ANI.py

import os
import tempfile
import shutil
import re
import sys

from argparse import ArgumentParser
from collections import deque
from itertools import islice
from numpy import mean, std
from subprocess import Popen, call
from sys import exit
from time import strftime
from Bio import SeqIO

def parseArgs():
	parser = ArgumentParser(description='Computes the average nucleotide '
	'identity (ANI) between two nucleic acid sequence sets', add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('subject', nargs=1, metavar='subject.fasta',
		help='Subject FastA sequence file')
	req.add_argument('query', nargs='+', metavar='query.fasta...',
		help='Query FastA sequence file')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-c', '--cpus', type=str, metavar='INT',
		default='1', help='number of CPUs [1]')
	opt.add_argument('-f', '--fraction', type=float, metavar='FLOAT',
		default=70.0, help='minimum alignment length percentage [70.0]')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-i', '--identity', type=float, metavar='FLOAT',
		default=30.0, help='minimum percent identity [30.0]')
	opt.add_argument('-j', '--fill', type=str, metavar='CHAR', default='', 
		help='character to add to end of fragments that are less than the '
		'window size to force them to be the same length; requires the -k '
		'switch as well [None]')
	opt.add_argument('-k', '--keep-small-frags', default=False,
		action='store_true', help='keep remaining nucleotide sequences less '
		' than the window size')
	opt.add_argument('-l', '--length', type=int, metavar='INT',
		default=0, help='minimum alignment character length (sum of all '
		'aligned segments and all gaps) [0]')
	opt.add_argument('-o', '--outpath', metavar='PATH',
		default=None, help='output directory [./ANI--<date>_<time>]') # TODO default=stdout
	opt.add_argument('-s', '--step-size', type=int, metavar='INT',
		default=200, help='nucleotide step size during sequence '
		'fragmentation prior to BLAST alignments; to turn off steps to '
		'speed up, set to same length as the fragment size [200]')
	opt.add_argument('-w', '--fragment-size', type=int, metavar='INT',
		default=1000, help='fragment lengths to slice nucleotide sequence '
		'sets into prior to BLAST alignments [1000]')
	return parser.parse_args()

def fragment(seq, win, step, fill):
	iters = iter(seq)
	q = deque(islice(iters, win), maxlen=win)
	q.extend(fill for _ in range(win-len(q)))
	while True:
		yield q
		q.append(next(iters))
		q.extend(next(iters, fill) for _ in range(step-1))

def logmsg(msg):
  sys.stderr.write(msg + "\n")

def main(tempdir):
	opts = parseArgs()
	subject = os.path.abspath(opts.subject[0]);
	subjectName = os.path.basename(subject).rsplit('.')[0]

  # Check for dependencies
	if not dependency("blastn"):
		logmsg("ERROR: Not found: blastn")
		sys.exit(1)
	if not dependency("makeblastdb"):
		logmsg("ERROR: Not found: makeblastdb")
		sys.exit(1)

	# Temporary directory
	subjectFrags = sequenceToFrags(subject, tempdir, opts)

	# Run pairwise ANI for each query
	for query in opts.query:
		queryName = os.path.basename(query).rsplit('.')[0]
		queryFrags = sequenceToFrags(query, tempdir, opts)
		ani = ANI(subjectFrags, queryFrags, tempdir, opts)
		print "\t".join([subjectName, queryName, str(ani)]) + "\n"


def ANI(sFrags, qFrags, tmpdir, opts):
	# Execute a bidirectional blast
	subjectHits=blastFrags(sFrags, qFrags, tmpdir, opts)
	queryHits  = blastFrags(qFrags, sFrags, tmpdir, opts)

	aniTotal=0
	rbbHitCount=0 # Reciprocal best blast hit counter
	for subjectId in subjectHits.keys():
		queryId = subjectHits[subjectId]["sseqid"]
		# if the subject has a hit that exists in the queries
		if queryId in queryHits:
			# if the subjectId is the best hit of the query
			if subjectId == queryHits[queryId]["sseqid"]:
				aniTotal    += float(queryHits[queryId]["pident"])
				rbbHitCount += 1

	ani = round(aniTotal/rbbHitCount, 2)
	return ani


# Returns a dictionary of blast hits
def blastFrags(sFrags, qFrags, tmpdir, opts):
	blastoutOsFp, blastout= tempfile.mkstemp(dir=tmpdir, suffix=".tsv", prefix="blastout_")
	blastdb = sFrags + ".blastdb"
	c1 = ['makeblastdb', '-in', sFrags, '-out', blastdb, '-dbtype', 'nucl']
	c2 = ['blastn', '-db', blastdb, '-query', qFrags, '-dust', 'no',
		'-max_hsps', '1', '-max_target_seqs', '1',
		'-num_threads', opts.cpus, '-task', 'blastn',
		'-outfmt', '6 qseqid sseqid pident length bitscore qcovhsp',
		'-out', blastout]

	for cmd in [c1, c2]:
		return_code = Popen(cmd, stdout=sys.stderr)
		if return_code.wait() != 0:
			exit('ERROR: failed sys call\n{}'.format(' '.join(cmd)))

	# Load the results into memory
	result={}
	header=("qseqid", "sseqid", "pident", "length", "bitscore", "qcovhsp")
	numHeaders=len(header)
	with open(blastout) as subjectFh:
		for line in subjectFh:
			# populate the result dictionary
			F = re.split("\t",line)
			Fdict = {}
			for i in range(0,numHeaders):
				Fdict[header[i]]=F[i]
			result[Fdict["qseqid"]]=Fdict
	
	return result

def dependency(dep):
	''' checks for binary or script availability and
	returns the path if it is available '''
	for path in os.environ.get('PATH', '').split(':'):
		if os.path.exists(os.path.join(path, dep)) and \
			not os.path.isdir(os.path.join(path, dep)):
			return os.path.join(path, dep)

	return None

# Turn a fasta file into a file of sequence fragments
# Returns a fasta file with the fragments.
def sequenceToFrags(inseqFile,tmpdir,opts):
	fragsOsFh,fragsFile=tempfile.mkstemp(dir=tmpdir, suffix=".fasta", prefix="frags_")
	fragsFileObj=open(fragsFile, 'w')
	fragmentCount=0
	for rec in SeqIO.parse(inseqFile, 'fasta'):
		sequence=rec.seq
		seqLength=len(sequence)
		for i in range(0, seqLength-opts.fragment_size, opts.step_size):
			fragmentCount+=1
			subseq=sequence[i:i+opts.fragment_size]
			fragsFileObj.write('>{}\n{}\n'.format(fragmentCount, subseq))
	
	fragsFileObj.close()
	return fragsFile

def logmsg(msg):
	sys.stderr.write(msg + "\n")
			
if __name__ == '__main__':
	try:
		tempdir  = tempfile.mkdtemp(dir=tempfile.gettempdir(), prefix="ANItree_")
		logmsg("Temporary directory is " + tempdir)
		main(tempdir)
	finally:
		shutil.rmtree(tempdir)

