#!/usr/bin/env python
from __future__ import with_statement

from sys import stdin, stderr
from optparse import OptionParser
from vfork.io.util import safe_rstrip
from vfork.util import exit, format_usage
import re
from math import log
#from random import shuffle
#from subprocess import Popen, PIPE
#from collections import defaultdict
#from vfork.io.colreader import Reader

class Transcript:
	def __init__(self,id,chr,source,b,e,score,strand,frame,attributes):
		self.id=id
		self.chr=chr
		self.source=source
		self.b=b
		self.e=e
		self.score=score
		self.strand=strand
		self.frame=frame
		self.attributes=attributes
		self.exons=[]

	def bed(self, options):
		id=self.id
		if options.add_to_id is not None:
			id = "%s~%.2f" % (id,float(self.attributes[options.add_to_id]))
		out = "\t".join((self.chr, str(self.b), str(self.e), id, str(self.score), self.strand, str(self.b), str(self.e), "0"))
		out += "\t%d" % len(self.exons)
		out += "\t" + ",".join((str(e[1]-e[0])   for e in self.exons))
		out += "\t" + ",".join((str(e[0]-self.b) for e in self.exons))
		return out

#chrom	chromStart	chromEnd	name	                score	strand	thickStart	thickEnd	itemRgb	blockCount 	blockSizes	chromStarts
#chr20	50267485	50279795	ENST00000425497.5	1000	+	50279795	50279795	0	4	      178,90,123,1619,	0,8160,10005,10691,
#chr20	50267498	50278427	ENST00000445003.5	1000	+	50278427	50278427	0	3	      165,123,251,	0,9992,10678,
#chr20	50275669	50279795	ENST00000622604.4	1000	+	50279795	50279795	0	3	      66,123,1619,	0,1821,2507,
#chr20	50276316	50279037	ENST00000422459.1	1000	+	50279037	50279037	0	3	      74,123,861,	0,1174,1860,
#chr20	50277533	50278415	ENST00000421019.1	1000	+	50278415	50278415	0	2	      447,239,	   	0,643,



def main():
	usage = format_usage('''
		%prog < some_file.gtf > some_file.bed
	''')
	parser = OptionParser(usage=usage)
	
	parser.add_option('-s', '--score', type=str, dest='score', default=None, help='use the attribute SCORE of the transcript to replace the score column in the outout bed [default: %default]', metavar='SCORE')
	parser.add_option('-a', '--add_to_id', type=str, dest='add_to_id', default=None, help='add to the id the value of the A attribute [default: %default]', metavar='A')
	parser.add_option('-n', '--normalize_score', type=float, dest="normalize_score", default=None, help="divide the score for N and multiply to 1000 (1000 is the maximum score value on ucsc borowser), often used with -s.\nIf this option si given all score that will result >1000 will be truncated =1000. [default: %default]", metavar="N")
	#parser.add_option('-l', '--log_score', action='store_true', dest="log_score", default=False, help="do not report the score but the log2 of the score, often used with -s [default: %default]")
	parser.add_option('-r', '--reference', type=str, dest='reference', default="transcript_id", help='use the attribute REFERENCE to group features [default: %default]', metavar='REFERENCE')
	parser.add_option('-i','--ignore',action='store_true', dest="ignore", default=False, help="ignore rows that do not have REFERENCE attribute")

	options, args = parser.parse_args()

	reference_feature = options.reference
	
	#if len(args) != 2:
	#	exit('Unexpected argument number.')
	
	#for id, sample, raw, norm in Reader(stdin, '0u,1s,2i,3f', False):
	transcripts={}
	warned_saturation=False
	for line in stdin:
		if line[0]=="#":
			print line,
			continue
		line = line.rstrip()

		(seqname,source,feature,b,e,score,strand,frame,attributes) = safe_rstrip(line).split('\t')
		b=int(b)-1 #the bed is 0 based
		e=int(e)
		attributes = [a.strip().split(" ") for a in attributes.split(";") if len(a)>0]
		attributes = [ (a[0], re.sub(r'^"',"", a[1])) for a in attributes] #removing "
		attributes = [ (a[0], re.sub(r'"$',"", a[1])) for a in attributes] #removing "
		attributes = dict(attributes)

		t_id =  attributes.get(reference_feature,None)
		if t_id is None: 
			if options.ignore:
				continue
			raise Exception("Found a row without REFERENCE attibute, use -i to ignore.")
		if options.reference is not None:
			try:
				t_id =  attributes[options.reference]
			except KeyError:
				pass
		if feature=="transcript":
			t = transcripts.get(t_id,None)
			if t is not None:
				raise Exception("ERROR: transcript (%s) already defined" % t_id)
			else:

				if options.score is not None:
					score = attributes[options.score]
				if options.normalize_score is not None:
					score = log(1+float(score))/log(2)
					score = int((score/options.normalize_score)*1000)
					if score > 1000 and not warned_saturation:
						stderr.write("WARNING: score > 1000 saturated at 1000\n")
						warned_saturation=True


				t=Transcript(t_id,seqname,source,b,e,score,strand,frame,attributes)
				transcripts[t_id]=t
		elif feature=="exon":
			transcripts[t_id].exons.append((b,e,attributes))

	for k,v in transcripts.iteritems():
		print v.bed(options)
		
	


if __name__ == '__main__':
	main()

#track name=sample15.stringtie_no_ref
#chr20	StringTie	transcript	50274722	50274954	1000	.	.	gene_id "STRG.32086"; transcript_id "STRG.32086.1"; cov "2.673820"; FPKM "0.384354"; TPM "0.787442";
#chr20	StringTie	exon	50274722	50274954	1000	.	.	gene_id "STRG.32086"; transcript_id "STRG.32086.1"; exon_number "1"; cov "2.673820";
#chr20	StringTie	transcript	50276209	50282279	1000	+	.	gene_id "STRG.32087"; transcript_id "STRG.32087.1"; cov "7.349761"; FPKM "1.056507"; TPM "2.164510";
#chr20	StringTie	exon	50276209	50276390	1000	+	.	gene_id "STRG.32087"; transcript_id "STRG.32087.1"; exon_number "1"; cov "11.019068";
#chr20	StringTie	exon	50277491	50277613	1000	+	.	gene_id "STRG.32087"; transcript_id "STRG.32087.1"; exon_number "2"; cov "29.767330";
#chr20	StringTie	exon	50278177	50282279	1000	+	.	gene_id "STRG.32087"; transcript_id "STRG.32087.1"; exon_number "3"; cov "6.514964";
#chr20	StringTie	transcript	50282335	50282650	1000	.	.	gene_id "STRG.32088"; transcript_id "STRG.32088.1"; cov "2.677215"; FPKM "0.384842"; TPM "0.788442";
#chr20	StringTie	exon	50282335	50282650	1000	.	.	gene_id "STRG.32088"; transcript_id "STRG.32088.1"; exon_number "1"; cov "2.677215";

