# python script to convert the bismark coverage file to a file with tiles of specified length:
#output
# chromosome windowstart windowend methylation nummethylatedreads numunmethylatedreads 
# input
# <chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>

# so the input data just needs to be binned

import sys

windowsize = 100
start = 1
end = 0 + windowsize



chrom = ""
methcount = 0
unmethcount = 0

for line in sys.stdin:

	if (end/windowsize)%100000 == 0:
		sys.stderr.write("processed " + str(end/windowsize) + " tiles \n")


	line = line.strip().split()
	chromosome = line[0]
	pos = int(line[1])
	#end = int(pos[2])
	meth = int(line[4])
	unmeth = int(line[5]) 

	# if new chromosome, start again
	if chromosome != chrom:
		if chrom == "": # first chromosome
			sys.stderr.write("first chromosome; setting chromosome to \n")
			chrom = chromosome
			sys.stderr.write(chrom+"\n")
		else:
			sys.stderr.write("initialising new chromosome\n")

			L = [chrom,start,end, frac, methcount, unmethcount]
			print(*L)
			chrom = chromosome
			start = 0
			end = 0 + windowsize
			methcount = 0
			unmethcount = 0


	if pos >= start and pos < end: # position is in tile, count values


		methcount += meth
		unmethcount += unmeth		
	elif pos >= end: # position in next tile, flush old tile, reset counters for new tile

		if methcount+unmethcount == 0:
			frac = 'NA'
		else:
			frac = float(methcount)/(methcount+unmethcount)
		

		L = [chrom,start,end, frac, methcount, unmethcount]
		print(*L)
		methcount = meth
		unmethcount = unmeth
		start = end +1
		end += windowsize
		
	else:
		sys.stderr.write("WARNING: unknown case\n")
		sys.stderr.write(" ".join(line))
		sys.exit()



