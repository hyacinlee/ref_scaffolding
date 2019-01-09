# ref_scaffolding
tool for scaffolding contigs genome bases on  reference align 

Before the first use, please modify the NUMER path in the script(line9) first.

try help:

python ref_scaffolding.py -h 

	optional arguments:
	  -h, --help  show this help message and exit
		-a          assemblyed contig/scaffold genome (default: None)
	  -r          reference chromsomes genome (default: None)
	  -o          outdir (default: ./output)
	  -t          thread of genome align by nucmer (default: 20)
	  -f          min size to draw dotplot (default: 10000)
	  -s  [ ...]  step choice [1 -> align, 2 -> creat new genome, 3 -> align, 4 ->
	                dotplot ] (default: [1, 2, 3, 4])


