import sys
import pysam
from collections import defaultdict
import math
import numpy as np


def binned_gc(fasta_path,contig,bin_size,n_cutoff):
    fasta=pysam.FastaFile(fasta_path)
    contig_length=fasta.get_reference_length(contig)
    elements=int(math.ceil(contig_length/bin_size))
	
    contig_gc=np.zeros(elements,dtype=np.int8)

    start=0
    gc_bases=set(["C","c","G","g"])
    for i in range(0,elements):
        slice=fasta.fetch(contig, start, start+bin_size)
        n=0
        gc=0

        for charachter in slice:
             if charachter in gc_bases:
                 gc+=1
             elif charachter in ("N","n"):
                 n+=1

        if n/bin_size > n_cutoff:
             contig_gc[i]=-1

        else:
             contig_gc[i]=round(100*gc/(len(slice)))
        start+=bin_size

    return(contig_gc)

def gc_hist(window_stats,gc_stats,bin_size):
    gc_dictionary={}
        
    for chrom in window_stats:
        if chrom in ("chrY","Y","X","chrX"):
            continue

        for i in range(0,len(window_stats[chrom])):
                if not gc_stats[chrom][i] in gc_dictionary:
                    gc_dictionary[ gc_stats[chrom][i] ]=[]
                
                gc_dictionary[ gc_stats[chrom][i] ].append(window_stats[chrom][i][0])
    hist={}
    for gc in gc_dictionary:
        hist[gc]=-1
        if len(gc_dictionary[gc]) > 50:
            hist[gc]=np.median(gc_dictionary[gc])
    return(hist)

def run(
    bam_path,
    ref_path,
    bin_size=1000,
    mapq_threshold=20,
    out=sys.stdout
):
    bam = pysam.AlignmentFile(bam_path, "rb")
    
    # window_stats[chrom][bin] = [highQ_reads, lowq_reads, total_reads]

    window_stats={}
    bam_header=bam.header
    for chrom in bam_header["SQ"]:
        bins= int(math.ceil(chrom["LN"]/float(bin_size)))
        window_stats[ chrom["SN"] ]=np.zeros([bins,3])

    gc_stats={}
    for chrom in bam_header["SQ"]:
        gc_stats[chrom["SN"]]=binned_gc(ref_path,chrom["SN"],bin_size,0.4)

    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_read2:
            continue

        chrom = bam.get_reference_name(read.reference_id)

        is_highq = read.mapping_quality >= mapq_threshold

        if read.next_reference_id != read.reference_id or abs(read.reference_start-read.next_reference_start) > 500: 
            is_highq = False
    
        # Track read counts (once per read)
        ref_start = read.reference_start
        ref_end = read.reference_end

        bin_start = ref_start // bin_size
        bin_end = (ref_end - 1) // bin_size

        for b in range(bin_start, bin_end + 1):
            window_stats[chrom][b][2] += 1
            if is_highq:
                window_stats[chrom][b][0] += 1
            else:
                window_stats[chrom][b][1] += 1

    gc_histogram=gc_hist(window_stats,gc_stats,bin_size)

    # Output BED
    for chrom in sorted(window_stats):
        for i in range(0,len(window_stats[chrom])):
            start = i * bin_size
            end = start + bin_size

            highq_bases, lowq_reads, total_reads = window_stats[chrom][i]


            mean_cov = highq_bases
            gc=gc_stats[chrom][i]
            if gc < 0:
                frac_cov=0.0
            else: 
                frac_cov = mean_cov / gc_histogram[gc_stats[chrom][i]] if gc_histogram[gc_stats[chrom][i]] > 0 else 0.0
            lowq_frac = lowq_reads / total_reads if total_reads > 0 else 0.0

            out.write(
                f"{chrom}\t{start}\t{end}\t"
                f"{mean_cov:.6f}\t{frac_cov:.6f}\t{lowq_frac:.6f}\n"
            )

    bam.close()
