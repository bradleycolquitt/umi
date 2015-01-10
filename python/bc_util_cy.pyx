import fasta_cy as fasta
import Levenshtein as lev

cdef int min_barcode(list barcodes, str sseq):
    cdef int val = 0
    cdef int distance_min = 100
    cdef int i = 0
    cdef int ind = 0

    for bc in barcodes:
        val = lev.distance(sseq, bc)
        if val == 0:
            ind = i
            
            return ind
        if val < distance_min:
            distance_min = val
            ind = i
        i += 1

    if distance_min > 0: ind = -1
    return ind

    
cpdef int get_barcode(str seq, char[:] qual, sub_range, nofilter, barcodes):
    cdef int squal_min
    cdef str sseq
    if not nofilter:
        #  First, check qual. If minimum PHRED is less than 20 anywhere in barcode discard
        squal_min = min(qual[sub_range[0]:sub_range[1]])
        if squal_min <= 20: return -1

    # Use Levenshtein distance to identify minimum scoring barcode
    sseq_p = seq[sub_range[0]:sub_range[1]]
    sseq = sseq_p
    return min_barcode(barcodes, sseq)

cpdef str get_umi(str seq, char[:] qual, sub_range, nofilter):
    cdef int min_qual
    cdef str sseq
    if not nofilter:
        min_qual = min(qual[sub_range[0]:sub_range[1]])
        if min_qual < 17:
            return 'N'
    return seq[sub_range[0]:sub_range[1]]