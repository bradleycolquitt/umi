import Levenshtein as lev

cpdef int min_barcode(list barcodes, char* sseq, int ind):
    cdef int val = 0
    cdef int distance_min = 100
    cdef int i = 0

    for bc in barcodes:
        val = lev.distance(sseq, bc)
        if val == 0:
            distance_min = val
            ind = i
            break
        if val < distance_min:
            distance_min = val
            ind = i
        i += 1

    if distance_min > 0: ind = 0
    return 0