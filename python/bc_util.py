import Levenshtein as lev

def min_barcode(barcodes, sseq, ind):
    val = 0
    distance_min = 100
    i = 0

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
    return ind
