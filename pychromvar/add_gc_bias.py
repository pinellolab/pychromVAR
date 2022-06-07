from pysam import Fastafile

def get_gc_bias(peaks, genome_file, sep="_"):
    fasta = Fastafile(genome_file)
    n_peaks = len(peaks)

    bias = [None] * n_peaks

    for i in range(n_peaks):
        peak = peaks[i].split(sep)
        chrom, start, end = peak[0], int(peak[1]), int(peak[2]) - 1
        seq = fasta.fetch(chrom, start, end).lower()

        freq_a = seq.count("a")
        freq_c = seq.count("c")
        freq_g = seq.count("g")
        freq_t = seq.count("t")

        if freq_a + freq_c + freq_g + freq_t == 0:
            bias[i] = 0.5
        else:
            bias[i] = (freq_g + freq_c) / (freq_a + freq_c + freq_g + freq_t)

    return bias