## Author: Derek Rothenheber
## Date: 10-12-18

"""
Python script to parse genome segments of Influenza A to remove smaller genes (potential isoforms)

1. The final output will be 8 genomic segments
2. Order of segments will be PB2, PB1, PA, HA, NP, NA, M, NS1  (haven't figured that aspect out yet...)

"""

import sys


def fasta_dict(filename):
    """
    Reads a fasta file and then creates a contigs dictionary from that file

    :param filename: input fasta file
    :return: dictionary, where headers are ">contig " and values are the corresponding sequence
    """

    ## dictionaries used for filtering out longest hit later.

    pb2_contigs = {}  ## 1: Polymerase basic protein 2
    pb1_contigs = {}  ## 2: RNA-directed RNA polymerase catalytic subunit
    pa_contigs = {}   ## 3: Polymerase acidic protein
    ha_contigs = {}   ## 4: Hemagglutinin
    np_contigs = {}   ## 5: Nucleoprotein
    na_contigs = {}   ## 6: Neuraminidase
    m_contigs = {}    ## 7: Matrix protein 1
    ns1_contigs = {}  ## 8: Non-structural protein 1

    filtered_contigs = {}
    final_contigs = {}

    contigs = {}                    ## empty dictionary
    header = ""                     ## empty string
    f = open(filename, "r")         ## opens file and reads it
    for line in f:                  ## loops through the lines of the file
        if line[0] == ">":                  ## boolean, makes sure the first line starts with ">"
            header = line.rstrip()          ## variable header equals the line strip of whitespace characters
            contigs[header] = ""            ## assigns header as key in contigs and value as empty string
        else:
            contigs[header] += line.rstrip()    ## adds sequence lines that are assoicated with the header, with
                                                    # with white space characters removed
    f.close()                       ## closes file and returns dictionary
    #print(contigs)


    key_to_value_lengths = {key: len(value) for key, value in contigs.items()}
    #print(key_to_value_lengths)


    for (k,v), (k2, v2) in zip(contigs.items(), key_to_value_lengths.items()):
        if "Polymerase basic protein 2" in k2:
            pb2_contigs[k2] = v2
        if "RNA-directed RNA polymerase catalytic subunit" in k2:
           pb1_contigs[k2] = v2
        if "Polymerase acidic protein" in k2:
            pa_contigs[k2] = v2
        if "Hemagglutinin" in k2:
            ha_contigs[k2] = v2
        if "Nucleoprotein" in k2:
            np_contigs[k2] = v2
        if "Neuraminidase" in k2:
            na_contigs[k2] = v2
        if "Matrix protein 1" in k2:
            m_contigs[k2] = v2
        if "Non-structural protein 1" in k2:
            ns1_contigs[k2] = v2


    print(pb2_contigs)
    print(pb1_contigs)
    print(pa_contigs)
    print(ha_contigs)
    print(np_contigs)
    print(na_contigs)
    print(m_contigs)
    print(ns1_contigs)


    maximum_pb2 = max(pb2_contigs.values())
    results_pb2 = filter(lambda x: x[1] == maximum_pb2, pb2_contigs.items())
    for x in results_pb2:
        filtered_contigs[x[0]] = x[1]

    maximum_pb1 = max(pb1_contigs.values())
    results_pb1 = filter(lambda x: x[1] == maximum_pb1, pb1_contigs.items())
    for x in results_pb1:
        filtered_contigs[x[0]] = x[1]

    maximum_pa = max(pa_contigs.values())
    results_pa = filter(lambda x: x[1] == maximum_pa, pa_contigs.items())
    for x in results_pa:
        filtered_contigs[x[0]] = x[1]

    maximum_ha = max(ha_contigs.values())
    results_ha = filter(lambda x: x[1] == maximum_ha, ha_contigs.items())
    for x in results_ha:
        filtered_contigs[x[0]] = x[1]

    maximum_np = max(np_contigs.values())
    results_np = filter(lambda x: x[1] == maximum_np, np_contigs.items())
    for x in results_np:
        filtered_contigs[x[0]] = x[1]

    maximum_na = max(na_contigs.values())
    results_na = filter(lambda x: x[1] == maximum_na, na_contigs.items())
    for x in results_na:
        filtered_contigs[x[0]] = x[1]

    maximum_m = max(m_contigs.values())
    results_m = filter(lambda x: x[1] == maximum_m, m_contigs.items())
    for x in results_m:
        filtered_contigs[x[0]] = x[1]

    maximum_ns1 = max(ns1_contigs.values())
    results_ns1 = filter(lambda x: x[1] == maximum_ns1, ns1_contigs.items())
    for x in results_ns1:
        filtered_contigs[x[0]] = x[1]


    for k, v in contigs.items():
        if k in filtered_contigs:
            print(True)
            final_contigs[k] = v
        else:
            print(False)

#    print(final_contigs)

    final_file = open(sys.argv[2], 'w')
    string = ""
    for a, b in final_contigs.items():
        fasta = a + '\n' + b + '\n'
        string += fasta

    final_file.write(string)
    final_file.close()

    #print(final_contigs)
    return contigs


if __name__ == '__main__':
    #filename = "test.fasta"
    filename = sys.argv[1]
    contigs = fasta_dict(filename)
