import biom
import sys
from collections import defaultdict
from operator import itemgetter


def get_inv_otu_map(lines):
    """Map a cluster member to its GG ID"""
    otu_map = {line[0]: line[1:] for line in lines}

    inv_otu_map = {}
    for gg_id, deblur_seqs in otu_map.items():
        for seq in deblur_seqs:
            inv_otu_map[seq] = gg_id

    return inv_otu_map


def next_index(current_index, id_):
    """Get the next index position starting at zero

    Modifies current_index in place
    """
    if id_ in current_index:
        idx = current_index[id_]
    else:
        idx = len(current_index)
        current_index[id_] = idx
    return idx


def index_to_order(index):
    """Returns the ascending order of the keys by the index values"""
    return tuple(id_ for id_, idx in sorted(index.items(), key=itemgetter(1)))


def perform_expansion(deblur_table, inv_otu_map):
    """Expand a deblurred table based on an OTU map

    ...perhaps this should be called collapse... algorithm is as follows:

    For each sequence in the deblurred table:
        If the sequence exists in the OTU map:
            get the representative ID
            for each sample:
                get the count of the deblurred sequence in the sammple
                increment the count of (representative ID, sample)
    """
    new_table_rcv = defaultdict(int)
    ggid_to_index = {}
    sample_to_index = {}

    for deblur_seq, sample in deblur_table.nonzero():
        ggid = inv_otu_map.get(deblur_seq, None)

        # if the deblurred sequence isn't in the OTU map (i.e., didn't recruit)
        if ggid is None:
            continue

        ggid_idx = next_index(ggid_to_index, ggid)
        sample_idx = next_index(sample_to_index, sample)

        value = deblur_table.get_value_by_ids(deblur_seq, sample)
        new_table_rcv[(ggid_idx, sample_idx)] += value

    ggid_order = index_to_order(ggid_to_index)
    sample_order = index_to_order(sample_to_index)
    return biom.Table(new_table_rcv, ggid_order, sample_order)


if __name__ == '__main__':
    deblur_table = biom.load_table(sys.argv[1])
    inv_otu_map = get_inv_otu_map([l.strip().split('\t')
                                   for l in open(sys.argv[2])])
    output = sys.argv[3]

    expanded = perform_expansion(deblur_table, inv_otu_map)

    import h5py
    with h5py.File(output, 'w') as fp:
        expanded.to_hdf5(fp, 'deblur remapped')
