from unittest import TestCase, main

import biom
import numpy as np

from expand import (get_inv_otu_map, next_index, index_to_order,
                    perform_expansion)


class ExpandTests(TestCase):
    def test_get_inv_otu_map(self):
        otu_map = [['foo', '1', '2'], ['bar', 'x']]
        exp = {'1': 'foo', '2': 'foo', 'x': 'bar'}
        obs = get_inv_otu_map(otu_map)
        self.assertEqual(obs, exp)

    def test_perform_expansion(self):
        otu_map = [['foo', 'O1'], ['bar', 'O5'], ['baz', 'O3', 'O4']]
        inv_map = get_inv_otu_map(otu_map)
        table = biom.Table(np.arange(15).reshape(5, 3),
                           ['O1', 'O2', 'O3', 'O4', 'O5'],
                           ['S1', 'S2', 'S3'])
        exp = biom.Table(np.array([[0, 1, 2],
                                   [12, 13, 14],
                                   [15, 17, 19]]),
                         ['foo', 'bar', 'baz'], ['S1', 'S2', 'S3'])
        obs = perform_expansion(table, inv_map)
        obs = obs.sort_order(exp.ids())
        obs = obs.sort_order(exp.ids(axis='observation'), axis='observation')
        self.assertEqual(obs, exp)

    def test_next_index(self):
        starting = {'foo': 0, 'bar': 1}

        exp1 = 0
        obs1 = next_index(starting, 'foo')

        exp2 = 1
        obs2 = next_index(starting, 'bar')

        exp3 = 0
        obs3 = next_index(starting, 'foo')

        exp4 = 2
        obs4 = next_index(starting, 'baz')

        exp5 = 0
        obs5 = next_index(starting, 'foo')

        self.assertEqual(obs1, exp1)
        self.assertEqual(obs2, exp2)
        self.assertEqual(obs3, exp3)
        self.assertEqual(obs4, exp4)
        self.assertEqual(obs5, exp5)

    def test_index_to_order(self):
        starting = {'foo': 1, 'bar': 0, 'baz': 2}
        exp = ('bar', 'foo', 'baz')
        obs = index_to_order(starting)
        self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()
