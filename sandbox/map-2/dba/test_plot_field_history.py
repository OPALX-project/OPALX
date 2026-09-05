"""Element-strip parsing preserves SDDS transitions, quoted names and overlaps."""
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from plot_field_history import read_element_spans


class ElementSpansTest(unittest.TestCase):
    def read_rows(self, rows, units='m'):
        text = ('SDDS1\n&parameter name=revision, type=string &end\n'
                f'&column name=s, type=double, units={units} &end\n'
                '&column name=element_names, type=string &end\n'
                '&data mode=ascii, no_row_counts=1 &end\nrevision with spaces\n' + rows)
        with TemporaryDirectory() as temporary:
            path = Path(temporary) / 'positions.sdds'
            path.write_text(text)
            return read_element_spans(path)

    def test_overlaps_and_duplicate_transition_rows(self):
        frame = self.read_rows('0 ""\n0 "B1"\n1 "B1"\n1 "D1, B1"\n2 "D1, B1"\n'
                               '2 ""\n2 "QACH"\n2.2 "QACH"\n2.2 ""\n')
        spans = frame.set_index('element')
        self.assertEqual(tuple(spans.loc['B1']), (0, 2, 0))
        self.assertEqual(tuple(spans.loc['D1']), (1, 2, 1))
        self.assertEqual(tuple(spans.loc['QACH']), (2, 2.2, 0))

    def test_separated_visits_are_not_bridged(self):
        frame = self.read_rows('0 "B1"\n1 "B1"\n1 ""\n2 ""\n2 "B1"\n3 "B1"\n')
        self.assertEqual(frame[['s_begin_m', 's_end_m']].values.tolist(), [[0, 1], [2, 3]])

    def test_touching_elements_share_one_lane(self):
        frame = self.read_rows('0 "B1"\n1 "B1"\n1 ""\n1 "D1"\n2 "D1"\n'
                               '2 ""\n2 "QACH"\n2.2 "QACH"\n2.2 ""\n')
        self.assertEqual(frame.element.tolist(), ['B1', 'D1', 'QACH'])
        self.assertEqual(frame.lane.tolist(), [0, 0, 0])
        self.assertEqual(frame.s_end_m.tolist()[:-1], frame.s_begin_m.tolist()[1:])

    def test_invalid_units_or_unsorted_positions_rejected(self):
        with self.assertRaises(ValueError):
            self.read_rows('0 "B1"\n1 "B1"\n', units='mm')
        with self.assertRaises(ValueError):
            self.read_rows('1 "B1"\n0 "B1"\n')


if __name__ == '__main__':
    unittest.main()
