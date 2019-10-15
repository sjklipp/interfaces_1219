"""
  Tests the varecof_io.writer functions
"""

import varecof_io.reader


def test__divsur_frag_geoms_reader():
    """ tests varecof_io.reader.divsur.frag_geoms_divsur_frame
    """
    with open('divsur.out', 'r') as divsur_file:
        divsur_string = divsur_file.read()

    frag_geoms = varecof_io.reader.divsur.frag_geoms_divsur_frame(
        divsur_string)
    print(frag_geoms)


if __name__ == '__main__':
    test__divsur_frag_geoms_reader()
