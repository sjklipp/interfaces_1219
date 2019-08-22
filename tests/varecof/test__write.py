"""
  Test writing the various input file
"""

import varecof_io.writer


def test__writer():
    """ test writing
    """

    # Write the tst input string
    nsamp_max = 2000
    nsamp_min = 500
    flux_err = 5
    pes_size = 1
    tst_inp_str = varecof_io.writer.write_tst_input(
        nsamp_max, nsamp_min, flux_err, pes_size)
    print('\ntst.inp:')
    print(tst_inp_str)

    # Write the divsur input file string; distances in Angstrom
    distances = [3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]
    divsur_inp_str = varecof_io.writer.write_divsur_input(
        distances)
    print('\ndivsur.inp:')
    print(divsur_inp_str)

    # Write the els input string
    exe_path = '/path/to/exe'
    base_name = 'mol'
    els_inp_str = varecof_io.writer.write_els_input(
        exe_path, base_name)
    print('\nels.inp:')
    print(els_inp_str)

    # Write the structure input string
    geom1 = (('H', (0.0, 0.0, 0.0)),)
    geom2 = (('C', (0.2115677758, -0.4050266480, 0.0238323931)),
             ('H', (0.1997580074, 0.1613210912, -0.9300443271)),
             ('H', (0.6278380682, 0.2349826345, 0.8287405717)),
             ('H', (0.8417991162, -1.3107013202, -0.0916314599)))
    struct_inp_str = varecof_io.writer.write_structure_input(
        geom1, geom2)
    print('\nstructure.inp:')
    print(struct_inp_str)

    # Write the *.tml input string
    memory = 4.0
    basis = 'cc-pvdz'
    wfn = """{uhf,maxit=300;wf,78,1,2}
   {multi,maxit=40;closed,38;occ,40;wf,78,1,0;orbprint,3}
   {multi,maxit=40;closed,37;occ,40;wf,78,1,0;state,2;orbprint,3}"""
    method = '{rs2c, shift=0.25}'
    inf_sep_energy = -654.3210123456
    tml_inp_str = varecof_io.writer.write_tml_input(
        memory, basis, wfn, method, inf_sep_energy)
    print('\nmol.tml:')
    print(tml_inp_str)


if __name__ == '__main__':
    test__writer()
