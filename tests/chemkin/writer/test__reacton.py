"""
 tests writers
"""

import chemkin_io


def test__plog_writer():
    """ test chemkin_io.writer.reaction.plog
    """

    reaction1 = 'CH2+H=CH3'
    reaction2 = 'CH3+H=CH4'
    rate_params_dct1 = {
        1: [111111.1111, 1.11, 111.111],
        10: [222222.2222, 2.22, 222.222],
        100: [333333.3333, 3.33, 333.333],
        'high': [444444.4444, 4.44, 444.444]
    }
    rate_params_dct2 = {
        1: [111111.1111, 1.11, 111.111, 555555.5555, 5.55, 555.555],
        10: [222222.2222, 2.22, 222.222, 666666.6666, 6.66, 666.666],
        100: [333333.3333, 3.33, 333.333, 777777.7777, 7.77, 777.777],
        'high': [444444.4444, 4.44, 444.444, 888888.8888, 8.88, 888.888]
    }
    err_dct = {
        1: [1.11, 11.111],
        10: [2.22, 22.222],
        100: [3.33, 33.333],
        'high': [4.44, 44.444]
    }

    plog_str1 = chemkin_io.writer.reaction.plog(
        reaction1, rate_params_dct1, err_dct)
    plog_str2 = chemkin_io.writer.reaction.plog(
        reaction2, rate_params_dct2, err_dct)
    print('\nplog_str1')
    print(plog_str1)
    print('\nplog_str2')
    print(plog_str2)


if __name__ == '__main__':
    test__plog_writer()
