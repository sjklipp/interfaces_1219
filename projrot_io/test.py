"""
test
"""

# import sys
import numpy as np

HESS = np.random.rand(33, 33)
NROWS, NCOLS = HESS.shape

if NROWS % 5 == 0:
    NCHUNKS = NROWS // 5
else:
    NCHUNKS = (NROWS // 5) + 1

print(HESS)
print('\n\n')

HESS_STR = '   ' + '    '.join([str(val) for val in range(1, 6)]) + '\n'
N = 0
while N+1 <= NCHUNKS:
    for i in range(NROWS):
        col_tracker = 1
        if i >= 5*N:
            HESS_STR += '{0}'.format(str(i+1))
            for j in range(5*N, NCOLS):
                if i >= j:
                    if col_tracker <= 5:
                        HESS_STR += '  {0:5.3f}'.format(HESS[i][j])
                        col_tracker += 1
                        if col_tracker == 6:
                            HESS_STR += '\n'
                    else:
                        continue
                elif i < j and col_tracker != 6:
                    HESS_STR += '\n'
                    break
                else:
                    break
        if i+1 == NROWS and N+1 < NCHUNKS-1:
            HESS_STR += '    ' + '    '.join([str(val) for val in range(5*(N+1) + 1, 5*(N+1) + 6)]) + '\n'
        if i+1 == NROWS and N+1 == NCHUNKS-1:
            HESS_STR += '    ' +  '    '.join([str(val) for val in range(5*(N+1) + 1, NROWS+1)]) + '\n'
    N += 1

print('\n\nhess string:\n')
print(HESS_STR)
print('\n\n')
