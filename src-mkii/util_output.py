import re
import numpy as np


out_nodes_per_line = 16
center_column = 0


def format_str(v):
    string = ' {: 5.4f}'.format(v)
    return re.sub(
        r'(?<=\.)(\d+?)(0+)(?=[^\d]|$)',
        lambda m: m.group(1) + ' ' * len(m.group(2)),
        string)


np.set_printoptions(precision=4,
                    linewidth=175,
                    suppress=True,
                    formatter={'float': format_str})


def print_vector(v, model_time=None, name=None):

    w = v.shape[0]

    half = int(out_nodes_per_line / 2)
    left = center_column - half
    right = center_column + half
    left = left if left >= 0 else 0
    right = right if right < w else w-1

    out = v[left:right]

    if model_time is None:
        if name is None:
            print(out)
        else:
            print('{} = '.format(name), out)
    else:
        if name is None:
            print('t = {:<6.2f} '.format(model_time), out)
        else:
            print('t = {:<6.2f}  {} = '.format(model_time, name), out)


def print_matrix(m, model_time=None, name=None):

    w = m.shape[1]

    half = int(out_nodes_per_line / 2)
    left = center_column - half
    right = center_column + half
    left = left if left >= 0 else 0
    right = right if right < w else w - 1

    out = m[left:right, left:right]

    if model_time is None:
        if name is None:
            print(out)
        else:
            print('{} =\n'.format(name), out)
    else:
        if name is None:
            print('t = {:6.2f}\n'.format(model_time), out)
        else:
            print('t = {:6.2f}\t{} =\n'.format(model_time, name), out)
