#!/usr/bin/env python
import h5py
import sys



def copy_input_to_output(input_file, output_file, final_file=None):
    input_h5 = h5py.File(input_file,'r')
    output_h5 = h5py.File(output_file)

    if final_file:
        final_h5 = h5py.File(final_file)
        for x in ('Output', 'Log'):
            if x in final_h5:
              del final_h5[x]
            output_h5.copy(x, final_h5)
    else:
        final_h5 = output_h5

    for x in ('Input', 'Galform Parameters'):
        if x in final_h5:
          del final_h5[x]
        input_h5.copy(x, final_h5)

    output_h5.close()
    input_h5.close()
    if final_file:
        final_h5.close()
    return


if __name__ == "__main__"  :
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    if len(sys.argv) > 3:
        final_file = sys.argv[3]
    else:
        final_file = None

    copy_input_to_output(input_file, output_file, final_file=final_file)
