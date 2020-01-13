import h5py, argparse
from magnetizer import Parameters

if __name__ == "__main__"  :

    parser = argparse.ArgumentParser(description='Extracts input parameter from'
                                     'an Magnetizer output hdf5 file')

    parser.add_argument("MAGNETIZER_OUTPUT",help="Name(s) of the Magnetizer "
                        "output file(s) to be examined.")

    parser.add_argument('-o', "--output", help="Global input parameters file "
                        "that will be generated", default=None)
    args = parser.parse_args()

    report = {}

    h5file = h5py.File(args.MAGNETIZER_OUTPUT, 'r')
    params = Parameters(h5file)
    print('Parameters extracted from ',args.MAGNETIZER_OUTPUT)
    print(params.txt)
    if args.output is not None:
      print('Exporting input parameters to file: ', args.output)
      params.dump(args.output)
