import argparse
from RLW_import_DELTAMED import import_deltamed

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--path", type=str, help="Path of the eeg file")

args = parser.parse_args()

out_header, out_data, message_string = import_deltamed(args.path)
