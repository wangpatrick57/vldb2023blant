#!/pkg/python/3.7.4/bin/python3
import sys
import datetime
from graph_helpers import *

def convert_btc(csv_path, is_tel):
    with open(csv_path, 'r') as btc_file:
        xel = []

        for line in btc_file:
            source, target, rating, time = line.strip().split(',')
            edge = [source, target]

            if is_tel:
                edge.append(str(int(float(time))))

            edge = tuple(edge)
            xel.append(edge)

        print_xel(xel)

def convert_reddit(tsv_path, is_tel):
    fmt = '%Y-%m-%d %H:%M:%S'

    with open(tsv_path, 'r') as reddit_file:
        xel = []

        for line in reddit_file:
            source, target, post, time_str, label, props = line.strip().split('\t')
            time = int(datetime.datetime.strptime(time_str, fmt).timestamp())
            edge = [source, target]

            if is_tel:
                edge.append(str(int(float(time))))

            edge = tuple(edge)
            xel.append(edge)

        print_xel(xel)

if __name__ == '__main__':
    tsv_path = sys.argv[1]
    convert_reddit(tsv_path, True)
