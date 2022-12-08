#/bin/sh
# Do whatever is necessary to load your favorite version of python on the line(s) below
module load python/3.6.8

# Now tell me the name of the python interpreter that you want to use, given the above
PYTHON=python3

# PYTHON GENIUSES DO NOT TOUCH ANYTHING BETWEEN HERE, AND WHERE IT SAYS "Put your fucking Python shit ..."
#############################################################################################################

USAGE="$0 bla bla bla
PURPOSE: Run a python script included at the end of the Bourne Shell script."

BASENAME=`basename "$0" .sh`
[ $BASENAME == "$BASENAME" ] || die "something weird with filename in '$BASENAME'"

# Temporary Filename + Directory (both, you can use either, note they'll have different random stuff in the XXXXXX part)
TMP=`mktemp /tmp/$BASENAME.XXXXXX`
TMPDIR=`mktemp -d /tmp/$BASENAME.XXXXXX`
trap "/bin/rm -rf $TMP $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N

sed '1,/^#Put/d' "$0" > $TMPDIR/fuck-python.py
chmod +x $TMPDIR/fuck-python.py

exec $PYTHON $TMPDIR/fuck-python.py "$@"
exit 1


#############################################################################################################
#Put your fucking Python shit below this line.

import numpy as np
import sys
from sklearn.metrics import average_precision_score as prc
from sklearn.metrics import roc_auc_score as auc
from sklearn.metrics import ndcg_score as ndcg

def generate_result():
    true_labels = []
    confidence = []
    filename = sys.argv[1]
    open_file = open(filename,"r")
    count = 0
    p_500 = 0
    for line in open_file:
        count += 1
        track = 0
        line = line.split("\t")
        for x in line:
            x = x.split()
            if track == 0:
                true_labels.append(int(x[0]))
                if count == 501:
                    p_500 = float(x[3])
            elif track == 1:
                confidence.append(float(x[0]))
            track += 1
    open_file.close()

    true_labels = np.array(true_labels, dtype = np.int32)
    confidence = np.array(confidence, dtype = np.float64)
    flabels = 1 - true_labels
    fconf = 1 - confidence
    all_labels = np.column_stack((flabels,true_labels))
    all_conf = np.column_stack((fconf,confidence))
    out_file = open("output.txt", 'a')
    result = ""
    result = result + str(auc(true_labels, confidence)) + "\t"
    result = result + str(prc(true_labels, confidence)) + "\t"
    result = result + str(p_500) + "\t"
    result = result + str(ndcg(all_labels, all_conf)) + "\t" 
    out_file.write(result)
    out_file.write("\n")
    # print(result)
    out_file.close()
if __name__ == "__main__":
    generate_result()
    # print("Added output from file ", filename)

