#!/bin/bash
PAUSE=0 NO8=1 NO7=1 make base
./blant -k6 -lDEG1 -mi networks/syeast.el 2>/dev/null >syeast-test.txt
diff -s syeast-test.txt indexes/syeast-k6-l1-dup.ind
rm syeast-test.txt
