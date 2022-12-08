#!/bin/sh
module unload gcc; module avail gcc 2>&1 | newlines | fgrep gcc/ | grep -v local | sed 's,.*/,,' |
    while read GCC; do
	if module load gcc/$GCC 2>/dev/null; then
	    echo "module unload gcc; module load gcc/$GCC 2>/dev/null; (../../libwayne/bin/wgcc -c blant-predict.c -o ../blant-predict.Linux.gcc$GCC.o)"
	fi
    done
