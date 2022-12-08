#/bin/sh
USAGE="$0 [-renormalize] testEdges.el blant-mp-outfile
PURPOSE: compute stats across orbit pairs like mean, stdDev, etc.
If the -renormalize option is given, then the re-normalized version is given on stdout,
and arranges that, for each orbit pair, the mean value across training edges is 1, and
that *all* non-edges in the training set are normalized so that edges in the *test* set
also have a mean of 1."

# Functions
die(){ (echo "USAGE: $USAGE"; echo "`basename $0`: FATAL ERROR: $@")>&2; exit 1; }

# generally useful Variables
NL='
'
TAB='	'
BASENAME=`basename "$0" .sh`
USAGE=${USAGE:?"$0 should have a USAGE message before sourcing skel.sh"}
echo "$BASENAME" | grep "[ $TAB]" && die "Shell script name ('$BASENAME') cannot contain spaces or tabs"
[ $BASENAME == "$BASENAME" ] || die "something weird with filename in '$BASENAME'"

# Temporary Filename + Directory (both, you can use either, note they'll have different random stuff in the XXXXXX part)
TMPDIR=`mktemp -d /tmp/$BASENAME.XXXXXX`
trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N

RENORM=0
case "$1" in
-renorm*) RENORM=1; shift;;
-*) die "unknown option '$1'";;
esac

[ $# -eq 2 ] || die "expecting exactly 2 filenames"
[ -f "$1" ] || die "cannot open test edge list file '$1'"
TEST="$1"; shift 1;

[ -f "$1" ] || die "cannot open blant-outfile '$1'"

SECOND=""

B=`basename "$1" | sed "s/[ $TAB]//g"`
wzcat "$1" > $TMPDIR/$B
[ $RENORM -eq 1 ] && SECOND=$TMPDIR/$B
hawk 'BEGIN{statout="/dev/stdout"}
    ARGIND==1{T[$1":"$2]=1}
    ARGIND==2{
	getAv=0;
	if($2){getAv=1; ASSERT(!($1 in T), $1" is "$2" in training but is also in the test set")}
	if($1 in T){getAv=1; ASSERT(!$2, $1" is in the test set but is "$2" in training")}
	for(i=3;i<NF;i+=2) {
	    pairs[$i]=1;
	    if(getAv) StatAddSample($i" "$2,$(i+1));
	}
    }
    ENDFILE{if(ARGIND==2) for(i in pairs) if(_statN[i" "0] && _statN[i" "1]) ratio[i] = StatMean(i" "0)/StatMean(i" "1)
	if(ARGIND==3) statout="/dev/stderr"
    }
    ARGIND==3{
	printf "%s %d", $1,$2
	for(i=3;i<NF;i+=2) {
	    op=$i;
	    if(op in ratio) {
		if(1*$2) factor = ratio[op]; # reduce the ratio of the true edges
		else factor = 1;
		norm = StatMean(i" "1); # normalize so that training edges have a mean of 1
		printf "\t%s %g", op, factor*$(i+1)/norm
	    } else ++noRatio[op];
	}
	print ""
    }
    END{
	for(i in pairs) if(_statN[i" "0] && _statN[i" "1])
	    printf "%8.6g %11s %8.6g %8.6g %8.6g %8.6g ratio\n", ratio[i], i,
		StatMean(i" "0), StatMean(i" "1), StatStdDev(i" "0), StatStdDev(i" "1) > statout
	for(op in noRatio)
	    printf "WARNING: %d values of %s were missed since no ratio could be computed\n",noRatio[op],op > statout
    }' "$TEST" "$TMPDIR/$B" $SECOND
