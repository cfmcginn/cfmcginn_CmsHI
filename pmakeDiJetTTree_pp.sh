if [ $# -ne 4 ]
then 
  echo "Usage: ./makeDiJetTTree_pp.sh <inputList> <MCBool> <outName> <outDir>"
  exit 1
fi

now="skimTreeJob_$(date +"%m_%d_%Y__%H_%M_%S")"
mkdir $now
mkdir -p $4
len=`wc -l $1 | awk '{print $1}'`
cp makeDiJetTTree_pp.sh $now
cp ptCorrDir_pp/*.tar.gz $now
cp $1 $now

NAME="makeDiJetTTree_pp.C"
g++ $NAME $(root-config --cflags --libs) -Werror -Wall -O2 -o "${NAME/%.C/}.exe"
cp makeDiJetTTree_pp.exe $now

cat pmakeDiJetTTree_pp.condor | sed "s@log_flag@$now@g" | sed "s@dir_flag@$PWD/$now@g" | sed "s@arg1@$1@g" | sed "s@arg2@$2@g" | sed "s@arg3@$3@g" | sed "s@arg4@$4@g"  | sed "s@njobs@$len@g" > $now/pmakeDiJetTTree_pp.condor
echo -=-
cat $now/pmakeDiJetTTree_pp.condor
echo -=-
