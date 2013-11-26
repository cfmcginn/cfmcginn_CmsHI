if [[ -z "$1" ]]
then
  echo "Usage: ./pmakeDiJetTTree.sh <inputList> <MCBool> <outputFile>"
  exit 1
fi
if [[ -z "$2" ]]
then
  echo "Usage: ./pmakeDiJetTTree.sh <inputList> <MCBool> <outputFile>"
  exit 1
fi
if [[ -z "$3" ]]
then
  echo "Usage: ./pmakeDiJetTTree.sh <inputList> <MCBool> <outputFile>"
  exit 1
fi

now="skimTreeJob_$(date +"%m_%d_%Y__%H_%M_%S")"
mkdir $now
mkdir -p $4
len=`wc -l $1 | awk '{print $1}'`
cp makeDiJetTTree.sh $now
cp $1 $now

NAME="makeDiJetTTree.C"
g++ $NAME $(root-config --cflags --libs) -Werror -Wall -O2 -o "${NAME/%.C/}.exe"
cp makeDiJetTTree.exe $now

cat pmakeDiJetTTree.condor | sed "s@log_flag@$now@g" | sed "s@dir_flag@$PWD/$now@g" | sed "s@arg1@$1@g" | sed "s@arg2@$2@g" | sed "s@arg3@$3@g" | sed "s@arg4@$4@g"  | sed "s@njobs@$len@g" > $now/pmakeDiJetTTree.condor
echo -=-
cat $now/pmakeDiJetTTree.condor
echo -=-