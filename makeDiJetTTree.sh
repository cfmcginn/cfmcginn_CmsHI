if [[ -z "$1" ]]
then
  echo "Usage: ./makeDiJetTTree.sh <inputList> <MCBool> <outputFile>"
  exit 1
fi
if [[ -z "$2" ]]
then
  echo "Usage: ./makeDiJetTTree.sh <inputList> <MCBool> <outputFile>"
  exit 1
fi
if [[ -z "$3" ]]
then
  echo "Usage: ./makeDiJetTTree.sh <inputList> <MCBool> <outputFile>"
  exit 1
fi


echo | awk -v inputList=$1 -v MCBool=$2 -v outputFile=$3 '{print "./makeDiJetTTree.exe \""inputList"\" \""MCBool"\" \""outputFile"\""}' | bash
mv *.root $4

echo "job done successfully"