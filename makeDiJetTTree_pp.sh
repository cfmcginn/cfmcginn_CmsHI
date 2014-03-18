if [ $# -ne 5 ]
then 
  echo "Usage: ./makeDiJetTTree_pp.sh <inputList> <MCBool> <outputFile> <outDir> <#>"
  exit 1
fi


tar -xzvf corrFilePP_20140317.tar.gz
echo | awk -v inputList=$1 -v MCBool=$2 -v outputFile=$3 -v num=$5 '{print "./makeDiJetTTree_pp.exe \""inputList"\" \""MCBool"\" \""outputFile"\" \""num"\""}' | bash
rm eff*.root
rm fake*.root
rm corrFilePP_20140317.tar.gz
mv $3_$5.root $4
rm *.root
rm *.txt

echo "job done successfully"