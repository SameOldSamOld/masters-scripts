echo HostTreeGen
java -jar jprime-0.3.6.jar HostTreeGen -min 8 -max 8 1 2 0 2017.10.30/1/host_tree
for i in {0..99}
  do
  echo GuestTreeGen
  java -jar jprime-0.3.6.jar GuestTreeGen -min 4 -max 20 2017.10.30/1/host_tree.pruned.tree 2 0 0 2017.10.30/2/run$i
  echo BranchRelaxer
  java -jar jprime-0.3.6.jar BranchRelaxer -o 2017.10.30/3/run$i 2017.10.30/2/run$i.pruned.tree IIDLogNormal 0.0 1.0
  echo Seq-gen
  r=$(( $RANDOM % 1300 + 200 ))
  seq-gen -l $r -n 1 -mPAM < 2017.10.30/3/run$i > 2017.10.30/4/run$i
  #echo RaXML
  #raxmlHPC -m PROTGAMMALGF -p 12345 -s 2017.10.30/4/run$i -n run$i.tree -w /home/sold/Documents/jprime/2017.10.30/5/
  echo run $i
done