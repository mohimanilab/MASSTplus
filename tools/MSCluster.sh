echo 'Running MSCluster...'
ts="$(date +%s%N)"
/projects/mohimanilab/MSCluster/MSCluster_bin --model LTQ_TRYP --model-dir /projects/mohimanilab/MSCluster/Models --list /projects/mohimanilab/ben/list4.txt --output-name MSCluster_output_test
echo "MSCluster ran in $((($(date +%s%N) - $ts)/1000000))ms" 
