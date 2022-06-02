REGISTRY=$1
ORG=$2

TESTORG=$3

echo "\n LOGIN \n"
docker login $REGISTRY -u $4 --password $5

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`

docker pull $REGISTRY/$TESTORG/wgs_remixt_hg38:$TAG
docker tag $REGISTRY/$TESTORG/wgs_remixt_hg38:$TAG $REGISTRY/$ORG/wgs_remixt_hg38:$TAG
docker push $REGISTRY/$ORG/wgs_remixt_hg38:$TAG