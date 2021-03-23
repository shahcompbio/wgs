REGISTRY=$1
ORG=$2

TESTORG=$3

echo "\n LOGIN \n"
docker login $REGISTRY -u $4 --password $5

docker pull $REGISTRY/$TESTORG/wgs_qc:$TAG
docker tag $REGISTRY/$TESTORG/wgs_qc:$TAG $REGISTRY/$ORG/wgs_qc:$TAG
docker push $REGISTRY/$ORG/wgs_qc:$TAG

