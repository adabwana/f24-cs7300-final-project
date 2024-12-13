# chmod +x docker-kill-all.sh

docker system prune -f
docker volume prune -f
docker rm -f -v $(docker ps -q -a)
