# chmod +x devpod-kill.sh

docker ps -a | grep devpod | awk '{print $1}' | xargs -r docker rm -f
docker volume ls | grep devpod | awk '{print $2}' | xargs -r docker volume rm
docker images | grep devpod | awk '{print $3}' | xargs -r docker rmi -f
docker system prune -af --volumes
