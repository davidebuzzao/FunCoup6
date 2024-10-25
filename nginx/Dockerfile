FROM nginx:1.21-alpine

RUN rm /etc/nginx/conf.d/default.conf
COPY nginx.conf /etc/nginx/conf.d
RUN mkdir /home/app
RUN mkdir /home/app/web
RUN mkdir /home/app/web/staticfiles

# https://geko.cloud/en/nginx-letsencrypt-certbot-docker-alpine/
RUN apk add python3 python3-dev py3-pip build-base libressl-dev musl-dev libffi-dev rust cargo
RUN pip3 install pip --upgrade
RUN pip3 install certbot-nginx
RUN mkdir /etc/letsencrypt
# RUN certbot certonly --nginx --email davide.buzzao@scilifelab.se -d funcoup.org --rsa-key-size 4096 --agree-tos
