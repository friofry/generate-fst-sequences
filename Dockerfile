FROM nginx:alpine

# Copy the contents of the local "html" directory to the Nginx document root
COPY html /usr/share/nginx/html

# Copy the updated default.conf file to the Nginx configuration directory
COPY nginx.conf /etc/nginx/conf.d/default.conf
