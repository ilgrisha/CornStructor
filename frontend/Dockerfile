# File: frontend/Dockerfile
# Version: v0.4.0
# Build uses npm install to tolerate a stale lock; switch back to `npm ci` once lock is clean.
FROM node:20-bookworm AS build
WORKDIR /app
COPY frontend/package*.json ./
RUN npm install
COPY frontend ./
# Force a predictable output dir regardless of project name
RUN npx ng build --configuration production --output-path /app/dist

FROM nginx:1.27-alpine AS runtime
COPY nginx/nginx.conf /etc/nginx/nginx.conf
COPY --from=build /app/dist /usr/share/nginx/html
RUN mkdir -p /usr/share/nginx/html/reports
EXPOSE 8080
CMD ["nginx", "-g", "daemon off;"]
