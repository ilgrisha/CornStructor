# File: Makefile
# Version: v0.7.0
.PHONY: help dev dev-down prod-build prod-up prod-down logs doctor ps

# You can override these on the command line, e.g.:
# make dev BACKEND_PORT=9000 FRONTEND_PORT_DEV=4300
# make prod-up NGINX_PORT=9090
BACKEND_PORT ?= 8016
FRONTEND_PORT_DEV ?= 4216
NGINX_PORT ?= 8098
CORS_ORIGINS ?= *
COMPOSE_FILE_DEV := docker-compose.dev.yml
COMPOSE_FILE_PROD := docker-compose.yml
# Optional: allow multiple concurrent stacks by project name
PROJECT ?= cornstructor

help:
	@echo "Targets:"
	@echo "  dev        - start live-reload dev stack"
	@echo "  dev-down   - stop dev stack and remove volumes"
	@echo "  prod-build - build production images"
	@echo "  prod-up    - start production stack (Nginx on $${NGINX_PORT})"
	@echo "  prod-down  - stop production stack and remove volumes"
	@echo "  logs       - tail logs from all services"
	@echo "  ps         - show running services"
	@echo "  doctor     - diagnose Docker/compose environment"
	@echo ""
	@echo "Port overrides (examples):"
	@echo "  make dev BACKEND_PORT=9001 FRONTEND_PORT_DEV=4301"
	@echo "  make prod-up NGINX_PORT=9090"

dev:
	@echo "== Starting dev stack =="
	@echo "BACKEND_PORT=$(BACKEND_PORT)"
	@echo "FRONTEND_PORT_DEV=$(FRONTEND_PORT_DEV)"
	@echo "CORS_ORIGINS=$(CORS_ORIGINS)"
	BACKEND_PORT=$(BACKEND_PORT) \
	FRONTEND_PORT_DEV=$(FRONTEND_PORT_DEV) \
	CORS_ORIGINS="$(CORS_ORIGINS)" \
	docker compose -p $(PROJECT)-dev -f $(COMPOSE_FILE_DEV) up --build

dev-down:
	docker compose -p $(PROJECT)-dev -f $(COMPOSE_FILE_DEV) down -v

prod-build:
	docker compose -p $(PROJECT)-prod -f $(COMPOSE_FILE_PROD) build

prod-up:
	@echo "== Starting prod stack =="
	@echo "NGINX_PORT=$(NGINX_PORT)"
	@echo "CORS_ORIGINS=$(CORS_ORIGINS)"
	NGINX_PORT=$(NGINX_PORT) \
	CORS_ORIGINS="$(CORS_ORIGINS)" \
	docker compose -p $(PROJECT)-prod -f $(COMPOSE_FILE_PROD) up -d --build
	@echo "Open: http://localhost:$(NGINX_PORT)/"

prod-down:
	docker compose -p $(PROJECT)-prod -f $(COMPOSE_FILE_PROD) down -v

logs:
	docker compose -p $(PROJECT)-dev -f $(COMPOSE_FILE_DEV) logs -f --tail=200 || \
	docker compose -p $(PROJECT)-prod -f $(COMPOSE_FILE_PROD) logs -f --tail=200

ps:
	docker compose -p $(PROJECT)-dev -f $(COMPOSE_FILE_DEV) ps || \
	docker compose -p $(PROJECT)-prod -f $(COMPOSE_FILE_PROD) ps

doctor:
	@echo "== Docker version ==" && docker version || true
	@echo "== Docker contexts ==" && docker context ls || true
	@echo "== Using dev compose file ==" && /bin/echo "$(COMPOSE_FILE_DEV)"
	@echo "== Using prod compose file ==" && /bin/echo "$(COMPOSE_FILE_PROD)"
	@echo "== docker ps test ==" && docker ps || true
