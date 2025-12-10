# File: Makefile
# Version: v0.8.0
.PHONY: help dev dev-down prod-build prod-up prod-down logs doctor ps \
	up down restart build build-nocache rebuild logs-backend logs-frontend \
	backend-shell frontend-shell db-shell print-ports db-info db-tables \
	db-describe migrate migrate-revision remove-db lint-backend \
	lint-backend-fix fmt-backend typecheck-backend test-backend lint-frontend \
	fmt-frontend test-frontend

MODE ?= dev
HOST_IP ?= 0.0.0.0
BACKEND_PORT ?= 8016
FRONTEND_PORT_DEV ?= 4216
NGINX_PORT ?= 8098
DB_PORT ?= 5432
DB_USER ?= appuser
DB_NAME ?= cornstructor-db
DB_VOLUME ?= cornstructor_dbdata
CORS_ORIGINS ?= *
COMPOSE_FILE_DEV := docker-compose.dev.yml
COMPOSE_FILE_PROD := docker-compose.yml
# Optional: allow multiple concurrent stacks by project name
PROJECT ?= cornstructor

DEV_ENV := HOST_IP=$(HOST_IP) BACKEND_PORT=$(BACKEND_PORT) FRONTEND_PORT_DEV=$(FRONTEND_PORT_DEV) CORS_ORIGINS="$(CORS_ORIGINS)"
PROD_ENV := HOST_IP=$(HOST_IP) NGINX_PORT=$(NGINX_PORT) CORS_ORIGINS="$(CORS_ORIGINS)"
DEV_COMPOSE := $(DEV_ENV) docker compose -p $(PROJECT)-dev -f $(COMPOSE_FILE_DEV)
PROD_COMPOSE := $(PROD_ENV) docker compose -p $(PROJECT)-prod -f $(COMPOSE_FILE_PROD)

ifeq ($(MODE),dev)
  COMPOSE_FILE := $(COMPOSE_FILE_DEV)
  COMPOSE_PROJECT := $(PROJECT)-dev
  COMPOSE_ENV := $(DEV_ENV)
else ifeq ($(MODE),prod)
  COMPOSE_FILE := $(COMPOSE_FILE_PROD)
  COMPOSE_PROJECT := $(PROJECT)-prod
  COMPOSE_ENV := $(PROD_ENV)
else
  $(error MODE must be either dev or prod)
endif

COMPOSE := $(COMPOSE_ENV) docker compose -p $(COMPOSE_PROJECT) -f $(COMPOSE_FILE)

help:
	@echo "Targets:"
	@echo "  dev / dev-down             - start/stop live-reload dev stack (foreground)"
	@echo "  prod-build / prod-up / prod-down - manage production stack"
	@echo "  up/down/restart/build*     - generic compose controls (use MODE=dev|prod)"
	@echo "  logs*/ps                   - inspect services for the selected MODE"
	@echo "  backend-shell/frontend-shell- open shells inside running containers"
	@echo "  db-* / migrate*             - database + Alembic helpers (require db service)"
	@echo "  lint-*/fmt-*/test-*         - quality commands (containerized)"
	@echo "  doctor                      - diagnose Docker/compose environment"
	@echo ""
	@echo "Port overrides (examples):"
	@echo "  make dev BACKEND_PORT=9001 FRONTEND_PORT_DEV=4301"
	@echo "  make prod-up NGINX_PORT=9090"
	@echo "  make up MODE=prod NGINX_PORT=9090"
	@echo "  make up MODE=dev BACKEND_PORT=9001 FRONTEND_PORT_DEV=4301"
	@echo ""
	@echo "Use print-ports to review the active values."

dev:
	@echo "== Starting dev stack =="
	@echo "BACKEND_PORT=$(BACKEND_PORT)"
	@echo "FRONTEND_PORT_DEV=$(FRONTEND_PORT_DEV)"
	@echo "CORS_ORIGINS=$(CORS_ORIGINS)"
	$(DEV_COMPOSE) up --build

dev-down:
	$(DEV_COMPOSE) down -v

prod-build:
	$(PROD_COMPOSE) build

prod-up:
	@echo "== Starting prod stack =="
	@echo "NGINX_PORT=$(NGINX_PORT)"
	@echo "CORS_ORIGINS=$(CORS_ORIGINS)"
	$(PROD_COMPOSE) up -d --build
	@echo "Open: http://localhost:$(NGINX_PORT)/"

prod-down:
	$(PROD_COMPOSE) down -v

up:
	@echo "== docker compose up ($(MODE)) =="
	$(COMPOSE) up -d

down:
	$(COMPOSE) down

restart:
	$(COMPOSE) restart

build:
	$(COMPOSE) build

build-nocache:
	$(COMPOSE) build --no-cache

rebuild: down build-nocache up

ps:
	$(COMPOSE) ps

logs:
	$(COMPOSE) logs -f --tail=200

logs-backend:
	$(COMPOSE) logs -f backend

logs-frontend:
	$(COMPOSE) logs -f frontend

backend-shell:
	$(COMPOSE) exec backend /bin/bash || \
	$(COMPOSE) exec backend /bin/sh

frontend-shell:
	$(COMPOSE) exec frontend /bin/bash || \
	$(COMPOSE) exec frontend /bin/sh

db-shell:
	@if ! $(COMPOSE) config --services | grep -q '^db$$'; then \
	  echo "db service is not defined in $(COMPOSE_FILE)."; \
	  exit 0; \
	fi
	$(COMPOSE) exec db psql -U $(DB_USER) -d $(DB_NAME)

print-ports:
	@echo "MODE             = $(MODE)"
	@echo "HOST_IP          = $(HOST_IP)"
	@echo "NGINX_PORT       = $(NGINX_PORT)"
	@echo "FRONTEND_PORT_DEV= $(FRONTEND_PORT_DEV)"
	@echo "BACKEND_PORT     = $(BACKEND_PORT)"
	@echo "DB_PORT          = $(DB_PORT)"
	@echo "DB_NAME          = $(DB_NAME)"
	@echo "DB_USER          = $(DB_USER)"
	@echo "DB_VOLUME        = $(DB_VOLUME)"
	@echo "CORS_ORIGINS     = $(CORS_ORIGINS)"

db-info:
	@if ! $(COMPOSE) config --services | grep -q '^db$$'; then \
	  echo "db service is not defined in $(COMPOSE_FILE)."; \
	  exit 0; \
	fi
	$(COMPOSE) exec db sh -c "psql -U $(DB_USER) -d $(DB_NAME) -c 'SELECT current_database() AS db, current_user AS user, inet_server_addr() AS host, inet_server_port() AS port, version();'"

db-tables:
	@if ! $(COMPOSE) config --services | grep -q '^db$$'; then \
	  echo "db service is not defined in $(COMPOSE_FILE)."; \
	  exit 0; \
	fi
	$(COMPOSE) exec db sh -c "psql -U $(DB_USER) -d $(DB_NAME) -c '\\dt'"

db-describe:
	@if ! $(COMPOSE) config --services | grep -q '^db$$'; then \
	  echo "db service is not defined in $(COMPOSE_FILE)."; \
	  exit 0; \
	fi
	$(COMPOSE) exec db sh -c "psql -U $(DB_USER) -d $(DB_NAME) -c '\\dt'"
	@read -p "Table name to describe: " tbl; \
	if [ -n "$$tbl" ]; then \
	  $(COMPOSE) exec db sh -c "psql -U $(DB_USER) -d $(DB_NAME) -c '\\d $$tbl'"; \
	else \
	  echo "No table selected."; \
	fi

migrate:
	$(COMPOSE) exec backend alembic upgrade head

migrate-revision:
	@if [ -z "$(MESSAGE)" ]; then \
	  echo "Usage: make migrate-revision MESSAGE=\"your message\""; \
	  exit 1; \
	fi
	$(COMPOSE) exec backend alembic revision --autogenerate -m "$(MESSAGE)"

remove-db:
	@echo "This will remove the Postgres volume '$(DB_VOLUME)' and ALL data."; \
	read -p "Are you sure? [y/N] " ans; \
	if [ "$$ans" = "y" ] || [ "$$ans" = "Y" ]; then \
	  $(COMPOSE) down; \
	  docker volume rm $(DB_VOLUME) || true; \
	else \
	  echo "Aborted."; \
	fi

lint-backend:
	$(COMPOSE) exec backend sh -c "cd /app && ruff check app"

lint-backend-fix:
	$(COMPOSE) exec backend sh -c "cd /app && ruff check app --fix"

fmt-backend:
	$(COMPOSE) exec backend sh -c 'cd /app && black app'

typecheck-backend:
	$(COMPOSE) exec backend sh -c "cd /app && mypy app"

test-backend:
	$(COMPOSE) exec backend sh -c "cd /app && pytest"

lint-frontend:
	$(COMPOSE) exec frontend npm run lint -- --quiet || \
	$(COMPOSE) exec frontend npm run lint

fmt-frontend:
	$(COMPOSE) exec frontend npm run format

test-frontend:
	$(COMPOSE) exec frontend npm test -- --watch=false --browsers=ChromeHeadless || true

doctor:
	@echo "== Docker version ==" && docker version || true
	@echo "== Docker contexts ==" && docker context ls || true
	@echo "== Using dev compose file ==" && /bin/echo "$(COMPOSE_FILE_DEV)"
	@echo "== Using prod compose file ==" && /bin/echo "$(COMPOSE_FILE_PROD)"
	@echo "== docker ps test ==" && docker ps || true
