# File: Makefile
# Version: v0.6.0
.PHONY: help dev dev-down prod-build prod-up prod-down logs doctor

help:
	@echo "Targets:"
	@echo "  dev        - start live-reload dev stack"
	@echo "  dev-down   - stop dev stack"
	@echo "  prod-build - build production images"
	@echo "  prod-up    - start production stack"
	@echo "  prod-down  - stop production stack"
	@echo "  logs       - tail logs from all services"
	@echo "  doctor     - diagnose Docker/compose environment"

dev:
	docker compose --env-file .compose.env -f docker-compose.dev.yml up --build

dev-down:
	docker compose --env-file .compose.env -f docker-compose.dev.yml down -v

prod-build:
	docker compose --env-file .compose.env -f docker-compose.yml build

prod-up:
	docker compose --env-file .compose.env -f docker-compose.yml up -d --build

prod-down:
	docker compose --env-file .compose.env -f docker-compose.yml down -v

logs:
	docker compose logs -f --tail=200

doctor:
	@echo "== Docker version ==" && docker version || true
	@echo "== Docker contexts ==" && docker context ls || true
	@echo "== Using env-file ==" && /bin/echo ".compose.env"
	@echo "== docker ps test ==" && docker ps || true
