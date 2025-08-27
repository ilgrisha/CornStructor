# Version: v4.0.0

# CornStructor

CornStructor is a modern, containerized full-stack app:

- **Frontend**: Angular (dev with `ng serve`, prod as static assets)
- **Backend**: FastAPI (OpenAPI docs at `/docs`)
- **Reverse Proxy (prod)**: Nginx
- **Containers**: Docker Compose (`dev` for live reload, `prod` for immutable builds)

This README tells you how to run **development**, **production**, test the **frontend**, test the **backend**, and exercise the **API**. A quick **troubleshooting** section is included at the end.

---

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Environment Setup](#environment-setup)
3. [Development (live reload)](#development-live-reload)
4. [Production (immutable images)](#production-immutable-images)
5. [Project Structure](#project-structure)
6. [Backend API](#backend-api)
   - [OpenAPI & Docs](#openapi--docs)
   - [Health Check](#health-check)
   - [Start a Design Job](#start-a-design-job)
   - [Stream Logs (SSE)](#stream-logs-sse)
7. [Testing](#testing)
   - [Backend tests (pytest)](#backend-tests-pytest)
   - [Frontend tests (Angular + Jasmine/Karma)](#frontend-tests-angular--jasminekarma)
8. [Common Commands](#common-commands)
9. [Troubleshooting](#troubleshooting)
10. [Design Decisions](#design-decisions)
11. [Next Steps](#next-steps)

---

## Prerequisites

- **Docker Desktop** (macOS/Windows) or **Docker Engine** (Linux)
- **Make** (optional, but the `Makefile` simplifies commands)
- (Optional for running frontend tests locally) **Node.js 20+** and a local Chrome browser

> If you prefer Colima on macOS: `brew install colima && colima start && docker context use colima`.

---

## Environment Setup

We keep a dedicated compose env file for port substitutions:

```bash
cp .compose.env.example .compose.env
```

Environment variables you can set there:

- **BACKEND_PORT** (default 8000) — dev FastAPI port
- **FRONTEND_PORT_DEV** (default 4200) — dev Angular port
- **NGINX_PORT** (default 8080) — prod Nginx port
- **CORS_ORIGINS** (default `*`) — restrict in prod, e.g. `https://your.domain`

> Note: The backend does not read .env files inside the container. Compose injects the vars it needs.

---

## Development (live reload)

Dev stack uses bind mounts and hot reloads:

```bash
make dev
# Frontend: http://localhost:4200
# Backend docs: http://localhost:8000/docs
# Backend API: http://localhost:8000/api
```

Frontend: `ng serve` runs inside the frontend container with a proxy:

- `/api` → backend:8000
- `/reports` → backend:8000

Backend: `uvicorn --reload` inside the backend container

Artifacts: generated files served at `/reports/<jobId>/...` (dev served by FastAPI)

Stop dev and remove volumes:

```bash
make dev-down
```

---

## Production (immutable images)

Build and run optimized images with Nginx reverse proxy:

```bash
make prod-up
# Open: http://localhost:8080
```

- Nginx serves built Angular assets and proxies `/api` → backend:8000
- `/reports` is served statically from a shared volume by Nginx

Stop & remove:

```bash
make prod-down
```

Rebuild only:

```bash
make prod-build
```

---

## Project Structure

```bash
backend/
  app/
    api/
      v1/
        design.py         # POST /design/start, GET /design/{jobId}/logs (SSE)
        health.py         # GET /health
    core/
      config.py           # FastAPI settings (Pydantic Settings)
      jobs/
        runner.py         # Async job runner + SSE queues
    main.py               # FastAPI app factory, routers, static /reports in dev
  requirements.txt
  Dockerfile

frontend/
  src/
    app/
      core/services/
        design.service.ts         # API wrapper for start + SSE logs
        design.service.spec.ts    # Unit tests
    environments/
      environment.ts
      environment.development.ts
  proxy.conf.json
  Dockerfile

nginx/
  nginx.conf            # prod reverse proxy + static reports

docker-compose.dev.yml  # dev stack (bind mounts + reload)
docker-compose.yml      # prod stack (immutable)
.compose.env.example     # sample env for compose
Makefile                # convenience commands
```

---

## Backend API

### OpenAPI & Docs

- Dev: http://localhost:8000/docs (schema at `/api/openapi.json`)
- Prod: http://localhost:8080/docs (proxied through Nginx)

### Health Check

```bash
curl -s http://localhost:8000/api/health   # dev
# or prod:
curl -s http://localhost:8080/api/health
```

Response:

```json
{"status":"ok"}
```

### Start a Design Job

`POST /api/design/start` (returns both job_id and jobId for compatibility)

```bash
curl -s -X POST http://localhost:8000/api/design/start   -H "Content-Type: application/json"   -d '{ "sequence": "ACGTACGTACGTACGTACGT" }'
# => {"job_id":"<id>","jobId":"<id>"}
```

### Stream Logs (SSE)

Use the jobId from the start response.

Dev (FastAPI directly):

```bash
curl -N http://localhost:8000/api/design/<jobId>/logs
```

Prod (Nginx):

```bash
curl -N http://localhost:8080/api/design/<jobId>/logs
```

Lines may include:

- `RUN: ...` (invoked CLI)
- normal logs and `ERR: ...` lines
- `RESULT: /reports/<jobId>/<file or directory>` artifact links
- `EXIT: <code>` (log stream completed sentinel)

---

## Testing

### Backend tests (pytest)

Run in the backend container (dev stack):

```bash
make dev
docker compose -f docker-compose.dev.yml exec backend pytest -q
```

Or target a specific test:

```bash
docker compose -f docker-compose.dev.yml exec backend pytest -q backend/tests/test_api_health.py
```

Local run (outside Docker) is fine too if you have Python 3.11 and `pip install -r backend/requirements.txt`.

### Frontend tests (Angular + Jasmine/Karma)

Recommended: run on your host (uses your local Chrome):

```bash
cd frontend
npm install
npm test
# or CI-style:
npm run test:ci
```

If you prefer running inside the dev container (headless only, additional deps may be required):

```bash
# Install test deps inside the container (first time):
docker compose -f docker-compose.dev.yml exec frontend sh -lc   "npm install -D jasmine-core @types/jasmine karma karma-chrome-launcher karma-jasmine karma-jasmine-html-reporter karma-coverage @angular-devkit/build-angular"

# Then:
docker compose -f docker-compose.dev.yml exec frontend sh -lc "npm test"
```

If Chrome is missing in the container, run tests on host or add a Chrome binary to the image.

---

## Common Commands

```bash
# Dev up / down
make dev
make dev-down

# Prod up / down / build
make prod-up
make prod-down
make prod-build

# Tail logs from running services
make logs

# Doctor (basic Docker diagnostics)
make doctor
```

---

## Troubleshooting

### Docker daemon not running

- Start Docker Desktop: `open -a "Docker"`
- Verify: `docker version`

### Makefile “missing separator”

- Recipe lines must start with a **TAB**. Use the provided Makefile.

### Frontend dev proxy errors (ECONNREFUSED)

- Usually means backend isn’t up. Check `backend-1` logs.
- Dev URLs:
  - API: http://localhost:8000/api
  - Frontend: http://localhost:4200

### npm ci fails due to lock mismatch

- In dev, we use `npm install` inside the container to reconcile.
- To restore strictness later, clean up lock locally and switch back to `npm ci`.

### Backend Pydantic “extra inputs are not permitted”

- We configured `config.py` to ignore unknown env keys and not read `.env` in-container.
- Ensure you’re using `.compose.env` for Compose substitution.

### SSE not updating in UI

- Check `CORS_ORIGINS` and make sure the frontend uses EventSource against `/api/design/<jobId>/logs`.
- In prod, confirm Nginx `location /api/` has `proxy_buffering off;`.

### Clean rebuilds

```bash
make dev-down
docker compose --env-file .compose.env -f docker-compose.dev.yml build --no-cache backend
make dev
```

---

## Design Decisions

- **FastAPI + SSE**: The API spawns the existing CLI via `asyncio.create_subprocess_exec` and streams combined stdout/stderr as Server-Sent Events for simple, reliable progress in the UI.
- **In-memory job manager**: Adequate for a single instance. Swap for Redis/Kafka when scaling horizontally.
- **Shared reports volume**: Artifacts written by backend are immediately web-accessible (served by FastAPI in dev, Nginx in prod).
- **Dev vs Prod**:
  - Dev uses bind mounts and hot reload for rapid iteration.
  - Prod builds immutable images (Angular build → Nginx static; FastAPI served by Uvicorn) for consistency and performance.
- **Environment handling**: Backend ignores unrelated env vars and doesn’t auto-load `.env` to avoid collisions with Compose substitution files.
