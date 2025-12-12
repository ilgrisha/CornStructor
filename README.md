# CornStructor
**Version: v4.1.0**

CornStructor is a modern, containerized full‑stack app for hierarchical DNA design and visualization. It helps you design construction trees, inspect oligos/fragments, track past runs (with notes), and re‑generate reports on demand — even after restarts.

- **Frontend**: Angular (dev via `ng serve`, prod as static assets)
- **Backend**: FastAPI (OpenAPI docs at `/docs`)
- **Reverse Proxy (prod)**: Nginx
- **Containers**: Docker Compose (`dev` for hot reload, `prod` for immutable images)
- **Artifacts**: Served under `/reports/<jobId>/...`

This README shows how to run **development** and **production**, pass **custom ports** using `make`, poke the **API**, run **tests**, and troubleshoot common issues.

---

## Overview

- **Construction tree design**: Launch jobs that build hierarchical assemblies with per‑level constraints and global parameters.
- **Live feedback**: Stream pipeline logs over SSE; “Open report” becomes available as soon as artifacts exist.
- **Visualization**: Sequence box overlays fragments/oligos (sense/antisense), shows overlaps/Tm, and keeps overlays in sync with the reference sequence.
- **Run history**: Browse previous runs (with stored notes/descriptions), reload a design into the Sequence box, and open the saved report in a new tab.
- **On‑the‑fly reports**: Dynamic `/reports/<jobId>/…` renderers pull tree/GA/params from the database so reports survive container restarts.

---

## Key Features & Flows

- **Start a design**: Provide sequence + parameters; backend persists Run/Design, merges tree configs, and runs the in‑process pipeline. Logs stream back; artifacts land under `/reports/<jobId>/`.
- **Load from history**: Selecting a completed run restores its reference sequence + tree into the Sequence box, updates the run description field, and keeps the report link targeting that run.
- **Sequence box tooling**:
  - Upload FASTA or paste sequence; tracks reference metadata (FASTA name/id or run note/id).
  - Toggle overlays by design level and mode (fragments/oligos), with sense/antisense coloring and per‑fragment details (length, coords, overlaps, Tm).
  - Range selection, per‑line rulers/ticks, and feature highlighting (GC bands, entropy, repeats, stems, etc.).
- **Persisted context**:
  - Runs store sequence length, params snapshot, note/description, status, exit code, report URL, and linked Design id.
  - Designs store sequence, params JSON (including locked globals/levels), tree JSON, GA progress JSON, and optional name.
- **Reports**: GA progress, tree, analysis, FASTA/CSV exports, and bundle download; all accessible via `/reports/<jobId>/index.html`.

---

## Table of Contents
1. [Prerequisites](#prerequisites)
2. [Environment Setup](#environment-setup)
3. [Quick Start](#quick-start)
4. [Development (live reload)](#development-live-reload)
5. [Production (immutable images)](#production-immutable-images)
6. [Custom Ports with Make](#custom-ports-with-make)
7. [Project Structure](#project-structure)
8. [Backend API](#backend-api)
   - [OpenAPI & Docs](#openapi--docs)
   - [Health Check](#health-check)
   - [Start a Design Job](#start-a-design-job)
   - [Stream Logs (SSE)](#stream-logs-sse)
9. [Testing](#testing)
   - [Backend tests (pytest)](#backend-tests-pytest)
   - [Frontend tests (Angular + Jasmine/Karma)](#frontend-tests-angular--jasminekarma)
10. [Troubleshooting](#troubleshooting)
11. [Design Decisions](#design-decisions)

---

## Prerequisites

- **Docker Desktop** (macOS/Windows) or **Docker Engine** (Linux)
- **Make** (optional but recommended)
- (Optional for local frontend tests) **Node.js 20+** and a local Chrome

> macOS + Colima users: `brew install colima && colima start && docker context use colima`

---

## Environment Setup

We use a Compose env file for convenient overrides:

```bash
cp .compose.env.example .compose.env
```

Environment variables you may set there (or pass via `make`, see below):

- **BACKEND_PORT** (default **8000**) — dev FastAPI port
- **FRONTEND_PORT_DEV** (default **4200**) — dev Angular port
- **NGINX_PORT** (default **8080**) — prod Nginx port
- **CORS_ORIGINS** (default `*`) — set to your domains in prod

> The backend doesn’t read `.env` files itself. Compose injects env vars into containers.

---

## Quick Start

```bash
# Development (hot reload)
make dev
# Frontend: http://localhost:4200
# Backend API/docs: http://localhost:8000/api , http://localhost:8000/docs

# Production (immutable builds)
make prod-up
# App: http://localhost:8080
```

Stop stacks:

```bash
make dev-down     # dev
make prod-down    # prod
```

---

## Development (live reload)

Dev uses bind mounts and hot reload for fast iteration.

```bash
make dev
```

- **Frontend**: Angular `ng serve` with a proxy
  - `/api` → backend:8000
  - `/reports` → backend:8000
- **Backend**: Uvicorn `--reload`
- **Artifacts**: generated under `/reports/<jobId>/...` (served by FastAPI)

Handy commands:

```bash
make logs   # tail logs
make ps     # list services
make dev-down
```

---

## Production (immutable images)

```bash
make prod-up
# Open http://localhost:8080
```

- Nginx serves built Angular assets
- Proxies `/api` → backend:8000
- Serves `/reports` statically from a shared `reports` volume

Other targets:

```bash
make prod-build   # build only
make prod-down    # stop + remove volumes
```

---

## Custom Ports with Make

Run multiple stacks side-by-side or reuse free ports by passing Make variables:

```bash
# Dev: change backend + frontend ports
make dev BACKEND_PORT=9000 FRONTEND_PORT_DEV=4300

# Prod: change Nginx external port
make prod-up NGINX_PORT=9090

# Optional: separate project names for concurrent stacks
make dev PROJECT=csA BACKEND_PORT=9001 FRONTEND_PORT_DEV=4301
make prod-up PROJECT=csB NGINX_PORT=9091
```

`CORS_ORIGINS` can also be provided at invocation time:

```bash
make dev CORS_ORIGINS="http://localhost:4300,http://127.0.0.1:4300"
```

---

## Project Structure

```text
backend/
  app/
    api/
      v1/
        design.py               # POST /design/start, GET /design/{jobId}/logs (SSE)
        health.py               # GET /health
        reports_dynamic.py      # /reports/<jobId>/... dynamic artifacts
    core/
      assembly/                 # GA + planners
      export/                   # exporters: FASTA/CSV/JSON/GenBank
      reporting/                # on-the-fly render + index
      models/                   # FragmentNode, etc.
      config.py                 # FastAPI settings
    cli/
      genbank_from_tree.py      # CLI to generate GenBank from tree.json + FASTA
    main.py                     # FastAPI app (dev serves reports)
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
  nginx.conf                    # prod reverse proxy + static reports

docker-compose.dev.yml          # dev stack (bind mounts + reload)
docker-compose.yml              # prod stack (immutable)
.compose.env.example            # sample env overrides for Compose
Makefile                        # convenience commands
```

---

## Backend API

### OpenAPI & Docs

- Dev: <http://localhost:8000/docs> (schema at `/api/openapi.json`)
- Prod: <http://localhost:8080/docs> (proxied through Nginx)

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

`POST /api/design/start` (returns both `job_id` and `jobId`)

```bash
curl -s -X POST http://localhost:8000/api/design/start   -H "Content-Type: application/json"   -d '{ "sequence": "ACGTACGTACGTACGTACGT" }'
# => {"job_id":"<id>","jobId":"<id>"}
```

### Stream Logs (SSE)

Use the `jobId` from the start response.

Dev:
```bash
curl -N http://localhost:8000/api/design/<jobId>/logs
```

Prod:
```bash
curl -N http://localhost:8080/api/design/<jobId>/logs
```

Log stream includes:
- `RUN:` invoked tools
- progress + `ERR:` lines
- `RESULT: /reports/<jobId>/<file>` artifact links
- `EXIT: <code>` sentinel

---

## Testing

### Backend tests (pytest)

```bash
make dev
docker compose -f docker-compose.dev.yml exec backend pytest -q
# specific test:
docker compose -f docker-compose.dev.yml exec backend pytest -q backend/tests/test_api_health.py
```

Local (outside Docker) works too with Python 3.11 and:
```bash
pip install -r backend/requirements.txt
pytest -q
```

### Frontend tests (Angular + Jasmine/Karma)

Recommended on host:
```bash
cd frontend
npm install
npm test
# or headless CI:
npm run test:ci
```

In container (headless only; may need extra deps):
```bash
docker compose -f docker-compose.dev.yml exec frontend sh -lc "npm test"
```

---

## Troubleshooting

**Docker not running**
- Start Docker Desktop / daemon; verify with `docker version`

**Makefile “missing separator”**
- Recipe lines must start with a TAB

**Frontend proxy errors (ECONNREFUSED)**
- Backend likely down; check logs
- Dev endpoints:
  - API: <http://localhost:8000/api>
  - Frontend: <http://localhost:4200>

**npm lock mismatch**
- Dev uses `npm install` in the container to reconcile quickly
- For strict builds later, fix the lock and switch to `npm ci`

**CORS/SSE issues**
- Set `CORS_ORIGINS` appropriately
- Nginx dev/prod configs set `proxy_buffering off` for SSE

**Clean rebuilds**
```bash
make dev-down
docker compose --env-file .compose.env -f docker-compose.dev.yml build --no-cache backend
make dev
```
