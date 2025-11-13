# CornStructor Frontend (Angular 18+)

Standalone Angular 18 app (signals, no NgModules). Card-based UI in a modern Apple-like style.

## Quick start

```bash
cd frontend
npm i
npm start
```

The app expects a backend with:
- `POST /api/design/start` -> `{ jobId }`
- `GET /api/design/:jobId/logs` -> SSE stream (text/event-stream)
- `GET /api/design/:jobId/result` -> `{ ok, done, jobId, outputDir, links: [...], fastaOut? }`
