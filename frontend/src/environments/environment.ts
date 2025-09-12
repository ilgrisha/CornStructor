// File: frontend/src/environments/environment.ts
// Version: v0.2.0
/**
 * Production environment.
 * In production behind Nginx, API is reverse-proxied under /api.
 */
export const environment = {
  production: true,
  apiBase: '/api',
  sseBase: '/api',
  reportsBase: '/reports'
};
