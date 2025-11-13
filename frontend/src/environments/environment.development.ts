// File: frontend/src/environments/environment.development.ts
// Version: v0.1.0
/**
 * Development environment.
 * ng serve + proxy.conf.json forwards /api to backend container.
 */
export const environment = {
  production: false,
  apiBase: '/api',
  sseBase: '/api',
  reportsBase: '/reports'
};
