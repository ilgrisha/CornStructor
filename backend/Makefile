# File: backend/Makefile
# Version: v0.1.0
# Handy shortcuts for migrations
.PHONY: db-rev db-up db-dn

db-rev:
	\talembic revision --autogenerate -m "$(m)"

db-up:
	\talembic upgrade head

db-dn:
	\talembic downgrade -1
