#!/bin/bash
set -e

psql -v ON_ERROR_STOP=1 -U "$POSTGRES_USER" -d "$POSTGRES_DB" -c "ALTER DATABASE ${POSTGRES_DB} SET citus.multi_shard_modify_mode = 'sequential';"