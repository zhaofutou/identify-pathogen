#!/bin/bash
set -e

case "${1:-web}" in
    web)
        echo "=== BioLab Web App ==="
        echo "Starting Flask on port 5050..."
        echo "Access at: http://localhost:5050"
        exec python3 /app/webapp/app.py
        ;;
    shell)
        echo "=== BioLab Debug Shell ==="
        exec /bin/bash
        ;;
    test)
        echo "=== Running Pipeline Test ==="
        exec python3 /app/pathogen_panel.py
        ;;
    *)
        echo "Usage: docker run biolab [web|shell|test]"
        exit 1
        ;;
esac
