@echo off
echo ============================================================
echo Starting BioLab Service (Docker Edition)
echo ============================================================
echo.
echo Launching WSL2 keep-alive process...
start /b wsl -d Ubuntu -u root bash -c "while true; do sleep 3600; done"
echo.
echo Starting docker compose...
echo Access the application at: http://localhost:5050
echo.
wsl -d Ubuntu docker compose -f /mnt/c/Users/rdpuser/bacterial_assemble/docker-compose.yml up
echo.
echo Stopping BioLab Service...
wsl -d Ubuntu docker compose -f /mnt/c/Users/rdpuser/bacterial_assemble/docker-compose.yml down
echo Done!
pause
