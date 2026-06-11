@echo off
:: Enable WSL2 features - run as Administrator
echo Enabling Windows Subsystem for Linux...
dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
echo.
echo Enabling Virtual Machine Platform...
dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
echo.
echo Done! Please reboot now.
echo After reboot, tell Hermes to continue the Salmonella assembly.
pause
