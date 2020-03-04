taskkill /F /IM SimExplorer.exe
taskkill /F /IM SimManager.exe
rem taskkill /F /IM ModelServer.exe
Start ""/min "C:\Program Files (x86)\Kongsberg\K-Spice\bin\SimManager.exe" -project BasicProcess -Minimized
Start ""/b "C:\Program Files (x86)\Kongsberg\K-Spice\bin\SimExplorer.exe" -d -search_path Operator -app Operator -instance 2