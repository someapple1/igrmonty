cmake --build build-riaf-release --target grmonty
.\build-riaf-release\grmonty.exe -par model\riaf\example.par
python tools\plspec.py spectrum.h5
for ($i = 1; $i -le 3; $i++) {
  $dir = "runs\test_$i"
  New-Item -ItemType Directory -Force $dir | Out-Null

  .\build-riaf-release\grmonty.exe -par model\riaf\example.par

  Copy-Item spectrum.h5 "$dir\spectrum.h5" -Force
  python tools\plspec.py spectrum.h5
  Copy-Item spectrum.png "$dir\spectrum.png" -Force
  Copy-Item spectrum-avg.png "$dir\spectrum-avg.png" -Force

  Write-Host "Finished run $i"
}