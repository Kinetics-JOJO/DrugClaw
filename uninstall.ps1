param(
  [string]$InstallDir,
  [switch]$CleanPath
)

$ErrorActionPreference = 'Stop'
$BinName = 'drugclaw.exe'

function Get-EnvAlias {
  param(
    [string]$Primary,
    [string]$Legacy,
    [string]$Fallback = ''
  )

  $primaryValue = (Get-Item -Path "Env:$Primary" -ErrorAction SilentlyContinue).Value
  if ($primaryValue) {
    return $primaryValue
  }
  $legacyValue = (Get-Item -Path "Env:$Legacy" -ErrorAction SilentlyContinue).Value
  if ($legacyValue) {
    return $legacyValue
  }
  return $Fallback
}

if (-not $PSBoundParameters.ContainsKey('InstallDir') -or [string]::IsNullOrWhiteSpace($InstallDir)) {
  $installFallback = if ($env:USERPROFILE) {
    Join-Path $env:USERPROFILE '.local\bin'
  } else {
    Join-Path (Get-Location) '.local\bin'
  }
  $InstallDir = Get-EnvAlias -Primary 'DRUGCLAW_INSTALL_DIR' -Legacy 'MICROCLAW_INSTALL_DIR' -Fallback $installFallback
}

function Write-Info([string]$msg) {
  Write-Host $msg
}

function Path-Contains([string]$pathValue, [string]$dir) {
  if ([string]::IsNullOrWhiteSpace($pathValue)) { return $false }
  $needle = $dir.Trim().TrimEnd('\\').ToLowerInvariant()
  foreach ($part in $pathValue.Split(';')) {
    if ([string]::IsNullOrWhiteSpace($part)) { continue }
    if ($part.Trim().TrimEnd('\\').ToLowerInvariant() -eq $needle) {
      return $true
    }
  }
  return $false
}

function Remove-UserPathEntry([string]$dir) {
  $userPath = [Environment]::GetEnvironmentVariable('Path', 'User')
  if ([string]::IsNullOrWhiteSpace($userPath)) { return $false }

  $needle = $dir.Trim().TrimEnd('\\').ToLowerInvariant()
  $parts = @()
  $changed = $false

  foreach ($part in $userPath.Split(';')) {
    if ([string]::IsNullOrWhiteSpace($part)) { continue }
    $normalized = $part.Trim().TrimEnd('\\').ToLowerInvariant()
    if ($normalized -eq $needle) {
      $changed = $true
      continue
    }
    $parts += $part.Trim()
  }

  if (-not $changed) { return $false }

  $newPath = ($parts -join ';')
  [Environment]::SetEnvironmentVariable('Path', $newPath, 'User')

  if (Path-Contains $env:Path $dir) {
    $procParts = @()
    foreach ($part in $env:Path.Split(';')) {
      if ([string]::IsNullOrWhiteSpace($part)) { continue }
      if ($part.Trim().TrimEnd('\\').ToLowerInvariant() -eq $needle) {
        continue
      }
      $procParts += $part.Trim()
    }
    $env:Path = ($procParts -join ';')
  }

  return $true
}

function Resolve-Targets([string]$installDir, [string]$binName) {
  $targets = [System.Collections.Generic.List[string]]::new()

  if (-not [string]::IsNullOrWhiteSpace($installDir)) {
    $targets.Add((Join-Path $installDir $binName))
  }

  $cmd = Get-Command drugclaw -ErrorAction SilentlyContinue
  if ($cmd -and -not [string]::IsNullOrWhiteSpace($cmd.Source)) {
    $targets.Add($cmd.Source)
  }

  $targets.Add((Join-Path $env:USERPROFILE '.local\bin\drugclaw.exe'))

  $seen = @{}
  foreach ($target in $targets) {
    if (-not [string]::IsNullOrWhiteSpace($target) -and -not $seen.ContainsKey($target)) {
      $seen[$target] = $true
      $target
    }
  }
}

$removed = 0
$failed = $false

Write-Info "Uninstalling drugclaw..."
foreach ($target in Resolve-Targets -installDir $InstallDir -binName $BinName) {
  if (Test-Path -LiteralPath $target) {
    try {
      Remove-Item -LiteralPath $target -Force
      Write-Info "Removed: $target"
      $removed++
    } catch {
      Write-Info "Failed to remove: $target"
      Write-Info "Reason: $($_.Exception.Message)"
      $failed = $true
    }
  }
}

if ($CleanPath) {
  if (Remove-UserPathEntry $InstallDir) {
    Write-Info "Removed '$InstallDir' from user PATH."
  } else {
    Write-Info "User PATH did not include '$InstallDir'."
  }
}

if ($failed) {
  exit 1
}

if ($removed -eq 0) {
  Write-Info "drugclaw binary not found. Nothing to uninstall."
  exit 0
}

Write-Info ""
Write-Info "drugclaw has been removed."
Write-Info "Optional cleanup (not removed automatically):"
Write-Info "  Remove-Item -Recurse -Force $HOME\\.drugclaw\\runtime"
Write-Info "  Remove-Item -Force .\\drugclaw.config.yaml,.\\drugclaw.config.yml"
