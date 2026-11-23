param(
    [string]$Name = 'scaffold',
    [string]$Prefix = '',
    [switch]$SkipMoleculeDownload
)
$ErrorActionPreference = 'Stop'
$RepoRoot = Split-Path -Parent $PSScriptRoot
if ($PSBoundParameters.ContainsKey('Name') -and $Prefix) {
    throw 'Use -Name or -Prefix, not both.'
}
if ($Prefix) {
    $Prefix = [IO.Path]::GetFullPath($Prefix)
    $PythonExe = Join-Path $Prefix 'python.exe'
    if (-not (Test-Path -LiteralPath $PythonExe)) {
        conda env create --prefix $Prefix --file (Join-Path $RepoRoot 'environment.yml') --yes
        if ($LASTEXITCODE -ne 0) { throw 'Conda environment creation failed' }
    }
    $ActivationTarget = $Prefix
} else {
    if ($Name -notmatch '^[A-Za-z0-9][A-Za-z0-9_.-]*$') { throw 'Invalid environment name.' }
    $CondaInfoText = conda info --json
    if ($LASTEXITCODE -ne 0) { throw 'Cannot read Conda configuration' }
    $CondaInfo = $CondaInfoText | ConvertFrom-Json
    $PythonExe = $null
    foreach ($EnvironmentDirectory in $CondaInfo.envs_dirs) {
        $Candidate = Join-Path (Join-Path $EnvironmentDirectory $Name) 'python.exe'
        if (Test-Path -LiteralPath $Candidate) { $PythonExe = $Candidate; break }
    }
    if (-not $PythonExe) {
        conda env create --name $Name --file (Join-Path $RepoRoot 'environment.yml') --yes
        if ($LASTEXITCODE -ne 0) { throw 'Conda environment creation failed' }
        foreach ($EnvironmentDirectory in $CondaInfo.envs_dirs) {
            $Candidate = Join-Path (Join-Path $EnvironmentDirectory $Name) 'python.exe'
            if (Test-Path -LiteralPath $Candidate) { $PythonExe = $Candidate; break }
        }
    }
    if (-not $PythonExe) { throw "Cannot locate Python in named Conda environment $Name" }
    $ActivationTarget = $Name
}
& $PythonExe -m pip install 'torch==2.6.0' --index-url https://download.pytorch.org/whl/cpu
if ($LASTEXITCODE -ne 0) { throw 'CPU PyTorch installation failed' }
# The lock includes the tested parser/annotation dependency closure.
& $PythonExe -m pip install --no-deps -r (Join-Path $RepoRoot 'requirements-win-cpu.lock')
if ($LASTEXITCODE -ne 0) { throw 'Locked dependency installation failed' }
& $PythonExe -m pip install --no-deps 'boltzgen==0.3.2'
if ($LASTEXITCODE -ne 0) { throw 'BoltzGen parser installation failed' }
& $PythonExe -m pip install --no-deps -e $RepoRoot
if ($LASTEXITCODE -ne 0) { throw 'Converter installation failed' }
if (-not $SkipMoleculeDownload) {
    & $PythonExe (Join-Path $PSScriptRoot 'fetch_boltzgen_mols.py')
    if ($LASTEXITCODE -ne 0) { throw 'Molecule dictionary download failed' }
}
& $PythonExe -c "import scaffold, anarcii, gemmi; from boltzgen.data.parse.schema import YamlDesignParser; print('Annotation and official parser imports OK')"
if ($LASTEXITCODE -ne 0) { throw 'Import verification failed' }
Write-Host "Ready. conda activate $ActivationTarget"
Write-Host 'Run: scaffold doctor'
