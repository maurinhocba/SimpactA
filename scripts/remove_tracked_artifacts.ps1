# PowerShell version: remove_tracked_artifacts.ps1
# Run from repository root. Removes tracked build artifacts (cached) and commits the deletions.
# WARNING: This will remove the matching files from the repository in the current branch.

param()

$branch = git rev-parse --abbrev-ref HEAD
if ($branch -ne 'try/borrarbins') {
    Write-Host "Warning: you are on branch $branch. It's recommended to run this script on branch try/borrarbins. Continue? (y/N)"
    $ans = Read-Host
    if ($ans -ne 'y') { Write-Host 'Aborted.'; exit 1 }
}

$patterns = @('*.mod','*.o','*.obj','.vs','*.suo','*.user','build','bin','*.log')
$matches = @()

foreach ($p in $patterns) {
    $list = git ls-files -- "$p" | ForEach-Object { $_.Trim() } | Where-Object { $_ -ne '' }
    if ($list) { $matches += $list }
}

if ($matches.Count -eq 0) { Write-Host 'No tracked files matched the patterns. Nothing to do.'; exit 0 }

Write-Host "The following tracked paths will be removed from the repo index (kept on disk):"
$matches | ForEach-Object { Write-Host $_ }

$confirm = Read-Host "Proceed to remove these from the index and commit? (y/N)"
if ($confirm -ne 'y') { Write-Host 'Aborted.'; exit 1 }

foreach ($f in $matches) {
    git rm --cached -r --ignore-unmatch -- "$f" | Out-Null
}

git add .gitignore

git commit -m "Remove tracked build artifacts and add .gitignore"
Write-Host "Committed removals. Push the branch with: git push origin try/borrarbins"
