# Publishing checklist

Order matters: OSF first (the pre-registration's whole function is proving
the predictions predate the data), then GitHub, then Zenodo via the GitHub
integration.

## 1. OSF pre-registration

Upload `preregistration/fanning_preregistration_v1_0.md` as a registration.
Record the resulting DOI.

## 2. GitHub

    gh repo create genealogical-topology-pilot --public \
      --description "Verified catalogue, analysis code, and pre-registration for a pilot test of genealogical topology in cosmological anomalies"
    git init && git add -A
    git commit -m "v1.0.0: verified catalogue, consolidated analysis, pre-registration"
    git branch -M main
    git remote add origin https://github.com/<user>/genealogical-topology-pilot.git
    git push -u origin main
    git tag -a v1.0.0 -m "Verified catalogue release"
    git push origin v1.0.0

Topics to add: cosmology, cosmological-principle, anisotropy, reproducible-research,
open-science, preregistration, philosophy-of-science, ai-collaboration

## 3. Zenodo

Enable the Zenodo GitHub integration for this repo, then cut a GitHub release
from the v1.0.0 tag. Zenodo mints a DOI automatically. Upload type: Software.

## 4. Fill in the placeholders

Replace `[DOI pending]` with the Zenodo concept DOI (the one that always
resolves to the latest version) and the OSF DOI in:

  - README.md
  - CITATION.cff (repository field)
  - cosmological_companion_paper.md (header, Section 9, Data and Code Availability)
  - fanning_preregistration_v1_0.md (header)
  - Nothing_Described_v4.md (Section 36 pre-registration reference)
