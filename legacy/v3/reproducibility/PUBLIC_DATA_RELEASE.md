# Public data release procedure

This document records the author-side publication procedure that converts the preserved Chapter 1 v2 inputs into a credential-free third-party reproducibility package.

## Published Zenodo v1

The prepared bundle has been deposited and published as:

```text
DOI:      10.5281/zenodo.22295791
record:   https://zenodo.org/records/22295791
version:  v1
bundle:   azami_ch1_v2_reproduction_inputs_2026-09-04.zip
size:     56,942,044 bytes
SHA-256:  50ec15b1280d4660839ca4bf0d55c970a5f49b4d4feaabb7a073b73500253677
```

The bundle contains the exact four minimum binary inputs recorded in `reproducibility/public_release_manifest.json`, plus a SHA ledger, public-bundle manifest, release README and Zenodo metadata. Each embedded binary was re-read from the finished bundle and verified against its frozen SHA-256 before publication.

The immutable code ref paired with this bundle is:

```text
584af97b050d15701f26ce1facea212d5b648d4d
```

## Publication status

The DOI publication and anonymous post-publication redownload gate are complete.

Current state:

```text
ready_for_independent_third_party_numerical_reproduction
```

This state means the public package can be obtained and byte-verified by an independent third party. It does not assert that every future reproduction run has already been executed by every third party; an actual reproduction must still satisfy the canonical numerical and figure criteria below.

## Independent post-publication verification

On 2026-09-04, a clean GitHub-hosted Ubuntu runner downloaded the public Zenodo file URL with `curl` without Zenodo author credentials or an access token. Verification run:

```text
https://github.com/zuizui0223/azami/actions/runs/33871722778
```

Observed public-download checks:

```text
anonymous download: PASS
bundle size:        56,942,044 bytes
bundle SHA-256:     50ec15b1280d4660839ca4bf0d55c970a5f49b4d4feaabb7a073b73500253677
embedded inputs:    PASS 4/4
```

The four extracted input SHA-256 values matched `reproducibility/public_release_manifest.json` exactly:

| File | SHA-256 |
|---|---|
| `artifact-9612943217-continuous.zip` | `101e996b638996a0c5ae79d358bf51293c3585f0e84c4a961b91dcbedf96211e` |
| `artifact-9632715852-multilevel.zip` | `51e7a26b5bd09e030b67b9342586699abaaf46e630f45b6bb4ee7bfc9152ced6` |
| `artifact-8227254443-historical.zip` | `499061e7a49f9455cf8c367fe26e313b7e0e33b2280d2354717e61a90ea8c6bc` |
| `artifact-8983877726-spatial.zip` | `151161f926dbc92d35832a0bd71622e4f1f3018c9047a7ba8b4ad94926273dca` |

For a fresh independent reproduction, the reader should still:

1. open `https://doi.org/10.5281/zenodo.22295791`;
2. download `azami_ch1_v2_reproduction_inputs_2026-09-04.zip` without author credentials;
3. require bundle size `56,942,044` bytes and SHA-256 `50ec15b1280d4660839ca4bf0d55c970a5f49b4d4feaabb7a073b73500253677`;
4. extract it;
5. clone `https://github.com/zuizui0223/azami` and checkout `584af97b050d15701f26ce1facea212d5b648d4d`;
6. run `python reproducibility/verify_local_materials.py --archive-dir <extracted-directory> --scope numerical`;
7. require all four embedded input SHA-256 values to match `reproducibility/public_release_manifest.json`;
8. follow `reproducibility/README.md` through the canonical rebuild;
9. require `PASS 24/24`, `PASS 7/7`, `PASS 11/11`, and figure build `status: ok` with 5 main and 7 supporting figures.

## Publication gate completion

The anonymous credential-free gate was satisfied on 2026-09-04: public download succeeded, the exact bundle size and SHA-256 matched, and all four embedded input SHA-256 values matched. The release manifest, material-availability ledger and recovery inventory therefore mark the package ready for independent third-party numerical reproduction. GitHub issue #88 is the completion record for this gate; repository integrity must remain green on `main`.

The five frozen CHELSA v2.1 process rasters may optionally be added to a future Zenodo version for a fully offline rerun. If omitted, their exact public URLs remain frozen in `analysis/ch1/chelsa_process_environment_sources.json`, and the reconstructed nine-predictor environment is still required to match SHA-256 `e242aa7ce69d12b11937c1335e84b9638799c50b42ef36b95725e77190df98e7`.
