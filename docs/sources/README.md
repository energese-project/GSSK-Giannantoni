# Sources: text, page renders and probes

`PLAN.md` cites equations as `[YY Eq n]`. This directory holds what those citations were checked
against: machine-readable text, renders of the equation pages, and the probe scripts that decide each
erratum. The PDFs themselves live in `docs/*.pdf` and are gitignored.

## Citation keys

| Key | Work | Licence | Text here |
|---|---|---|---|
| `[02]` | Giannantoni, C. (2002) *The Maximum Em-Power Principle as the Basis for Thermodynamics of Quality*. SGE, Padova. ISBN 88-86281-76-5 | All rights reserved | `local/` (OCR; scanned two book pages per PDF sheet, book page ≈ 2·sheet − 24) |
| `[06]` | Giannantoni, C. (2006) Mathematics for generative processes. *J. Comput. Appl. Math.* 189, 324–340 | Elsevier | `local/` |
| `[06b]` | Giannantoni, C. (2006) Emergy Analysis as the First Ordinal Theory of Complex Systems. *4th Emergy Conf.*, 15.1–15.14 | Proceedings | `local/` |
| `[09]` | Giannantoni, C. & Zoli, M. (2009) The Derivative "Drift" in Complex Systems Dynamics. *ECOS 2009* | Proceedings | `local/` |
| `[10]` | Giannantoni, C. (2010) The Maximum Ordinality Principle. A Harmonious Dissonance. *6th Emergy Conf.*, 55–72 (two printings in `docs/`, same text) | Proceedings | `local/` |
| `[10p]` | Giannantoni, C. (2010) Protein Folding, Molecular Docking, Drug Design. *BIOINFORMATICS 2010* | Proceedings | `local/` |
| `[22]` | Giannantoni, C. (2022) The MOP and Its Formal Language, IDC. *J. Appl. Math. Phys.* 10, 2649–2689. doi:10.4236/jamp.2022.109178 | **CC BY 4.0** | `cc-by/` |
| `[23]` | Giannantoni, C. (2023) Generativity of Self-Organizing Processes and Their Maximum Ordinality. *J. Appl. Math. Phys.* 11. doi:10.4236/jamp.2023.1110206 | **CC BY 4.0** | `cc-by/` |

## What is committed and what is not

- **`cc-by/`** is committed. It holds text extracted from `[22]` and `[23]` (`pdftotext -layout`) and
  renders of `[23]` pp. 3183, 3186, 3189 and 3190, redistributed under CC BY 4.0 with the attribution
  above.
- **`probes/`** is committed. These are this project's own scripts, and each decides one erratum or
  resolution in `PLAN.md` §4–5. Run each with `python3 -I docs/sources/probes/<name>.py`.
- **`local/`** is gitignored. It holds text, OCR and page renders of works that are not openly
  licensed, kept beside the PDFs they came from, for the same reason those PDFs are not in a public
  repository. To regenerate it:

  ```sh
  pdftotext -layout docs/<paper>.pdf docs/sources/local/text/<name>.txt
  pdftoppm -r 200 -gray -png docs/2002_giannantoni_*.pdf /tmp/g2002/p      # book: no text layer
  for f in /tmp/g2002/p-*.png; do tesseract "$f" "${f%.png}" -l eng --psm 6; done
  ```

  The OCR is reliable for prose and unreliable for equations. Every equation `PLAN.md` relies on was
  read from a page render, not from OCR or a text layer.
