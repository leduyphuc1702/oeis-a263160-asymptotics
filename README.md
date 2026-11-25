# Exact Diagonal Asymptotics for Multidimensional Lattice Paths (OEIS A263160 & A262810)

This repository contains the LaTeX manuscript and independent numerical verification code proving two diagonal asymptotic formulas in the On-Line Encyclopedia of Integer Sequences (OEIS):

1. **OEIS [A263160](https://oeis.org/A263160)** (first-zero lattice paths starting at $(n,\dots,n)\in\mathbb{Z}_{>0}^n$):
   $$a(n) \sim C\,\frac{n^{n(n-3/2)+1}}{\pi^{n/2}(\log 2)^{n(n-1)}}, \qquad C = \frac{2\sqrt{2\pi}}{\log 2}\exp\!\left(\frac{5}{12}-\frac{(\log 2)^2}{24}\right) = 10.753662835613289269\ldots$$
   This proves Vaclav Kotesovec's November 13, 2026 asymptotic assertion and determines the exact closed form of the constant $C$.

2. **OEIS [A262810](https://oeis.org/A262810)** (origin-reaching lattice paths from $(n,\dots,n)$ to $(0,\dots,0)$, equivalently many-to-many alignments of $n$ sequences of length $n$):
   $$b(n) \sim \frac{n^{n^2 - n/2 + 1}}{e^{1/12}\,2^{n + (\log 2)/24}\,\pi^{(n-1)/2}\,(\log 2)^{n^2+1}},$$
   proving Vaclav Kotesovec's March 23, 2016 asymptotic formula within the same unified framework.

## Contents

- `paper.tex` — Self-contained LaTeX manuscript (standard packages only: `amsmath`, `amssymb`, `amsthm`, `booktabs`, `geometry`, `hyperref`).
- `verify.py` — Zero-dependency Python 3 verification script (direct combinatorial DP for $n\le 4$, 50-digit evaluation of $C$, comparison with OEIS b-files, and cancellation-free large-$n$ evaluation up to $n=1000$).

## Reproducing the Numerical Verification

```bash
python3 verify.py
```

To compare directly against the official OEIS b-files:

```bash
curl -sL https://oeis.org/A263160/b263160.txt -o /tmp/b263160.txt
curl -sL https://oeis.org/A262810/b262810.txt -o /tmp/b262810.txt
python3 verify.py --b263160 /tmp/b263160.txt --b262810 /tmp/b262810.txt
```

## Compiling the Manuscript

Compile `paper.tex` on Overleaf (upload `paper.tex` directly) or locally with any standard TeX engine:

```bash
pdflatex paper.tex && pdflatex paper.tex
# or:
tectonic paper.tex
```

## License

Manuscript text: CC BY 4.0. Verification code (`verify.py`): MIT License.
